#include "function.cc"
#include "TaAccumulator.cc"
#include "lib/TaRunInfo_v2.cc"
#include "LoadRunInfoMap.C"

vector<STAT> solve(vector<STAT*> fVar_det, 
		   vector<STAT*> fCovar_det,
		   vector<STAT*> fCovar_bpm);

void regression(Int_t slug=94){
  map<Int_t,TaRunInfo> fRunInfoMap = LoadRunInfoMap();
  TFile *cov_file = TFile::Open(Form("rootfiles/slug%d_covv.root",slug));
  TTree *cov_tree = (TTree*)cov_file->Get("covv");
  TString output_filename;
  output_filename = Form("rootfiles/slug%d_reg.root",slug);
  TFile *out_file = TFile::Open(output_filename,"RECREATE");
  TString tree_name = "reg";
  TTree *cor_tree = new TTree(tree_name,"correction tree");

  vector<TString> det_list={"asym_usl","asym_usr","asym_us_avg"};
  vector<TString> bpm_list;
  vector<TString> bpm_set1={"diff_bpm4aX","diff_bpm4eX","diff_bpm4aY","diff_bpm4eY",
			    "diff_bpm11X","diff_bpm12X"};
  vector<TString> bpm_set2={"diff_bpm4aX","diff_bpm4eX","diff_bpm4aY","diff_bpm4eY",
			    "diff_bpm12X"};
  if(slug>=4)
    bpm_list = bpm_set1;
  else{
    bpm_list = bpm_set2;
  }
  int ndet = det_list.size();
  int nbpm = bpm_list.size();
  vector<STAT*> fVar_det(ndet);
  vector<STAT*> fCovar_det(ndet*nbpm);
  vector<STAT*> fCovar_bpm(nbpm*nbpm);
  for(int idet=0;idet<ndet;idet++){
    cov_tree->SetBranchAddress(det_list[idet],&fVar_det[idet]);
    for(int ibpm=0;ibpm<nbpm;ibpm++){
      TString det_base = get_base(det_list[idet]);
      TString bpm_base = get_base(bpm_list[ibpm]);
      TString branch_name = Form("%s_%s",det_base.Data(),bpm_base.Data());
      cov_tree->SetBranchAddress(branch_name,&fCovar_det[idet*nbpm+ibpm]);
    }
  }
  for(int ibpm=0;ibpm<nbpm;ibpm++){
    TString bpm_base1 = get_base(bpm_list[ibpm]);
    for(int jbpm=0;jbpm<nbpm;jbpm++){
      TString bpm_base2 = get_base(bpm_list[jbpm]);
      TString branch_name = Form("%s_%s",bpm_base1.Data(),bpm_base2.Data());
      cov_tree->SetBranchAddress(branch_name,&fCovar_bpm[ibpm*nbpm+jbpm]);
    }
  }
  
  Int_t run_id;
  Double_t fSamples;
  cov_tree->SetBranchAddress("run",&run_id);
  cov_tree->SetBranchAddress("num_samples",&fSamples);
  Int_t arm_flag,ihwp,wien,sign,good,fSlug;
  cor_tree->Branch("arm",&arm_flag);
  cor_tree->Branch("ihwp",&ihwp);
  cor_tree->Branch("wien",&wien);
  cor_tree->Branch("sign",&sign);
  cor_tree->Branch("good",&good);
  cor_tree->Branch("slug",&fSlug);
  vector<STAT> fRegCorrection(ndet);
  for(int idet=0;idet<ndet;idet++)
    cor_tree->Branch("reg_"+get_base(det_list[idet]),&fRegCorrection[idet]);

  typedef struct{Double_t ppm,ppb,mm,um,nm;} UNIT;
  UNIT parity_unit;
  parity_unit.ppm = 1e-6;
  parity_unit.ppb = 1e-9;
  parity_unit.mm = 1;
  parity_unit.um = 1e-3;
  parity_unit.nm = 1e-6;
  cor_tree->Branch("unit",&parity_unit,"ppm/D:ppb:mm:um:nm");
  //FIXME raw detector and raw bpm tree for residual correlation check
  cor_tree->Branch("run",&run_id);

  Int_t nrun  = cov_tree->GetEntries();
  for(int ievt=0;ievt<nrun;ievt++){
    cov_tree->GetEntry(ievt);
    vector<STAT> fRegRes = solve(fVar_det,fCovar_det,fCovar_bpm);
    for(int idet=0;idet<ndet;idet++)
      fRegCorrection[idet] = fRegRes[idet];
    
    arm_flag=-1;
    sign = 0.0;
    good = -1;
    if(fRunInfoMap.find(run_id)==fRunInfoMap.end())
      continue;
    TaRunInfo myRunInfo = fRunInfoMap[run_id];
    arm_flag = myRunInfo.GetArmFlag();
    sign = myRunInfo.GetSign();
    fSlug = myRunInfo.GetSlugNumber(); 
    if(myRunInfo.GetIHWPStatus()=="IN") 
      ihwp=1;
    else
      ihwp=-1;

    if(myRunInfo.GetWienMode()=="FLIP-RIGHT") 
      wien=1;
    else
      wien=-1;

    if(myRunInfo.GetRunFlag()=="Good") 
      good=1;
    else
      good=0;
    cout << run_id << endl;
    cor_tree->Fill();
  }
  out_file->cd();
  cor_tree->Write();
  cout << " -- Closing " << endl;
  out_file->Close();
}

vector<STAT> solve(vector<STAT*> fVar_det, 
		   vector<STAT*> fCovar_det,
		   vector<STAT*> fCovar_bpm){

  Int_t ndet = fVar_det.size();
  Int_t nbpm = fCovar_det.size() / ndet;
  Bool_t kUsedCombined=kFALSE;
  if(nbpm==6){
    nbpm = 5; // forced to 11X+0.4*12X
    kUsedCombined=kTRUE;
    cout << " -- Using Combo BPM for regression " << endl;
  }
  vector<STAT> fRegDet(ndet);  
  TMatrixD m_cov(nbpm,nbpm);
  TMatrixD d_cov(ndet,nbpm);
  TMatrixD d_cov_trans(nbpm,ndet);

  if(kUsedCombined){
    for(int ibpm=0;ibpm<nbpm;ibpm++){
      for(int jbpm=0;jbpm<nbpm;jbpm++){
	m_cov[ibpm][jbpm] = fCovar_bpm[ibpm*6+jbpm]->m2;
	if(ibpm==4 && ibpm==jbpm)
	  m_cov[ibpm][jbpm] += 2*0.4*(fCovar_bpm[4*6+5]->m2) + 0.16*(fCovar_bpm[5*6+5]->m2);
	else if (ibpm==4)
	  m_cov[ibpm][jbpm]+=0.4*(fCovar_bpm[5*6+jbpm]->m2);
	else if (jbpm==4)
	  m_cov[ibpm][jbpm]+=0.4*(fCovar_bpm[ibpm*6+5]->m2);

      }
    }
    for(int idet=0;idet<ndet;idet++){
      for(int ibpm=0;ibpm<nbpm;ibpm++){
	d_cov[idet][ibpm] = fCovar_det[idet*6+ibpm]->m2;
	if(ibpm==4)
	  d_cov[idet][ibpm]+= 0.4*(fCovar_det[idet*6+5]->m2);

	d_cov_trans[ibpm][idet] = d_cov[idet][ibpm];
      }
    }
  }else{
    for(int ibpm=0;ibpm<nbpm;ibpm++)
      for(int jbpm=0;jbpm<nbpm;jbpm++)
	m_cov[ibpm][jbpm] = fCovar_bpm[ibpm*nbpm+jbpm]->m2;

    for(int idet=0;idet<ndet;idet++){
      for(int ibpm=0;ibpm<nbpm;ibpm++){
	d_cov[idet][ibpm] = fCovar_det[idet*nbpm+ibpm]->m2;
	d_cov_trans[ibpm][idet] = d_cov[idet][ibpm];
      }
    }
  }

  m_cov.Invert();
  TMatrixD sol = d_cov*m_cov;
  // sol.Print();
  for(int idet=0;idet<ndet;idet++){
    fRegDet[idet].mean = fVar_det[idet]->mean;
    fRegDet[idet].m2= fVar_det[idet]->m2;
    if(kUsedCombined){
      for(int ibpm=0;ibpm<nbpm;ibpm++){
	fRegDet[idet].mean -= sol[idet][ibpm]*fCovar_bpm[ibpm*6+ibpm]->mean;
	fRegDet[idet].m2 -= 2*sol[idet][ibpm]*fCovar_det[idet*6+ibpm]->m2;
	for(int jbpm=0;jbpm<nbpm;jbpm++)
	  fRegDet[idet].m2 += sol[idet][ibpm]*sol[idet][jbpm]*fCovar_bpm[ibpm*6+jbpm]->m2;
      }
      fRegDet[idet].mean -= 0.4*sol[idet][4]*fCovar_bpm[35]->mean;
      fRegDet[idet].m2 -= 2*0.4*sol[idet][4]*fCovar_det[idet*6+5]->m2;
      for(int jbpm=0;jbpm<nbpm;jbpm++)
	fRegDet[idet].m2 += 2*0.4*sol[idet][jbpm]*sol[idet][4]*fCovar_bpm[jbpm*6+5]->m2;
      fRegDet[idet].m2 += 0.16*sol[idet][4]*sol[idet][4]*fCovar_bpm[35]->m2;

    }else{ // else if not use combo BPM
      for(int ibpm=0;ibpm<nbpm;ibpm++){
	fRegDet[idet].mean -= sol[idet][ibpm]*fCovar_bpm[ibpm*nbpm+ibpm]->mean;
	fRegDet[idet].m2 -= 2*sol[idet][ibpm]*fCovar_det[idet*nbpm+ibpm]->m2;
	for(int jbpm=0;jbpm<nbpm;jbpm++)
	  fRegDet[idet].m2 += sol[idet][ibpm]*sol[idet][jbpm]*fCovar_bpm[ibpm*nbpm+jbpm]->m2;
      }
    }

    fRegDet[idet].rms= sqrt(fRegDet[idet].m2/fVar_det[idet]->num_samples);
  }
  return fRegDet;
}

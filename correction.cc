#include "function.cc"
#include "TaAccumulator.cc"

void correction(Int_t slug=1){
  TFile *cov_file = TFile::Open(Form("rootfiles/slug%d_covv.root",slug));
  TTree *cov_tree = (TTree*)cov_file->Get("covv");
  TFile *slope_file = TFile::Open(Form("slopes/slug%d_dit_slope.root",slug));
  TTree *slope_tree = (TTree*)slope_file->Get("slope");

  TFile *out_file = TFile::Open(Form("rootfiles/slug%d_correct.root",slug),"RECREATE");
  TTree *cor_tree = new TTree("cor","correction tree");

  
  vector<TString> det_list={"asym_usl","asym_usr","asym_us_avg"};
  vector<TString> iv_list = {"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm11X12X"};
  vector<TString> iv_old = {"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm12X"};
  vector<TString> bpm_list;
  vector<TString> bpm_set1={"diff_bpm4aX","diff_bpm4eX","diff_bpm4aY","diff_bpm4eY",
			    "diff_bpm11X","diff_bpm12X"};
  vector<TString> bpm_set2={"diff_bpm4aX","diff_bpm4eX","diff_bpm4aY","diff_bpm4eY",
			    "diff_bpm12X"};
  if(slug>=4)
    bpm_list = bpm_set1;
  else{
    bpm_list = bpm_set2;
    iv_list = iv_old;
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
  
  map<Int_t, vector<Double_t> > fSlopeMap;
  fSlopeMap = LoadSlopeMap(slope_tree, det_list, iv_list);
  
  Int_t run_id;
  cov_tree->SetBranchAddress("run",&run_id);
  Double_t fSamples;
  cov_tree->SetBranchAddress("num_samples",&fSamples);
  vector<STAT> fDitCorrection(ndet);
  for(int idet=0;idet<ndet;idet++)
    cor_tree->Branch("dit_"+get_base(det_list[idet]),&fDitCorrection[idet]);

  typedef struct{Double_t ppm,ppb,mm,um,nm;} UNIT;
  UNIT parity_unit;
  parity_unit.ppm = 1e-6;
  parity_unit.ppb = 1e-9;
  parity_unit.mm = 1;
  parity_unit.um = 1e-3;
  parity_unit.nm = 1e-6;
  cor_tree->Branch("unit",&parity_unit,"ppm/D:ppb:mm:um:nm");
  //FIXME raw detector and raw bpm tree for residual correlation check

  Int_t nrun  = cov_tree->GetEntries();
  for(int ievt=0;ievt<nrun;ievt++){
    cov_tree->GetEntry(ievt);
    if(fSlopeMap.find(run_id)==fSlopeMap.end()){
      cout << " -- run " << run_id
	   << " slope not found " << endl;
      continue;
    }
    vector<Double_t> my_slope = fSlopeMap[run_id];
    for(int idet=0;idet<ndet;idet++){
      fDitCorrection[idet].num_samples = fVar_det[idet]->num_samples;
      fDitCorrection[idet].m2 = fVar_det[idet]->m2;
      fDitCorrection[idet].mean = fVar_det[idet]->mean;
      for(int ibpm=0;ibpm<nbpm;ibpm++){
      	fDitCorrection[idet].mean -=my_slope[idet*nbpm+ibpm]*(fCovar_bpm[ibpm*nbpm+ibpm]->mean);

      	fDitCorrection[idet].m2 += 2*my_slope[idet*nbpm+ibpm]*(fCovar_det[idet*nbpm+ibpm]->m2);
	
      	for(int jbpm=0;jbpm<nbpm;jbpm++){
      	  fDitCorrection[idet].m2 += my_slope[idet*nbpm+ibpm]*my_slope[idet*nbpm+jbpm]*(fCovar_bpm[ibpm*nbpm+jbpm]->m2);
      	}
      }

      fDitCorrection[idet].rms = sqrt(fDitCorrection[idet].m2/fDitCorrection[idet].num_samples);
      fDitCorrection[idet].err = fDitCorrection[idet].rms /sqrt(fDitCorrection[idet].num_samples);
    }
    
    cor_tree->Fill();
  }
  out_file->cd();
  cor_tree->Write();
  out_file->Close();
}



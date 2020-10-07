/*
   author: Tao Ye
   Last update: Aug 2019
 */
#include "utilities.cc"
void CalculateResidual(Int_t slug,Int_t kSwitch);
void CalculateResidual(){
  for(int i=1; i<=94;i++){
    if(i==7)
      continue;
    CalculateResidual(i,1);
    CalculateResidual(i,2);
  }
}

void CalculateResidual(Int_t slug,Int_t kSwitch){
  // kSwitch: 1-Lagrange, 2-Regression
  TString tag, full_tag;
  if(kSwitch==1){
    tag="lagr";
    full_tag = "lagrange";
  }else if(kSwitch==2){
    tag="reg";
    full_tag = "regression";
  }

  TFile *output = TFile::Open(Form("rootfiles/slug%d_%s_residual.root",
				   slug,full_tag.Data()),"RECREATE");
  TTree *res_tree  = new TTree("res"," Residual sensitivities");
  
  Int_t nDet  = 4;
  Int_t nCoil  = 7;
  Int_t nBPM = 12;
  if(slug<=2)
    nBPM = 10;
  vector<TString> fDet_name = { "usl","usr","us_avg","us_dd"};
  vector<TString> fBPM_name;
  vector<TString> fBPM_name1 = {"11X","11Y","12X","12Y","16X","16Y",
				"1X","1Y","4aX","4aY","4eX","4eY"};
  vector<TString> fBPM_name2 = {"8X","8Y","12X","12Y",
				"1X","1Y","4aX","4aY","4eX","4eY"};
  if(slug<=2)
    fBPM_name = fBPM_name2;
  else
    fBPM_name = fBPM_name1;
  
  vector<Double_t > fDet_res(nDet*nCoil,0);  // residual(corrected) sensitivities
  vector<Double_t > fDet_sens(nDet*nCoil,0); // raw sensitivities
  vector<Int_t > fDet_flag(nDet*nCoil,0);
  Int_t arm_flag;
  Int_t run_flag;
  Int_t kCycle, kRun, kBurstCounter,kSlug;
  kSlug = slug;
  res_tree->Branch("arm_flag",&arm_flag);
  res_tree->Branch("kGood",&run_flag);
  res_tree->Branch("bmwcycnum",&kCycle);
  res_tree->Branch("run",&kRun);
  res_tree->Branch("slug",&kSlug);
  res_tree->Branch("BurstCounter",&kBurstCounter);
  for(int idet=0;idet<nDet;idet++){
    for(int icoil=0;icoil<nCoil;icoil++){
      res_tree->Branch(Form("%s_coil%d_res",fDet_name[idet].Data(),icoil+1), 
		       &fDet_res[idet*nCoil+icoil]);
      res_tree->Branch(Form("%s_coil%d_flag",fDet_name[idet].Data(),icoil+1), 
		       &fDet_flag[idet*nCoil+icoil]);
      res_tree->Branch(Form("%s_coil%d_sens",fDet_name[idet].Data(),icoil+1), 
		       &fDet_sens[idet*nCoil+icoil]);

    }
  }

  TFile* sens_file = TFile::Open(Form("ditcoeffs/slug%d.root",slug));
  /* Note : the first two slugs are from respin2, because of channel map */
  TTree* sens_tree = (TTree*)sens_file->Get("sens");
  vector<Double_t> fBPM_sens(nBPM*nCoil,0);
  vector<Double_t> fBPM_err(nBPM*nCoil,0);
  vector<Double_t > fDet_err(nDet*nCoil,0);
  Double_t fRun, fCycle;
  sens_tree->SetBranchAddress("run",&fRun);
  sens_tree->SetBranchAddress("cycID",&fCycle);
  
  for(int idet=0;idet<2;idet++){ // Only Load the first two
    for(int icoil=0;icoil<nCoil;icoil++){
      sens_tree->SetBranchAddress(Form("%s_coil%d",fDet_name[idet].Data(),icoil+1),
				  &fDet_sens[idet*nCoil+icoil]);
      sens_tree->SetBranchAddress(Form("%s_coil%d_err",fDet_name[idet].Data(),icoil+1),
				  &fDet_err[idet*nCoil+icoil]);
    }
  }

  for(int ibpm=0;ibpm<nBPM;ibpm++){
    for(int icoil=0;icoil<nCoil;icoil++){
      sens_tree->SetBranchAddress(Form("bpm%s_coil%d",fBPM_name[ibpm].Data(),icoil+1),
				  &fBPM_sens[ibpm*nCoil+icoil]);
      sens_tree->SetBranchAddress(Form("bpm%s_coil%d_err",fBPM_name[ibpm].Data(),icoil+1),
				  &fBPM_err[ibpm*nCoil+icoil]);
    }
  }

  map<Int_t, Int_t> fArmFlagMap = LoadArmFlag(slug);
  map<Int_t, Int_t> fRunFlagMap = LoadRunFlag(slug);
  map< pair<Int_t,Int_t>, vector<Double_t> > fSlopeMap;
  if(kSwitch==1)
    fSlopeMap= GetLagrangeSlope(slug);
  else if(kSwitch==2)
    fSlopeMap= GetRegressionSlope(slug);
  map< pair<Int_t,Int_t>, Int_t > fC2BMap;
  map< Int_t, vector<Int_t> > fBadCycleMap = LoadBadCycleList();
  Int_t nEntries = sens_tree->GetEntries();
  Int_t last_run = -1;
  for(int ievt=0;ievt<nEntries;ievt++){
    sens_tree->GetEntry(ievt);
    kRun = (Int_t) fRun;
    kCycle = (Int_t) fCycle;
    // 1. Load Arm Flags from run info
    if( fArmFlagMap.find(kRun) == fArmFlagMap.end())  
      continue; // if a run has not found in run info,skip.
    if(kRun!=last_run){
      last_run=kRun;
      fC2BMap =  LoadCyc2BurstMap(slug,kRun);
      kBurstCounter = 0; // for weird slug 2
    }
    arm_flag = fArmFlagMap[kRun];
    run_flag = fRunFlagMap[kRun];
    for(int icoil=0;icoil<nCoil;icoil++){
      // 2. Check if it is in bad cycle list 
      Bool_t kBadCyceList = InBadCycleList(fBadCycleMap,kCycle,icoil+1);
      Bool_t fBPM_flag = 1;
      for(int ibpm=0;ibpm<nBPM;ibpm++){
	if( fBPM_err[ibpm*nCoil+icoil]<=0)
	  fBPM_flag = 0 ; // if any BPM fails in this coil, then failed all detectors
      }

      // 3. Find the burst-counter and map out Lagrange/Regression slopes
      /* 
	 What if a bmod coil has no mapping for burst counter ? 
	 This is a rare case. It will use previous burst counter. 
      */
      if( fC2BMap.find( make_pair(kCycle,icoil+1)) != fC2BMap.end())
	  kBurstCounter = fC2BMap[ make_pair( kCycle, icoil+1)  ];
      vector<Double_t> fSlopes = fSlopeMap[ make_pair(kRun,kBurstCounter) ] ;
      // 4. Combined detector sensitivities and error flag;
      fDet_sens[2*nCoil + icoil] = (fDet_sens[0*nCoil+icoil] + fDet_sens[1*nCoil+icoil])*0.5;
      fDet_sens[3*nCoil + icoil] = (fDet_sens[0*nCoil+icoil] - fDet_sens[1*nCoil+icoil])*0.5;

      if(fDet_err[0*nCoil+icoil] > 0 && fDet_err[1*nCoil+icoil]> 0 ){
	fDet_err[2*nCoil+icoil] = 1;
	fDet_err[3*nCoil+icoil] = 1; // so that it is good;
      }else{  // any of it fails
	fDet_err[2*nCoil+icoil] = -1;
	fDet_err[3*nCoil+icoil] = -1; // so that it is bad;
      }
      
      for(int idet=0;idet<nDet;idet++){
	// 5. Check if it is a good dithering coil  
	if(fDet_err[idet*nCoil+icoil] <= 0 || (fBPM_flag == 0) || kBadCyceList)
	  fDet_flag[idet*nCoil+icoil] = 0; // 0 : bad
	else
	  fDet_flag[idet*nCoil+icoil] = 1; // 1 : good

	// 6. Calculate residual sensitivities 
	fDet_res[idet*nCoil+icoil] = fDet_sens[idet*nCoil+icoil];
	for(int ibpm=0;ibpm<nBPM;ibpm++)
	  fDet_res[idet*nCoil+icoil]  -= fSlopes[idet*nBPM+ibpm] * fBPM_sens[ibpm*nCoil+icoil];
	
      } // end of loop over detectors
    } // end of loop over dithering coils
    res_tree->Fill();    
  } // end of loop over dithering cycles
  sens_file->Close();
  output->cd();
  res_tree->Write();
  output->Close();
			      
}

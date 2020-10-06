void Extractor(Int_t run_number){
  TString normal_cut = "ErrorFlag==0";
  TString combine_cut = "(ErrorFlag&0xfbfe6fff)==0";
  TString bmw_only_cut = "ErrorFlag==0x4019000";
  /*
    kBModFFBErrorFlag + kBModErrorFlag + kEventCutMode3 + kGlobalCut = 0x4019000
    0xFFFFFFFF - 0x4019000  = 0xfbfe6fff;
  */
  vector<TString> device_list={"asym_bcm_target","yield_bcm_target",
			       "yield_bpm4eX","yield_bpm4eY",
			       "yield_bpm4aX","yield_bpm4aY",
			       "yield_bpm1X","yield_bpm1Y",
			       "yield_bpm8X","yield_bpm8Y",
			       "yield_bpm11X","yield_bpm11Y",
			       "yield_bpm12X","yield_bpm12Y",
			       "yield_bpm16X","yield_bpm16Y",
			       "diff_bpm4eX","diff_bpm4eY",
			       "diff_bpm4aX","diff_bpm4aY",
			       "diff_bpm1X","diff_bpm1Y",
			       "diff_bpm8X","diff_bpm8Y",
			       "diff_bpm11X","diff_bpm11Y",
			       "diff_bpm12X","diff_bpm12Y",
			       "diff_bpm16X","diff_bpm16Y",
			       "asym_us_avg","asym_us_dd","asym_usl","asym_usr",
			       "yield_usl","yield_usr",
			       "lagr_asym_us_avg","lagr_asym_us_dd",
			       "lagr_asym_usl","lagr_asym_usr",
			       "reg_asym_us_avg","reg_asym_us_dd",
			       "reg_asym_usl","reg_asym_usr",
			       "lagr_asym_us_avg","lagr_asym_us_dd",
			       "lagr_asym_usl","lagr_asym_usr",
			       "dit_asym_us_avg","dit_asym_us_dd",
			       "dit_asym_usl","dit_asym_usr" };
  for(int i=0;i<12;i++)
    device_list.push_back(Form("diff_evMon%d",i) );
  
  Int_t nDev = device_list.size();

  TFile *japan_file = TFile::Open(Form("japanOutput/prexPrompt_pass2_%d.000.root",run_number));
  TTree *mul_tree = (TTree*)japan_file->Get("mul");
  mul_tree->AddFriend( (TTree*)japan_file->Get("mulc") );
  mul_tree->AddFriend( (TTree*)japan_file->Get("mulc_dit") );
  mul_tree->AddFriend( (TTree*)japan_file->Get("mulc_dit_combo") );

  TFile *lagr_file = TFile::Open(Form("LagrangeOutput/prexRespin2_lagrange_eigen_%d.000.root",run_number));
  TTree *lagr_tree = (TTree*)lagr_file->Get("lagrall");
  TTree *reg_tree = (TTree*)lagr_file->Get("regall");
  
  mul_tree->AddFriend(lagr_tree);
  mul_tree->AddFriend(reg_tree);
  
  TString output_filename = Form("output/prexRespin2_bmw_extra_%d.root",run_number);
  TFile* output = TFile::Open(output_filename,"RECREATE");
  TTree *mini_tree = new TTree("mini","Normal");
  TTree *mini_plus_tree = new TTree("mini_plus","Normal + BMW");
  TTree *mini_bmw_tree =  new TTree("mini_bmw","BMW only");
  Int_t nCoil = 7;
  vector<TTree*> mini_bmw_tree_by_coil(nCoil);
  for(int icoil=0;icoil<nCoil;icoil++){
    mini_bmw_tree_by_coil[icoil] = new TTree(Form("mini_bmw_%d",icoil+1),
					     "BMW Only By Coil");
  }
  Int_t mini_counter=0;
  
  mini_tree->Branch("run",&run_number);
  mini_tree->Branch("minirun",&mini_counter);
  mini_plus_tree->Branch("run",&run_number);
  mini_plus_tree->Branch("minirun",&mini_counter);
  mini_bmw_tree->Branch("run",&run_number);
  mini_bmw_tree->Branch("minirun",&mini_counter);
  
  typedef struct {Double_t mean, error,rms,num_samples;} STAT;
  STAT fStat_zero;
  fStat_zero.mean=0.0;
  fStat_zero.error=-1;
  fStat_zero.rms=0.0;
  fStat_zero.num_samples=0.0;
  vector<STAT> fStat_mini(nDev);
  vector<STAT> fStat_plus(nDev);
  vector<STAT> fStat_bmw(nDev);
  vector<STAT> fStat_bmw_by_coil(nDev*nCoil);

  for(int i=0;i<nDev;i++){
    mini_tree->Branch(device_list[i], &fStat_mini[i],"mean/D:err:rms:num_samples");
    mini_plus_tree->Branch(device_list[i], &fStat_plus[i],"mean/D:err:rms:num_samples");
    mini_bmw_tree->Branch(device_list[i], &fStat_bmw[i],"mean/D:err:rms:num_samples");
    for(int icoil=0;icoil<nCoil;icoil++){
      mini_bmw_tree_by_coil[icoil]->Branch(device_list[i], 
					   &fStat_bmw_by_coil[i*nCoil+icoil],
					   "mean/D:err:rms:num_samples");
    }
  }
  cout << " -- Run " << run_number << endl;
  while( mul_tree->GetEntries(combine_cut + Form("&& BurstCounter==%d",mini_counter)) !=0){
    cout << " -- Burst " << mini_counter << endl;
    TString burst_cut = Form(" && BurstCounter==%d",mini_counter);
    for(int idev=0;idev<nDev;idev++){
      if( mul_tree->GetBranch(device_list[idev])==NULL){
	cout << " ** skip " << device_list[idev] << endl;
	fStat_mini[idev] = fStat_zero;
	fStat_plus[idev] = fStat_zero;
	fStat_bmw[idev] = fStat_zero;
	continue;
      }
      TH1D *h1dptr;
      int npt = mul_tree->Draw(device_list[idev],normal_cut+burst_cut,"goff");
      if(npt!=0){
	h1dptr= (TH1D*)gDirectory->FindObject("htemp");
	fStat_mini[idev].mean = h1dptr->GetMean();
	fStat_mini[idev].error = h1dptr->GetMeanError();
	fStat_mini[idev].rms = h1dptr->GetRMS();
	fStat_mini[idev].num_samples = h1dptr->GetEntries();
      }else{
	cout << " ** Empty channel " << device_list[idev] << " with normal cut. "  << endl;
	fStat_mini[idev] = fStat_zero;
      }

      npt = mul_tree->Draw(device_list[idev],combine_cut+burst_cut,"goff");
      if(npt!=0){
	h1dptr= (TH1D*)gDirectory->FindObject("htemp");
	fStat_plus[idev].mean = h1dptr->GetMean();
	fStat_plus[idev].error = h1dptr->GetMeanError();
	fStat_plus[idev].rms = h1dptr->GetRMS();
	fStat_plus[idev].num_samples = h1dptr->GetEntries();
      }else{
	cout << " ** Empty channel " << device_list[idev] << " with combined cut. "  << endl;
	fStat_plus[idev] = fStat_zero;
      }
      npt = mul_tree->Draw(device_list[idev],bmw_only_cut+burst_cut,"goff");
      if(npt!=0){
	h1dptr= (TH1D*)gDirectory->FindObject("htemp");
	fStat_bmw[idev].mean = h1dptr->GetMean();
	fStat_bmw[idev].error = h1dptr->GetMeanError();
	fStat_bmw[idev].rms = h1dptr->GetRMS();
	fStat_bmw[idev].num_samples = h1dptr->GetEntries();
      }else{
	cout << " ** Empty channel " << device_list[idev] << " with bmw-only cut. "  << endl;
	fStat_bmw[idev] = fStat_zero;
      }
      
      for(int icoil=0;icoil<nCoil;icoil++){
	TString bmwobj_cut = Form("&&yield_bmwobj==%d",icoil+1);
	npt = mul_tree->Draw(device_list[idev],bmw_only_cut+burst_cut+bmwobj_cut,"goff");
	if(npt!=0){
	  h1dptr= (TH1D*)gDirectory->FindObject("htemp");
	  fStat_bmw_by_coil[idev*nCoil+icoil].mean = h1dptr->GetMean();
	  fStat_bmw_by_coil[idev*nCoil+icoil].error = h1dptr->GetMeanError();
	  fStat_bmw_by_coil[idev*nCoil+icoil].rms = h1dptr->GetRMS();
	  fStat_bmw_by_coil[idev*nCoil+icoil].num_samples = h1dptr->GetEntries();
	}else{
	  cout << " ** Empty channel " << device_list[idev] << " with bmw-only cut. "  << endl;
	  fStat_bmw_by_coil[idev*nCoil+icoil] = fStat_zero;
	}
	
      } // end of coil loop

    } // end of idev loop
    mini_tree->Fill();
    mini_plus_tree->Fill();
    mini_bmw_tree->Fill();
    for(int icoil=0;icoil<nCoil;icoil++)
      mini_bmw_tree_by_coil[icoil]->Fill();
    mini_counter++;
  } // end of Burst loop
    
  
  output->cd();
  mini_tree->Write();
  mini_plus_tree->Write();
  mini_bmw_tree->Write();
  for(int icoil=0;icoil<nCoil;icoil++)
    mini_bmw_tree_by_coil[icoil]->Write();

  output->Close();
}

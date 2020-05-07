void CheckBeamPosition(Int_t slug);

void CheckBeamPosition(){
  for(int i=1;i<=94;i++)
    CheckBeamPosition(i);
}

void CheckBeamPosition(Int_t slug){
  TString label = Form("slug%d",slug);
  TString filename  = label+".list";
  TString path = "./prex-runlist/simple_list/";
  FILE *runlist = fopen((path+filename).Data(),"r");
  
  if(runlist==NULL){
    cout << " -- Error: runlist " 
	 << filename
	 << " is not found! " << endl;
    return;
  }

  TFile *output = TFile::Open("./rootfiles/BPMCheck_"+label+".root","RECREATE");
  TString fChannelList[] = { "yield_bpm4aX","yield_bpm4eX","yield_bpm4aY","yield_bpm4eY",
			     "yield_bpm12X","yield_usl","yield_usr","yield_dsl","yield_dsr"};
  Int_t nch = sizeof(fChannelList)/sizeof(*fChannelList);
  typedef struct {Double_t mean,rms,ll,ul,num_samples;} DATA;
  vector<DATA> fData(nch);
  typedef struct {Double_t ppm,ppb,nm,um,mm;} UNIT;
  UNIT aUnit; 
  aUnit.ppm = 1e-6;
  aUnit.ppb = 1e-9;
  aUnit.mm = 1.0;
  aUnit.um = 1e-3;
  aUnit.nm = 1e-6;
  TTree *output_tree = new TTree("T","");
  output_tree->Branch("unit",&aUnit,"ppm/D:ppb:nm:um:mm");
  Int_t fRun;
  output_tree->Branch("run",&fRun,"run/I");
  for(int ich=0;ich<nch;ich++){
    output_tree->Branch(fChannelList[ich],&fData[ich],"mean/D:rms:ll:ul:num_samples");
  }
  // ===== First-Pass
  while(!feof(runlist)){
    Int_t run_number=0;
    Int_t seg_number=0;
    fscanf(runlist,"%d/n",&run_number);
    if(run_number==0)
      continue;
    TString qw_path = TString(gSystem->Getenv("QW_ROOTFILES") );
    TString rootfile_name = Form("prexPrompt_pass1_%d.%03d.root",
				 run_number,seg_number);
    TFile *this_file;
    while(gSystem->AccessPathName(qw_path+rootfile_name)==0){
      this_file = TFile::Open(qw_path+rootfile_name);
      cout << " -- Found ROOT file ";
      cout << this_file->GetName() << endl;
      
      TTree *mul_tree = (TTree*)this_file->Get("mul");
      for(int ich=0;ich<nch;ich++){
	TString device_name = fChannelList[ich];
	int npt=mul_tree->Draw(device_name,"ErrorFlag==0","goff");
	if(npt==0){
	  fData[ich].mean = -9999;
	  fData[ich].rms = -9999;
	  fData[ich].ll = -9999;
	  fData[ich].ul = -9999;
	}else{
	  TH1D *h1d_ptr = (TH1D*)gDirectory->FindObject("htemp");
	  fData[ich].mean = h1d_ptr->GetMean();
	  fData[ich].rms = h1d_ptr->GetRMS();
	  fData[ich].ll = h1d_ptr->GetBinCenter(h1d_ptr->FindFirstBinAbove(0));
	  fData[ich].ul = h1d_ptr->GetBinCenter(h1d_ptr->FindLastBinAbove(0));
	}
	fData[ich].num_samples = npt;
      } // end of channel list loop
      fRun = run_number;
      output_tree->Fill();
      this_file->Close();
      seg_number++;
      rootfile_name = Form("prexPrompt_pass1_%d.%03d.root",
			   run_number,seg_number);
    } // end of file split search
  } // end of runlist loop


  output->cd();
  cout << " -- Writing Tree " << endl;
  output_tree->Write();
  cout << " -- Writing Tree : Done " << endl;
  cout << " -- Closing " << output->GetName() << endl;
  output->Close();
}


void GetAsymmetryByPolarityFromCut(Int_t slug);
void GetAsymmetryByPolarityFromCut();
void GetAsymmetryByPolarityFromCut(){
  for(int i=1;i<=94;i++)
    GetAsymmetryByPolarityFromCut(i);
}
void GetAsymmetryByPolarityFromCut(Int_t slug){
  TString qwrootfile_path = "/media/yetao/prex/PREXII-respin1/";

  TString list_name = Form("./beamoff_cut/slug%d.txt",slug);
  ifstream prex_runlist;
  prex_runlist.open(list_name.Data());
  if(!prex_runlist.is_open())
    return;
  TString output_filename = Form("./rootfiles/slug%d_beamoff.root",slug);
  TFile *output = TFile::Open(output_filename,"RECREATE");
  TTree *mini_pos = new TTree("pos","");
  TTree *mini_neg = new TTree("neg","");
  TTree *mini_null = new TTree("neutral","");
  TTree *mini_norm = new TTree("normal","");
  
  vector<TString> device_list={"yield_usl","yield_usr","yield_dsl","yield_dsr",
			       "yield_atl1","yield_atl2","yield_atr1","yield_atr2",
			       "yield_bpm4aWS","yield_bpm4eWS",
			       "yield_bpm11WS","yield_bpm12WS","yield_bpm1WS",
			       "yield_bcm_an_us","yield_bcm_an_ds",
			       "yield_bcm_dg_us","yield_bcm_dg_ds",
			       "yield_cav4cQ","yield_bcm_an_ds3",
			       "diff_usl","diff_usr","diff_dsl","diff_dsr",
			       "diff_atl1","diff_atl2","diff_atr1","diff_atr2",
			       "diff_bpm4aWS","diff_bpm4eWS",
			       "diff_bpm11WS","diff_bpm12WS","diff_bpm1WS",
			       "diff_bcm_an_us","diff_bcm_an_ds",
			       "diff_bcm_dg_us","diff_bcm_dg_ds",
			       "diff_cav4cQ","diff_bcm_an_ds3"};

  Int_t ndev = device_list.size();
  Int_t run_number = 0;  
  Int_t mini_id = 0;
  TString line_string;
  mini_pos->Branch("run",&run_number);
  mini_neg->Branch("run",&run_number);
  mini_null->Branch("run",&run_number);
  mini_norm->Branch("run",&run_number);
  mini_pos->Branch("mini",&mini_id);
  mini_neg->Branch("mini",&mini_id);
  mini_null->Branch("mini",&mini_id);
  mini_norm->Branch("mini",&mini_id);
  typedef struct {Double_t mean, error,rms;} STAT;
  vector<STAT> fStat_pos(ndev);
  vector<STAT> fStat_neg(ndev);
  vector<STAT> fStat_null(ndev);
  vector<STAT> fStat_norm(ndev);
  STAT fStat_zero;
  fStat_zero.mean=0.0;
  fStat_zero.error=-1;
  fStat_zero.rms=0.0;
  for(int i=0;i<ndev;i++){
    mini_pos->Branch(device_list[i],&fStat_pos[i],"mean/D:err:rms");
    mini_neg->Branch(device_list[i],&fStat_neg[i],"mean/D:err:rms");
    mini_null->Branch(device_list[i],&fStat_null[i],"mean/D:err:rms");
    mini_norm->Branch(device_list[i],&fStat_norm[i],"mean/D:err:rms");
  }
  Int_t counts=0;
  while(line_string.ReadLine(prex_runlist)){
    Ssiz_t first_col = line_string.First(':');
    run_number = ((TString)line_string(0,first_col)).Atoi();
    cout << run_number << endl;
    TString ped_cut = line_string(first_col+1,line_string.Length()-first_col-1);
    ped_cut = "("+ped_cut+")";
    TFile *input = TFile::Open(qwrootfile_path+Form("prexPrompt_pass1_%d.000.root",run_number));
    if(input==NULL)
      continue;
    TTree *mul_tree = (TTree*)input->Get("mul");
    
    for(int idev=0;idev<ndev;idev++){
      if(device_list[idev].Contains("diff_")){
	TString ch_name = device_list[idev];
	ch_name.ReplaceAll("diff_","");
	if(mul_tree->GetBranch("yield_"+ch_name)!=NULL)
	  mul_tree->SetAlias(device_list[idev],
			     "asym_"+ch_name+"*yield_"+ch_name);
      }
    }
    
    TString device_name;
    TH1D *htemp;
    for(int idev=0;idev<ndev;idev++){
      device_name = device_list[idev];
      TString asym_name = device_name;
      asym_name.ReplaceAll("diff_","asym_");
      asym_name.ReplaceAll("yield_","asym_");
      TString not_blinder = Form("&&((%s.Device_Error_Code & 512)!=512)",asym_name.Data());
      if(mul_tree->GetBranch(device_name)==NULL
	 && mul_tree->GetAlias(device_name)==NULL){
	fStat_pos[idev] = fStat_zero;
	fStat_neg[idev] = fStat_zero;
	fStat_null[idev] = fStat_zero;
	fStat_norm[idev] = fStat_zero;
	continue;
      }
      mul_tree->Draw(device_name,"actual_pattern_polarity==1 && "+ped_cut+not_blinder,"goff");
      htemp  = (TH1D*)gDirectory->FindObject("htemp");
      htemp->SetName(Form("htemp%d",counts++));
      fStat_pos[idev].mean = htemp->GetMean();
      fStat_pos[idev].error = htemp->GetMeanError();
      fStat_pos[idev].rms = htemp->GetRMS();
	
      mul_tree->Draw(device_name,"actual_pattern_polarity==0 && "+ped_cut+not_blinder,"goff");
      htemp  = (TH1D*)gDirectory->FindObject("htemp");
      htemp->SetName(Form("htemp%d",counts++));
      fStat_neg[idev].mean = htemp->GetMean();
      fStat_neg[idev].error = htemp->GetMeanError();
      fStat_neg[idev].rms = htemp->GetRMS();

      mul_tree->Draw("2*(actual_pattern_polarity-0.5)*"+device_name,ped_cut+not_blinder,"goff");
      htemp  = (TH1D*)gDirectory->FindObject("htemp");
      htemp->SetName(Form("htemp%d",counts++));
      fStat_null[idev].mean = htemp->GetMean();
      fStat_null[idev].error = htemp->GetMeanError();
      fStat_null[idev].rms = htemp->GetRMS();
      
      mul_tree->Draw(device_name,ped_cut+not_blinder,"goff");
      htemp  = (TH1D*)gDirectory->FindObject("htemp");
      htemp->SetName(Form("htemp%d",counts++));
      fStat_norm[idev].mean = htemp->GetMean();
      fStat_norm[idev].error = htemp->GetMeanError();
      fStat_norm[idev].rms = htemp->GetRMS();
    }
    mini_norm->Fill();
    mini_null->Fill();
    mini_pos->Fill();
    mini_neg->Fill();

    input->Close();
  }
  
  prex_runlist.close();
  output->cd();
  mini_neg->Write();
  mini_pos->Write();
  mini_null->Write();
  mini_norm->Write();
  output->Close();
}

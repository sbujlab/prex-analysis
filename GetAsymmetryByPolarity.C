void GetAsymmetryByPolarity(Int_t slug);
void GetAsymmetryByPolarity();
void GetAsymmetryByPolarity(){
  for(int i=1;i<=94;i++)
    GetAsymmetryByPolarity(i);
}
void GetAsymmetryByPolarity(Int_t slug){
  TString qwrootfile_path = "/media/yetao/prex/PREXII-respin1/";
  TString postpan_path = "/media/yetao/prex/postpan_respin/";

  TString list_name = Form("/home/yetao/workarea/prex-runlist/good_list/slug%d.list",slug);
  FILE* prex_runlist = fopen(list_name.Data(),"r");
  if(prex_runlist==NULL)
    return;
  TString output_filename = Form("./rootfiles/slug%d.root",slug);
  TFile *output = TFile::Open(output_filename,"RECREATE");
  TTree *mini_pos = new TTree("pos","");
  TTree *mini_neg = new TTree("neg","");
  TTree *mini_null = new TTree("neutral","");

  
  vector<TString> device_list={"reg_asym_us_avg","reg_asym_us_dd","reg_asym_usl","reg_asym_usr",
  			       "asym_us_avg","asym_us_dd","asym_usl","asym_usr",
			       "reg_asym_ds_avg","reg_asym_ds_dd","reg_asym_dsl","reg_asym_dsr",
  			       "asym_ds_avg","asym_ds_dd","asym_dsl","asym_dsr",
			       "diff_bpm4aX","diff_bpm4eX","diff_bpm4aY","diff_bpm4eY",
			       "diff_bpm11X","diff_bpm12X","diff_bpmE",
  			       "asym_bcm_an_us","asym_bcm_an_ds","asym_bcm_dg_us","asym_bcm_dg_ds",
  			       "asym_battery1l","asym_battery2l","asym_battery1r","asym_battery2r",
  			       "asym_ch_battery_1","asym_ch_battery_2",
			       "diff_battery1l","diff_battery2l","diff_battery1r","diff_battery2r",
  			       "diff_ch_battery_1","diff_ch_battery_2" };

  // vector<TString> device_list={"asym_battery1l","asym_battery2l","asym_battery1r","asym_battery2r",
  // 			       "asym_ch_battery_1","asym_ch_battery_2"};

  Int_t ndev = device_list.size();
  Int_t run_number = 0;  
  Int_t mini_id = 0;

  mini_pos->Branch("run",&run_number);
  mini_neg->Branch("run",&run_number);
  mini_null->Branch("run",&run_number);
  mini_pos->Branch("mini",&mini_id);
  mini_neg->Branch("mini",&mini_id);
  mini_null->Branch("mini",&mini_id);
  typedef struct {Double_t mean, error,rms;} STAT;
  vector<STAT> fStat_pos(ndev);
  vector<STAT> fStat_neg(ndev);
  vector<STAT> fStat_null(ndev);
  STAT fStat_zero;
  fStat_zero.mean=0.0;
  fStat_zero.error=-1;
  fStat_zero.rms=0.0;
  for(int i=0;i<ndev;i++){
    mini_pos->Branch(device_list[i],&fStat_pos[i],"mean/D:err:rms");
    mini_neg->Branch(device_list[i],&fStat_neg[i],"mean/D:err:rms");
    mini_null->Branch(device_list[i],&fStat_null[i],"mean/D:err:rms");
  }
  Int_t counts=0;
  while(!feof(prex_runlist)){
    run_number =0;
    fscanf(prex_runlist,"%d\n",&run_number);
    cout << run_number << endl;
    TFile *input = TFile::Open(qwrootfile_path+Form("prexPrompt_pass1_%d.000.root",run_number));
    TFile *redfile;
    if(slug<=3)
      redfile = TFile::Open(postpan_path+Form("prexPrompt_%d_000_regress_postpan.root",run_number));
    else
      redfile = TFile::Open(postpan_path+Form("prexPrompt_%d_000_regress_comboBPM.root",run_number));
    
    if(input==NULL || redfile==NULL)
      continue;
    TTree *mul_tree = (TTree*)input->Get("mul");
    TTree *reg_tree = (TTree*)redfile->Get("reg");
    mul_tree->AddFriend(reg_tree);
    
    vector<TString> user_define={"battery1l","battery2l","battery1r","battery2r",
				 "ch_battery_1","ch_battery_2"};
    Int_t nDefine = user_define.size();
    for(int i=0;i<nDefine;i++){
      mul_tree->SetAlias("diff_"+user_define[i],
			 "asym_"+user_define[i]+"*yield_"+user_define[i]);
    }
    
    if(mul_tree->GetBranch("diff_bpm11X")!=NULL && slug>3)
      mul_tree->SetAlias("diff_bpmE","diff_bpm11X+0.4*diff_bpm12X");
    else
      mul_tree->SetAlias("diff_bpmE","diff_bpm12X");
    
    Int_t nMini  = ((TTree*)redfile->Get("mini"))->GetEntries();

    TString device_name;
    TH1D *htemp;
    TString mini_cut;
    for(int imini=0;imini<nMini;imini++){
      mini_cut = Form("&& minirun==%d",imini);
      for(int idev=0;idev<ndev;idev++){
	device_name = device_list[idev];
	
	if(mul_tree->GetBranch(device_name)==NULL
	  && mul_tree->GetAlias(device_name)==NULL){
	  fStat_pos[idev] = fStat_zero;
	  fStat_neg[idev] = fStat_zero;
	  fStat_null[idev] = fStat_zero;
	  continue;
	}
	mul_tree->Draw(device_name,"ok_cut && actual_pattern_polarity==1"+mini_cut,"goff");
	htemp  = (TH1D*)gDirectory->FindObject("htemp");
	htemp->SetName(Form("htemp%d",counts++));
	fStat_pos[idev].mean = htemp->GetMean();
	fStat_pos[idev].error = htemp->GetMeanError();
	fStat_pos[idev].rms = htemp->GetRMS();
	
	mul_tree->Draw(device_name,"ok_cut && actual_pattern_polarity==0"+mini_cut,"goff");
	htemp  = (TH1D*)gDirectory->FindObject("htemp");
	htemp->SetName(Form("htemp%d",counts++));
	fStat_neg[idev].mean = htemp->GetMean();
	fStat_neg[idev].error = htemp->GetMeanError();
	fStat_neg[idev].rms = htemp->GetRMS();

	mul_tree->Draw("2*(actual_pattern_polarity-0.5)*"+device_name,"ok_cut"+mini_cut,"goff");
	htemp  = (TH1D*)gDirectory->FindObject("htemp");
	htemp->SetName(Form("htemp%d",counts++));
	fStat_null[idev].mean = htemp->GetMean();
	fStat_null[idev].error = htemp->GetMeanError();
	fStat_null[idev].rms = htemp->GetRMS();

      }
      
      mini_id = imini;
      mini_null->Fill();
      mini_pos->Fill();
      mini_neg->Fill();
    }

    input->Close();
    redfile->Close();
  }
  
  fclose(prex_runlist);
  output->cd();
  mini_neg->Write();
  mini_pos->Write();
  mini_null->Write();
  output->Close();
}

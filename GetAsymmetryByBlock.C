#include "device_list.hh"
void GetAsymmetryByBlock(Int_t slug);
void GetAsymmetryByBlock();
void GetAsymmetryByBlock(){
  for(int i=1;i<=94;i++)
    GetAsymmetryByBlock(i);
}
void GetAsymmetryByBlock(Int_t slug){
  TString qwrootfile_path = "$QW_ROOTFILES/";
  TString postpan_path = " /lustre/expphy/volatile/halla/parity/postpan_respin/";

  TString list_name = Form("./prex-runlist/simple_list/slug%d.list",slug);
  FILE* prex_runlist = fopen(list_name.Data(),"r");
  if(prex_runlist==NULL)
    return;
  TString output_filename = Form("./rootfiles/slug%d_by_block.root",slug);
  TFile *output = TFile::Open(output_filename,"RECREATE");
  TTree *mini_tree = new TTree("norm","");
  
  Int_t ndev = device_list.size(); // >> "device_list.hh"
  Int_t run_number = 0;  
  Int_t mini_id = 0;

  mini_tree->Branch("run",&run_number);
  mini_tree->Branch("mini",&mini_id);
  typedef struct {Double_t mean, error,rms;} STAT;
  TString block_format[4]={".block0",".block1",".block2",".block3"};
  vector<STAT> fStat_mini(ndev*4);
  STAT fStat_zero;
  fStat_zero.mean=0.0;
  fStat_zero.error=-1;
  fStat_zero.rms=0.0;
  for(int i=0;i<ndev;i++)
    for(int iblock=0;iblock<4;iblock++)
      mini_tree->Branch(device_list[i]+block_format[iblock],
			&fStat_mini[4*i+iblock],"mean/D:err:rms");

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
    
    if(input==NULL)
      continue;
    TTree *mul_tree = (TTree*)input->Get("mul");
    TTree *mulc_tree = (TTree*)input->Get("mulc");
    TTree *mulc_dit_tree = (TTree*)input->Get("mulc_dit");
    TTree *mulc_dit_combo_tree = (TTree*)input->Get("mulc_dit_combo");
    TTree *reg_tree = (TTree*)redfile->Get("reg");
    mul_tree->AddFriend(reg_tree);
    mul_tree->AddFriend(mulc_tree);
    mul_tree->AddFriend(mulc_dit_tree);
    mul_tree->AddFriend(mulc_dit_combo_tree);
    // ==== Set Up Alias for new combination 
    vector<TString> user_define={"battery1l","battery2l","battery1r","battery2r",
				 "ch_battery_1","ch_battery_2"};
    Int_t nDefine = user_define.size();
    for(int i=0;i<nDefine;i++){
      for(int iblock=0;iblock<4;iblock++){
	mul_tree->SetAlias("diff_"+user_define[i]+block_format[iblock],
			   "asym_"+user_define[i]+block_format[iblock]+"*yield_"+user_define[i]+block_format[iblock]);
      }
    }
    // FIXME: remove this part after respin2
    // =======>>
    for(int iblock=0;iblock<4;iblock++){
      mul_tree->SetAlias("dit_asym_us_dd_div2"+block_format[iblock],"0.5*dit_asym_us_dd"+block_format[iblock]);
      mul_tree->SetAlias("dit_asym_ds_dd_div2"+block_format[iblock],"0.5*dit_asym_ds_dd"+block_format[iblock]);
    }
    // =======<<
    vector<TString> det_list={"asym_usl","asym_usr","asym_us_avg","asym_us_dd",
			      "asym_dsl","asym_dsr","asym_ds_avg","asym_ds_dd",
			      "asym_atl1","asym_atr1","asym_atl2","asym_atr2",
			      "asym_atl_avg","asym_atl_dd","asym_atr_avg","asym_atr_dd",
			      "asym_at1_avg","asym_at1_dd","asym_at2_avg","asym_at2_dd",
			      "asym_atl1r2_avg","asym_atl1r2_dd",
			      "asym_atr1l2_avg","asym_atr1l2_dd",
			      "asym_atl_dd_atr_dd_avg","asym_atl_dd_atr_dd_dd",
			      "asym_atl_avg_atr_avg_avg","asym_atl_avg_atr_avg_dd"};
    Int_t nDet = det_list.size();
    for(int idet=0;idet<nDet;idet++){
      if(mul_tree->GetBranch(det_list[idet])==NULL){
	cout << "SetAlias(): " <<det_list[idet] << " is not found in mul tree " << endl;
	continue;
      }
      // FIXME : uncomment later
      // if(mul_tree->GetBranch("reg_"+det_list[idet])!=NULL){
      // 	for(int iblk=0;iblk<4;iblk++)
      // 	  mul_tree->SetAlias("reg_corr_"+det_list[idet]+block_format[iblk],
      // 			     Form("%s-reg_%s",(det_list[idet]+block_format[iblk]).Data(),
      // 				  (det_list[idet]+block_format[iblk]).Data()));
      // }

      if(mul_tree->GetBranch("dit_"+det_list[idet])!=NULL){
      // FIXME: remove this part after respin2
      // =======>>
	for(int iblk=0;iblk<4;iblk++){
	  if(det_list[idet]=="asym_us_dd"
	     || det_list[idet]=="asym_ds_dd")
	    mul_tree->SetAlias("dit_corr_"+det_list[idet]+block_format[iblk],
			       Form("%s-0.5*dit_%s",
				    (det_list[idet]+block_format[iblk]).Data(),
				    (det_list[idet]+block_format[iblk]).Data()));
	  else      // =======<<
	    mul_tree->SetAlias("dit_corr_"+det_list[idet]+block_format[iblk],
			       Form("%s-dit_%s",
				    (det_list[idet]+block_format[iblk]).Data(),
				    (det_list[idet]+block_format[iblk]).Data()));
	}
      }
    }
    //---- BPM
    vector<TString> bpm_list={"diff_bpm4a","diff_bpm4e","diff_bpm1","diff_bpm4ac","diff_bpm4ec",
			      "diff_bpm11","diff_bpm12"};
    Int_t nbpm = bpm_list.size();
    for(int i=0;i<nbpm;i++){
      if(mul_tree->GetBranch(bpm_list[i]+"X")==NULL)
	continue;
      for(int iblk=0;iblk<4;iblk++){
	mul_tree->SetAlias(bpm_list[i]+"X_unrotated"+block_format[iblk],
			   Form("sqrt(2)/2*(%sX%s +%sY%s)",
				bpm_list[i].Data(),block_format[iblk].Data(),
				bpm_list[i].Data(),block_format[iblk].Data()));
	mul_tree->SetAlias(bpm_list[i]+"Y_unrotated"+block_format[iblk],
			   Form("sqrt(2)/2*(%sY%s -%sX%s)",
				bpm_list[i].Data(),block_format[iblk].Data(),
				bpm_list[i].Data(),block_format[iblk].Data()));
      }
    }
    

    if(mul_tree->GetBranch("diff_bpm11X")!=NULL && slug>3){
      for(int iblk=0;iblk<4;iblk++){
	mul_tree->SetAlias("diff_bpmE"+block_format[iblk],
			   "diff_bpm11X"+block_format[iblk]+"+0.4*diff_bpm12X"+block_format[iblk]);
      }
    }else{
      for(int iblk=0;iblk<4;iblk++){
	mul_tree->SetAlias("diff_bpmE"+block_format[iblk],
			   "diff_bpm12X"+block_format[iblk]);
      }
    }
    // --- Position DD and AVG
    vector<TString> bpmAE_list={"4a","4e","4ac","4ec"};
    Int_t nbpmAE = bpmAE_list.size();
    for(int i=0;i<nbpmAE;i++){
      for(int j=i+1;j<nbpmAE;j++){
	for(int iblk=0;iblk<4;iblk++){
	  mul_tree->SetAlias(Form("dd_bpm%sX_bpm%sX",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sX%s-diff_bpm%sX%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));
			   
	  mul_tree->SetAlias(Form("dd_bpm%sY_bpm%sY",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sY%s-diff_bpm%sY%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));

	  mul_tree->SetAlias(Form("avg_bpm%sX_bpm%sX",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sX%s+diff_bpm%sX%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));
			   
	  mul_tree->SetAlias(Form("avg_bpm%sY_bpm%sY",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sY%s+diff_bpm%sY%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));

	  mul_tree->SetAlias(Form("dd_bpm%sX_bpm%sX_unrotated",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sX_unrotated%s-diff_bpm%sX_unrotated%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));
			   
	  mul_tree->SetAlias(Form("dd_bpm%sY_bpm%sY_unrotated",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sY_unrotated%s-diff_bpm%sY_unrotated%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));

	  mul_tree->SetAlias(Form("avg_bpm%sX_bpm%sX_unrotated",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sX_unrotated%s+diff_bpm%sX_unrotated%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));
			   
	  mul_tree->SetAlias(Form("avg_bpm%sY_bpm%sY_unrotated",
				  bpmAE_list[i].Data(),bpmAE_list[j].Data())+block_format[iblk],
			     Form("0.5*(diff_bpm%sY_unrotated%s+diff_bpm%sY_unrotated%s)",
				  bpmAE_list[i].Data(),block_format[iblk].Data(),
				  bpmAE_list[j].Data(),block_format[iblk].Data()));

	}
      }
    }
    // --- WireSum DD  and AVG
    vector<TString> bpmws_list={"4a","4e","4ac","4ec","1","11","12"};
    Int_t nWS = bpmws_list.size();
    for(int i=0;i<nWS;i++){
      for(int j=i+1;j<nWS;j++){
	if( mul_tree->GetBranch("asym_bpm"+bpmws_list[i]+"WS")==NULL||
	    mul_tree->GetBranch("asym_bpm"+bpmws_list[j]+"WS")==NULL)
	  continue;
	for(int iblk=0;iblk<4;iblk++){
	  mul_tree->SetAlias(Form("dd_bpm%sWS_bpm%s_WS",
				  bpmws_list[i].Data(),bpmws_list[j].Data())+block_format[iblk],
			     Form("0.5*(asym_bpm%sWS%s-asym_bpm%sWS%s)",
				  bpmws_list[i].Data(),block_format[iblk].Data(),
				  bpmws_list[j].Data(),block_format[iblk].Data()));
	
	  mul_tree->SetAlias(Form("avg_bpm%sWS_bpm%s_WS",
				  bpmws_list[i].Data(),bpmws_list[j].Data())+block_format[iblk],
			     Form("0.5*(asym_bpm%sWS%s+asym_bpm%sWS%s)",
				  bpmws_list[i].Data(),block_format[iblk].Data(),
				  bpmws_list[j].Data(),block_format[iblk].Data()));

	}
      }
    }
    //---- BCM
    for(int iblk=0;iblk<4;iblk++)
      mul_tree->SetAlias("asym_bcm_an_new"+block_format[iblk],
			 "2.0/3.0*asym_bcm_an_us"+block_format[iblk]+"+1.0/3.0*asym_bcm_an_ds"+block_format[iblk]);

    vector<TString> bcm_list={"asym_bcm_an_us","asym_bcm_an_ds","asym_bcm_dg_us","asym_bcm_dg_ds",
			      "asym_bcm_an_new","asym_cav4cQ"};
    Int_t nbcm = bcm_list.size();
    for(int i=0;i<nbcm;i++){
      TString chA_name = bcm_list[i];
      chA_name.ReplaceAll("asym_","");
      for(int j=i+1;j<nbcm;j++){
	TString chB_name = bcm_list[j];
	chB_name.ReplaceAll("asym_","");
	for(int iblk=0;iblk<4;iblk++){
	  mul_tree->SetAlias(Form("dd_%s_%s",chA_name.Data(),chB_name.Data())+block_format[iblk],
			     Form("0.5*(%s-%s)",
				  (bcm_list[i]+block_format[iblk]).Data(),
				  (bcm_list[j]+block_format[iblk]).Data()));
	  mul_tree->SetAlias(Form("avg_%s_%s",chA_name.Data(),chB_name.Data())+block_format[iblk],
			     Form("0.5*(%s+%s)",
				  (bcm_list[i]+block_format[iblk]).Data(),
				  (bcm_list[j]+block_format[iblk]).Data()));
	}
      }
    }
    // ==== Done with alias for new combination
    
    Int_t nMini  = ((TTree*)redfile->Get("mini"))->GetEntries();

    TString device_name;
    TH1D *htemp;
    TString mini_cut;
    for(int imini=0;imini<nMini;imini++){
      mini_cut = Form("&& minirun==%d",imini);
      for(int idev=0;idev<ndev;idev++){
	device_name = device_list[idev];

	// FIXME: Skip postpan block output at this point 
	if(device_name.Contains("reg")){
	  for(int iblk=0;iblk<4;iblk++)
	    fStat_mini[4*idev+iblk] = fStat_zero;
	  // cout << " -- Warning: " << device_name << " is not found " <<endl;
	  continue;
	}

	if(mul_tree->GetBranch(device_name)==NULL
	  && mul_tree->GetAlias(device_name)==NULL
	   && mul_tree->GetAlias(device_name+".block0")==NULL){
	  for(int iblk=0;iblk<4;iblk++)
	    fStat_mini[4*idev+iblk] = fStat_zero;
	  cout << " -- Warning: " << device_name << " is not found " <<endl;
	  continue;
	}
	
	// mul_tree->Draw(device_name,"ok_cut && actual_pattern_polarity==1"+mini_cut,"goff");
	// htemp  = (TH1D*)gDirectory->FindObject("htemp");
	// htemp->SetName(Form("htemp%d",counts++));
	// fStat_pos[idev].mean = htemp->GetMean();
	// fStat_pos[idev].error = htemp->GetMeanError();
	// fStat_pos[idev].rms = htemp->GetRMS();
	
	// mul_tree->Draw(device_name,"ok_cut && actual_pattern_polarity==0"+mini_cut,"goff");
	// htemp  = (TH1D*)gDirectory->FindObject("htemp");
	// htemp->SetName(Form("htemp%d",counts++));
	// fStat_neg[idev].mean = htemp->GetMean();
	// fStat_neg[idev].error = htemp->GetMeanError();
	// fStat_neg[idev].rms = htemp->GetRMS();

	// mul_tree->Draw("2*(actual_pattern_polarity-0.5)*"+device_name,"ok_cut"+mini_cut,"goff");
	// htemp  = (TH1D*)gDirectory->FindObject("htemp");
	// htemp->SetName(Form("htemp%d",counts++));
	// fStat_null[idev].mean = htemp->GetMean();
	// fStat_null[idev].error = htemp->GetMeanError();
	// fStat_null[idev].rms = htemp->GetRMS();
	
	for(int iblk=0;iblk<4;iblk++){
	  mul_tree->Draw(device_name+block_format[iblk],"ok_cut"+mini_cut,"goff");
	  htemp  = (TH1D*)gDirectory->FindObject("htemp");
	  htemp->SetName(Form("htemp%d",counts++));
	  fStat_mini[4*idev+iblk].mean = htemp->GetMean();
	  fStat_mini[4*idev+iblk].error = htemp->GetMeanError();
	  fStat_mini[4*idev+iblk].rms = htemp->GetRMS();
	}

      }
      
      mini_id = imini;
      mini_tree->Fill();
    }

    input->Close();
    redfile->Close();
  }
  
  fclose(prex_runlist);
  output->cd();
  mini_tree->Write();
  output->Close();
}

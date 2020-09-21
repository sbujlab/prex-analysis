#include "LoadRunInfoMap.cc"
#include "TaRunInfo.cc"

map<TString,TString> fLinks={ {"mini","PostpanMerged"},
			      {"burst_mulc_dit","MergedSum"},
			      {"mini_lagr6bpm","MergedLagrange"},
			      {"mini_reg6bpm","MergedLagrange"},
			      {"mini_lagrall","MergedLagrange"},
			      {"mini_regall","MergedLagrange"},
			      {"mini_eiglagrall_tr","MergedLagrange_trunc"}};

map<TString,TString > fMaps ={ {"dit","burst_mulc_dit"},
			       {"reg_5bpm","mini"},
			       {"reg_6bpm","mini_reg6bpm"},
			       {"lagr_6bpm","mini_lagr6bpm"},
			       {"reg_all","mini_regall"},
			       {"lagr_all","mini_lagrall"},
			       {"lagr_all_trunc","mini_eiglagrall_tr"}};
  
vector< pair<TString,TString> > fSets ={{"lagr_all_trunc","lagr_all"},
					{"dit","reg_5bpm"},
					{"dit","lagr_6bpm"},
					{"dit","lagr_all"},
					{"dit","reg_all"},
					{"reg_5bpm","reg_6bpm"},
					{"reg_5bpm","reg_all"},
					{"lagr_6bpm","reg_6bpm"},
					{"lagr_all","reg_all"},
					{"lagr_6bpm","lagr_all"}};

void MinirunDiffBySlug(){
  auto iter_set = fSets.begin();
  while(iter_set!=fSets.end()){
    TString set_name1 = (*iter_set).first;
    TString set_name2 = (*iter_set).second;
    TString tree_name1 = fMaps[set_name1];
    TString tree_name2 = fMaps[set_name2];
    TString file_stem1 = fLinks[tree_name1];
    TString file_stem2 = fLinks[tree_name2];
    TString output_filename = Form("prex_minidiff_%s_%s.root",
				   set_name1.Data(),set_name2.Data());
    TFile *output_file = TFile::Open(output_filename,"RECREATE");
    TTree *output_tree = new TTree("slug","");
    TString prefix1,prefix2;
    if(set_name1.Contains("reg")){
      prefix1 = "reg";
    }else if(set_name1.Contains("lagr")){
      prefix1 = "lagr";
    }else if(set_name1.Contains("dit")){
      prefix1 = "dit";
    }
    if(set_name2.Contains("reg")){
      prefix2 = "reg";
    }else if(set_name2.Contains("lagr")){
      prefix2 = "lagr";
    }else if(set_name2.Contains("dit")){
      prefix2 = "dit";
    }

    Int_t fSlug;
    Int_t fSign,fArm_flag;
    output_tree->Branch("slug",&fSlug,"slug/I");
    output_tree->Branch("sign",&fSign,"sign/I");
    output_tree->Branch("arm_flag",&fArm_flag);
    vector<Double_t> delta(4,0);
    vector<Double_t> distance(4,0);
    vector<Double_t> delta_error(4,0);
    TString det_array[]={"us_avg","usr","usl","us_dd"};
    Int_t nDet=4;
    for(int idet=0;idet<nDet;idet++){
      output_tree->Branch("delta_"+det_array[idet],&delta[idet]);
      output_tree->Branch("dist_"+det_array[idet],&distance[idet]);
      output_tree->Branch("delta_"+det_array[idet]+"_err",&delta_error[idet]);
    }
    Double_t delta_md, distance_md,delta_md_err;
    output_tree->Branch("delta_Adet",&delta_md);
    output_tree->Branch("dist_Adet",&distance_md);
    output_tree->Branch("delta_Adet_err",&delta_md_err);
    TF1 *fgaus = new TF1("fgaus","gaus(0)",-1e6,1e6);
    TCanvas *chist= new TCanvas ("chist","chist",1200,600);
    chist->Divide(4,2);
    chist->Print(Form("prex_minidiff_hist_%s_%s.pdf[",set_name1.Data(),set_name2.Data()));
    for(int islug=1;islug<=94;islug++){

      // Load Input Trees
      TFile *burst_file = TFile::Open(Form("./treeMergeOutput/MergedSum_slug%d.root",islug));
      TTree *burst_tree  = (TTree*)burst_file->Get("burst");
      TFile *input_file1 = TFile::Open(Form("./treeMergeOutput/%s_slug%d.root",
					    file_stem1.Data(),islug));
      TFile *input_file2 = TFile::Open(Form("./treeMergeOutput/%s_slug%d.root",
					    file_stem2.Data(),islug));
      TTree *input_tree1 = (TTree*)input_file1->Get(tree_name1);
      TTree *input_tree2 = (TTree*)input_file2->Get(tree_name2);
      burst_tree->AddFriend(input_tree1,set_name1);
      burst_tree->AddFriend(input_tree2,set_name2);
      if(tree_name1=="burst_mulc_dit" || tree_name2=="burst_mulc_dit"){
	TTree * dit_combo_tree = (TTree*)burst_file->Get("burst_mulc_dit_combo");
	burst_tree->AddFriend(dit_combo_tree,"dit");
      }
      for(int arm_switch=0;arm_switch<3;arm_switch++){
	map< Int_t, TaRunInfo > fRunInfoMap = LoadRunInfoMapBySlug(islug,arm_switch);
	if(fRunInfoMap.size()==0)
	  continue;
	for(int idet=0;idet<nDet;idet++){
	  // Load Run cuts
	  TString myCut =  LoadArmFlagCuts(fRunInfoMap,det_array[idet]);
	  TH1D *hbuff;
	  Int_t npt = 0;
	  npt = burst_tree->Draw(Form("(%s.%s_asym_%s  - %s.%s_asym_%s)/ppb",
				      set_name1.Data(),prefix1.Data(),det_array[idet].Data(),
				      set_name2.Data(),prefix2.Data(),det_array[idet].Data()),
				 myCut, "goff");
	  
	  chist->cd(idet+1+4);
	  gPad->Clear();
	  chist->cd(idet+1);
	  gPad->Clear();
	  gStyle->SetStatH(0.3);
	  gStyle->SetStatW(0.3);
	  gStyle->SetOptFit(1);
	  hbuff = (TH1D*)gDirectory->FindObject("htemp");
	  if(hbuff!=NULL && npt!=0){
	    double mean = hbuff->GetMean();
	    double rms = hbuff->GetRMS();
	    delta[idet] = mean;
	    double min_bin = hbuff->GetXaxis()->GetXmin();
	    double max_bin = 5*(mean-min_bin);
	    double bin_width = 2*hbuff->GetXaxis()->GetBinWidth(1);
	    int nbins = (max_bin-min_bin)/bin_width;

	    TH1D *hnew = new TH1D("hnew","",nbins,min_bin,max_bin);
	    burst_tree->Draw(Form("(%s.%s_asym_%s  - %s.%s_asym_%s)/ppb>>hnew",
				  set_name1.Data(),prefix1.Data(),det_array[idet].Data(),
				  set_name2.Data(),prefix2.Data(),det_array[idet].Data()),
			     myCut, "goff");
	    hnew->SetTitle(Form("slug%d:  ",islug)+det_array[idet]+";ppb;");
	    hnew->Draw();
	    fgaus->SetParameter(1,mean);
	    fgaus->SetParameter(2,rms);
	    hnew->Fit("fgaus","QIR","",mean-3*rms,mean+3*rms);
	    hnew->SetName("delta_A");
	    
	    // Pull plots
	    chist->cd(idet+1+4);
	    burst_tree->Draw(Form("(%s.%s_asym_%s  - %s.%s_asym_%s)/ppb:sqrt(pow(%s.%s_asym_%s.err,2)- pow(%s.%s_asym_%s.err,2))/ppb",
				  set_name1.Data(),prefix1.Data(),det_array[idet].Data(),
				  set_name2.Data(),prefix2.Data(),det_array[idet].Data(),
				  set_name1.Data(),prefix1.Data(),det_array[idet].Data(),
				  set_name2.Data(),prefix2.Data(),det_array[idet].Data()),
			     myCut, "goff");

	    TH1D *hpull =  new TH1D("hpull","",34,-5,12);
	    double *yval = burst_tree->GetV1();
	    double *errval = burst_tree->GetV2();
	    for(int i=0;i<npt;i++){
	      hpull->Fill( yval[i]/errval[i]) ;
	    }
	    hpull->Draw();
	    hpull->SetName("delta_A:PullFit");
	    // hpull->Fit("gaus","QR","",-3,3);
	  } else
	    delta[idet] = -1;
	  
	  npt = burst_tree->Draw(Form("pow(%s.%s_asym_%s-%s.%s_asym_%s,2)",
				      set_name1.Data(),prefix1.Data(),det_array[idet].Data(),
				      set_name2.Data(),prefix2.Data(),det_array[idet].Data()),
				 myCut, "goff");
	  hbuff = (TH1D*)gDirectory->FindObject("htemp");
	  if(hbuff!=NULL && npt!=0){
	    distance[idet] = TMath::Sqrt(hbuff->GetMean())*1e9;
	    delta_error[idet]= distance[idet]/sqrt(npt);
	    hbuff->Draw();
	    hbuff->SetName("delta_A_sq");
	  }else{
	    distance[idet] = -1;
	    delta_error[idet]=-1;
	  }

	} // loop over detector channels
	// Fill output tree with slug info
	fSlug = islug;
	fArm_flag = arm_switch;
	fSign = (*(fRunInfoMap.begin())).second.GetSign();
	delta_md = delta[arm_switch];
	distance_md = distance[arm_switch];
	delta_md_err = delta_error[arm_switch];
	chist->Print(Form("prex_minidiff_hist_%s_%s.pdf",set_name1.Data(),set_name2.Data()));
	output_tree->Fill();
      } // end of arm switch loop;

      
      input_file1->Close();
      input_file2->Close();
      burst_file->Close();
    } // end of loop over slugs
    chist->Print(Form("prex_minidiff_hist_%s_%s.pdf]",set_name1.Data(),set_name2.Data()));
    iter_set++;

    output_file->cd();
    output_tree->Write();
    output_file->Close();
  } // iteration of comparison sets

}

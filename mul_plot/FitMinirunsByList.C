vector<Int_t> LoadRunList(TString key);
void FitMinirunsByList(){
  TString list_array[]={"50uA","50uA_1arm",
			"70uA_120Hz","70uA_120Hz_1arm",
			"70uA_240Hz","70uA_240Hz_1arm",
			"85uA"};
  Int_t npt=  sizeof(list_array)/sizeof(*list_array);
  struct STAT{ double mean, error, rms,num_samples;};
  double* fGrandMean_reg =  new double[npt];
  double* fGrandError_reg =  new double[npt];
  double* fGrandChisq_reg =  new double[npt];
  int* fGrandNDF_reg =  new int[npt];

  double* fGrandMean_lagr =  new double[npt];
  double* fGrandError_lagr =  new double[npt];
  double* fGrandChisq_lagr =  new double[npt];
  int* fGrandNDF_lagr =  new int[npt];

  double* fnpt = new double[npt];

  TFile *merged_file = TFile::Open("../tree_merge/rootfiles/MergedLagrange_all_slugs.root");
  TTree *lagr_tree = (TTree*)merged_file->Get("mini_lagrall");
  TTree *reg_tree = (TTree*)merged_file->Get("mini_regall");
  lagr_tree->AddFriend(reg_tree);
  lagr_tree->AddFriend("mini_info","../tree_merge/rootfiles/runinfo_all_slugs.root");
    
  Int_t run_number;
  Int_t arm_flag;
  Int_t run_flag;
  Int_t kSign;
  lagr_tree->SetBranchAddress("kGood",&run_flag);
  lagr_tree->SetBranchAddress("run",&run_number);
  lagr_tree->SetBranchAddress("arm_flag",&arm_flag);
  lagr_tree->SetBranchAddress("spin",&kSign);
  vector<STAT> fRegStatArray(3);
  vector<STAT> fLagrStatArray(3);
  vector<TString> det_name ={"us_avg","usr","usl"};
  for(int idet=0;idet<3;idet++){
    lagr_tree->SetBranchAddress("reg_asym_"+det_name[idet],&fRegStatArray[idet]);
    lagr_tree->SetBranchAddress("lagr_asym_"+det_name[idet],&fLagrStatArray[idet]);
  }
  gStyle->SetOptFit(1111);
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->SaveAs("mini_fit_sorted.pdf[");

  for(int ipt=0;ipt<npt;ipt++){
    vector<Int_t> fRunlist = LoadRunList(list_array[ipt]);
    vector<Double_t> fMean_reg;
    vector<Double_t> fError_reg;
    vector<Double_t> fMean_lagr;
    vector<Double_t> fError_lagr;
    vector<Double_t> fIndex;
    double fCounter =0;
    Int_t kEntries = lagr_tree->GetEntries();
    for(int ievt=0;ievt<kEntries;ievt++){
      lagr_tree->GetEntry(ievt);
      if(run_flag!=1)
	continue;
      if(find(fRunlist.begin(),fRunlist.end(),run_number)==fRunlist.end())
	continue;
      fMean_reg.push_back( kSign*fRegStatArray[arm_flag].mean * 1e9);
      fError_reg.push_back( fRegStatArray[arm_flag].error * 1e9);
      fMean_lagr.push_back( kSign*fLagrStatArray[arm_flag].mean * 1e9);
      fError_lagr.push_back( fLagrStatArray[arm_flag].error * 1e9);
      fIndex.push_back(fCounter++);
    }
    Int_t nmini = fIndex.size();
    double* fMean_reg_ptr = new double[nmini];
    double* fError_reg_ptr = new double[nmini];
    double* fMean_lagr_ptr = new double[nmini];
    double* fError_lagr_ptr = new double[nmini];
    double* fIndex_ptr = new double[nmini];
    for(int imini=0;imini<nmini;imini++){
      fMean_reg_ptr[imini] = fMean_reg[imini];
      fError_reg_ptr[imini] = fError_reg[imini];
      fMean_lagr_ptr[imini] = fMean_lagr[imini];
      fError_lagr_ptr[imini] = fError_lagr[imini];
      fIndex_ptr[imini] = fIndex[imini];
    }
    TGraphErrors *ger_reg = new TGraphErrors(nmini,fIndex_ptr,fMean_reg_ptr,0,fError_reg_ptr);
    TGraphErrors *ger_lagr = new TGraphErrors(nmini,fIndex_ptr,fMean_lagr_ptr,0,fError_lagr_ptr);
    
    ger_reg->SetMarkerStyle(20);
    ger_reg->SetTitle(list_array[ipt]+": Regression ; mini# ; ppb");
    ger_reg->Draw("AP");
    ger_reg->Fit("pol0");
    
    c1->SaveAs("mini_fit_sorted.pdf");

    ger_lagr->SetMarkerStyle(20);
    ger_lagr->SetTitle(list_array[ipt]+": Lagrange ; mini# ; ppb");
    ger_lagr->Draw("AP");
    ger_lagr->Fit("pol0");
    
    c1->SaveAs("mini_fit_sorted.pdf");
    fnpt[ipt]= ipt;
    fGrandMean_reg[ipt] = ger_reg->GetFunction("pol0")->GetParameter(0);
    fGrandError_reg[ipt] = ger_reg->GetFunction("pol0")->GetParError(0);
    fGrandChisq_reg[ipt] = ger_reg->GetFunction("pol0")->GetChisquare();
    fGrandNDF_reg[ipt] = ger_reg->GetFunction("pol0")->GetNDF();

    fGrandMean_lagr[ipt] = ger_lagr->GetFunction("pol0")->GetParameter(0);
    fGrandError_lagr[ipt] = ger_lagr->GetFunction("pol0")->GetParError(0);
    fGrandChisq_lagr[ipt] = ger_lagr->GetFunction("pol0")->GetChisquare();
    fGrandNDF_lagr[ipt] = ger_lagr->GetFunction("pol0")->GetNDF();

  }
  c1->SaveAs("mini_fit_sorted.pdf]");

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0,1,0.33);
  TPad *pad2 = new TPad("pad2","pad2",0,0.33,1,1);
  pad1->Draw();
  pad2->Draw();
  pad2->SetBottomMargin(0.0);
  pad1->SetTopMargin(0.0);
  pad1->SetBottomMargin(0.15);
  c2->SaveAs("mini_grand_fit_sorted.pdf[");

  TGraphErrors *ger_reg = new TGraphErrors(npt,fnpt,fGrandMean_reg,0,fGrandError_reg);
  TGraphErrors *ger_lagr = new TGraphErrors(npt,fnpt,fGrandMean_lagr,0,fGrandError_lagr);
  TH1D *hsig_reg = new TH1D("hsig_reg","hsig_reg",npt,-0.5,npt-0.5);
  TH1D *hsig_lagr = new TH1D("hsig_lagr","hsig_lagr",npt,-0.5,npt-0.5);
  
  pad2->cd();
  ger_reg->SetTitle("Regression: Mini-run Weighted Average (Upstream Main);;ppb");
  ger_reg->SetMarkerStyle(20);
  ger_reg->Draw("AP");
  ger_reg->Fit("pol0");

  double fPar0_reg = ger_reg->GetFunction("pol0")->GetParameter(0);
  for(int ipt=0;ipt<npt;ipt++){
    hsig_reg->SetBinContent(ipt+1,
			    (fGrandMean_reg[ipt]-fPar0_reg)/fGrandError_reg[ipt]);

    TLatex *data_text = new TLatex(ipt, fGrandMean_reg[ipt], 
				 Form(" %.1f#pm %.1f", fGrandMean_reg[ipt], fGrandError_reg[ipt]));
    // TLatex *data_text2 = new TLatex(ipt, fGrandMean_reg[ipt]-50, 
    // 				 Form(" #chi^{2}/NDF %.1f/ %d", fGrandChisq_reg[ipt], fGrandNDF_reg[ipt]));
    // data_text2->SetTextColor(kBlue);
    data_text->Draw("same");
    // data_text2->Draw("same");
  }
  
  TH1D *hger_reg = (TH1D*)ger_reg->GetHistogram();
  hger_reg->GetXaxis()->Set(npt,-0.5,npt-0.5);
  pad1->cd();
  gStyle->SetOptStat(0);
  hsig_reg->SetTitle(";;Normalized Residual");
  hsig_reg->SetFillColor(kGreen+2);
  hsig_reg->Draw("b");
  for(int ipt=0;ipt<npt;ipt++){
    hsig_reg->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
  }
  hsig_reg->GetYaxis()->SetTitleSize(0.07);
  hsig_reg->GetYaxis()->SetTitleOffset(0.5);
  hsig_reg->GetYaxis()->SetLabelSize(0.07);
  hsig_reg->GetXaxis()->SetLabelSize(0.13);

  c2->SaveAs("mini_grand_fit_sorted.pdf");

  pad1->Clear();
  pad2->Clear();

  pad2->cd();
  ger_lagr->SetTitle("Lagrange Mini-run Weighted Average (Upstream Main);;ppb");
  ger_lagr->SetMarkerStyle(20);
  ger_lagr->Draw("AP");
  ger_lagr->Fit("pol0");

  double fPar0_lagr = ger_lagr->GetFunction("pol0")->GetParameter(0);
  for(int ipt=0;ipt<npt;ipt++){
    hsig_lagr->SetBinContent(ipt+1,
			    (fGrandMean_lagr[ipt]-fPar0_lagr)/fGrandError_lagr[ipt]);

    TLatex *data_text = new TLatex(ipt, fGrandMean_lagr[ipt], 
				 Form(" %.1f#pm %.1f", fGrandMean_lagr[ipt], fGrandError_lagr[ipt]));
    data_text->Draw("same");
  }
  
  TH1D *hger_lagr = (TH1D*)ger_lagr->GetHistogram();
  hger_lagr->GetXaxis()->Set(npt,-0.5,npt-0.5);
  pad1->cd();
  gStyle->SetOptStat(0);
  hsig_lagr->SetTitle(";;Normalized Residual");
  hsig_lagr->SetFillColor(kGreen+2);
  hsig_lagr->Draw("b");
  for(int ipt=0;ipt<npt;ipt++){
    hsig_lagr->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
  }
  hsig_lagr->GetYaxis()->SetTitleSize(0.07);
  hsig_lagr->GetYaxis()->SetTitleOffset(0.5);
  hsig_lagr->GetYaxis()->SetLabelSize(0.07);
  hsig_lagr->GetXaxis()->SetLabelSize(0.13);
  c2->SaveAs("mini_grand_fit_sorted.pdf");

  c2->SaveAs("mini_grand_fit_sorted.pdf]");

}


vector<Int_t> LoadRunList(TString key){
  FILE *runlist = fopen("./list/"+key,"r");
  vector<Int_t> fRet;
  while(!feof(runlist)){
    int run_number=0;
    fscanf(runlist,"%d\n",&run_number);
    if(run_number!=0)
      fRet.push_back(run_number);
  }
  fclose(runlist);
  return fRet;
}

void CompareCentralValues(){

  map<TString,TString> fLinks={ {"reg_5bpm","prex_grand_average_regress_cbpm.root"},
				{"dit","prex_grand_average_dither_respin2.root"},
				{"reg_6bpm","prex_grand_average_reg_6bpm.root"},
				{"lagr_6bpm","prex_grand_average_lagr_6bpm.root"},
				{"reg_all","prex_grand_average_reg_allbpm.root"},
				{"lagr_all","prex_grand_average_lagr_allbpm.root"},
				{"lagr_all_trunc","prex_grand_average_lagr_allbpm_truncated.root"}};

  
  vector< pair<TString,TString> > fSets ={ {"dit","reg_5bpm"},
					   {"dit","lagr_6bpm"},
					   {"dit","lagr_all"},
					   {"dit","reg_all"},
					   {"reg_5bpm","reg_6bpm"},
					   {"reg_5bpm","reg_all"},
					   {"lagr_6bpm","reg_6bpm"},
					   {"lagr_all","reg_all"},
					   {"lagr_6bpm","lagr_all"},
					   {"lagr_all_trunc","lagr_all"}  };

  vector<Double_t> fDeltaMu_slug;
  vector<Double_t> fSigmaDeltaMu_slug;
  vector<Double_t> fChiSquare_slug;
  vector<Int_t> fNDF_slug;
  vector<Double_t> fProb_slug;
  
  vector<Double_t> fDeltaMu_pitt;
  vector<Double_t> fSigmaDeltaMu_pitt;
  vector<Double_t> fChiSquare_pitt;
  vector<Int_t> fNDF_pitt;
  vector<Double_t> fProb_pitt;
  
  vector<Double_t> fDeltaMu_slugfix;
  vector<Double_t> fSigmaDeltaMu_slugfix;
  vector<Double_t> fChiSquare_slugfix;
  vector<Int_t> fNDF_slugfix;
  vector<Double_t> fProb_slugfix;

  vector<Double_t> fDeltaMu_pittfix;
  vector<Double_t> fSigmaDeltaMu_pittfix;
  vector<Double_t> fChiSquare_pittfix;
  vector<Int_t> fNDF_pittfix;
  vector<Double_t> fProb_pittfix;
  
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  TPad *pad1 = new TPad("","",0,0,0.66,1.0);
  TPad *pad2 = new TPad("","",0.66,0,1.0,1.0);
  pad1->Draw();
  pad2->Draw();
  pad1->SetGridx();
  pad1->SetGridy();

  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->Divide(1,2);
  auto iter= fSets.begin();

  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;

    cout << set1 <<": " << fLinks[set1] << endl;
    cout << set2 <<": " << fLinks[set2] << endl;
    
    TFile* file1 = TFile::Open(fLinks[set1]);
    TTree *slug_tree1 = (TTree*)file1->Get("slug");
    slug_tree1->AddFriend(set2+"=slug",fLinks[set2]);
    TTree *pitt_tree1 = (TTree*)file1->Get("pitt");
    pitt_tree1->AddFriend(set2+"=pitt",fLinks[set2]);

    TF1 *ffix = new TF1("ffix","[0]",-100,100);
    ffix->FixParameter(0,0);
    TF1 *ffloat = new TF1("ffloat","[0]",-100,100);
    c1->cd();
    slug_tree1->Draw(Form("sign*(Adet.mean-%s.Adet.mean)*1e9:sqrt(abs(Adet.error**2-%s.Adet.error**2))*1e9:slug",set2.Data(),set2.Data()),"","goff");

    pad1->cd();
    TGraphErrors *ger_slug = new TGraphErrors(slug_tree1->GetSelectedRows(),
					      slug_tree1->GetV3(),slug_tree1->GetV1(),
					      0,slug_tree1->GetV2());
    ger_slug->SetMarkerStyle(20);
    ger_slug->Draw("AP");
    ger_slug->SetTitle(set1 + " vs "+set2 + ":      sign corrected (#mu_{1} - #mu_{2}) #pm #sqrt{#sigma_{1}^{2} - #sigma_{2}^{2} }; Slug ; (ppb) "  );
    ger_slug->GetYaxis()->SetTitleSize(0.05);
    ger_slug->GetXaxis()->SetTitleSize(0.05);
    ger_slug->Fit("ffloat","Q");
    ger_slug->GetXaxis()->SetTitleSize(0.07);
    ger_slug->GetXaxis()->SetTitleOffset(0.5);
    ger_slug->GetYaxis()->SetTitleSize(0.07);
    ger_slug->GetYaxis()->SetTitleOffset(0.5);
    pad2->cd();
    TH1D *hpull_free = new TH1D("hpull_free","Pull Fit",40,-5,8);
    double* yval1_free = slug_tree1->GetV1();
    double* yerr1_free = slug_tree1->GetV2();
    double ymean1_free = ffloat->GetParameter(0);
    for(int ipt=0;ipt<slug_tree1->GetSelectedRows();ipt++){
      hpull_free->Fill( (yval1_free[ipt] -ymean1_free)/yerr1_free[ipt]);
    }
    hpull_free->Draw("");
    hpull_free->Fit("gaus","Q");
    hpull_free->SetName("Pull");
    pad2->Update();
    TPaveStats *pspull_free = (TPaveStats*)hpull_free->FindObject("stats");
    pspull_free->SetOptStat(111111);
    pspull_free->SetY1NDC(0.5);
    pspull_free->SetX1NDC(0.5);

    c1->SaveAs(Form("compare_%s_vs_%s_by_slug_freepar.png",set1.Data(),set2.Data()) );

    fDeltaMu_slug.push_back(ffloat->GetParameter(0));
    fSigmaDeltaMu_slug.push_back(ffloat->GetParError(0));
    fChiSquare_slug.push_back(ffloat->GetChisquare());
    fNDF_slug.push_back(ffloat->GetNDF());
    fProb_slug.push_back(ffloat->GetProb());


    pad1->cd();
    ger_slug->Fit("ffix","Q");
    TH1D *hpull_fixed = new TH1D("hpull_fixed","Pull Fit",40,-5,8);
    for(int ipt=0;ipt<slug_tree1->GetSelectedRows();ipt++){
      hpull_fixed->Fill( (yval1_free[ipt])/yerr1_free[ipt]);
    }
    pad2->cd();
    hpull_fixed->Draw();
    hpull_fixed->Fit("gaus","Q");
    hpull_fixed->SetName("Pull");
    pad2->Update();
    TPaveStats *pspull_fixed = (TPaveStats*)hpull_fixed->FindObject("stats");
    pspull_fixed->SetOptStat(111111);
    pspull_fixed->SetY1NDC(0.5);
    pspull_fixed->SetX1NDC(0.5);
    
    fDeltaMu_slugfix.push_back(ffix->GetParameter(0));
    fSigmaDeltaMu_slugfix.push_back(ffix->GetParError(0));
    fChiSquare_slugfix.push_back(ffix->GetChisquare());
    fNDF_slugfix.push_back(ffix->GetNDF());
    fProb_slugfix.push_back(ffix->GetProb());

    
    c1->SaveAs(Form("compare_%s_vs_%s_by_slug_fixpar.png",set1.Data(),set2.Data()) );

    c2->cd(1);
    gPad->SetGridx();
    slug_tree1->Draw("Adet.rms*1e6:slug","","goff");
    
    TGraph *grms1 = new TGraph(slug_tree1->GetSelectedRows(),
			       slug_tree1->GetV2(),slug_tree1->GetV1());
    grms1->SetMarkerStyle(20);
    grms1->SetMarkerColor(kBlue);
    grms1->SetLineColor(kBlue);

    slug_tree1->Draw(Form("%s.Adet.rms*1e6:slug",set2.Data()),"","goff");
    TGraph *grms2 = new TGraph(slug_tree1->GetSelectedRows(),
			       slug_tree1->GetV2(),slug_tree1->GetV1());
    grms2->SetMarkerStyle(20);
    grms2->SetMarkerColor(kRed);
    grms2->SetLineColor(kRed);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(grms1,"lp");
    mg->Add(grms2,"lp");
    mg->Draw("A");
    mg->SetTitle(set1 + " vs "+set2 + ":  RMS(ppm); Slug ; (ppm) "  );
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetTitleSize(0.05);
    TLegend* legrms = new TLegend(0.9,0.7,0.99,0.9);
    legrms->AddEntry(grms1,set1,"lp");
    legrms->AddEntry(grms2,set2,"lp");
    legrms->Draw("same");
    c2->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    slug_tree1->Draw(Form("sqrt(abs(Adet.rms**2-%s.Adet.rms**2))*1e6:slug",set2.Data()),"","goff");
    // slug_tree1->Draw(Form("(Adet.rms-%s.Adet.rms)*1e6:slug",set2.Data()),"","goff");
    TGraph *grms_diff = new TGraph(slug_tree1->GetSelectedRows(),
				   slug_tree1->GetV2(),slug_tree1->GetV1());

    grms_diff->SetMarkerStyle(20);
    grms_diff->Draw("ALP");
    grms_diff->SetTitle(set1 + " vs "+set2 + ": #sqrt{RMS_{1}^{2} - RMS_{2}^{2} }; Slug ; (ppm) "  );
    grms_diff->GetYaxis()->SetTitleSize(0.05);
    grms_diff->GetXaxis()->SetTitleSize(0.05);
    c2->SaveAs(Form("rms_%s_vs_%s_by_slug.png",set1.Data(),set2.Data()) );
    
    /////////////////// PITT 
    pad1->cd();
    pitt_tree1->Draw(Form("(Adet.mean-%s.Adet.mean)*1e9:sqrt(abs(Adet.error**2-%s.Adet.error**2))*1e9:pitt",set2.Data(),set2.Data()),"","goff");

    TGraphErrors *ger_pitt = new TGraphErrors(pitt_tree1->GetSelectedRows(),
					      pitt_tree1->GetV3(),pitt_tree1->GetV1(),
					      0,pitt_tree1->GetV2());
    ger_pitt->SetMarkerStyle(20);
    ger_pitt->Draw("AP");
    ger_pitt->SetTitle(set1 + " vs "+set2 + ":       (#mu_{1} - #mu_{2}) #pm #sqrt{#sigma_{1}^{2} - #sigma_{2}^{2} }; Pitt ; (ppb) "  );
    ger_pitt->GetYaxis()->SetTitleSize(0.05);
    ger_pitt->GetXaxis()->SetTitleSize(0.05);
    ger_pitt->Fit("ffloat","Q");
    double *yval_pitt = pitt_tree1->GetV1();
    double *yerr_pitt = pitt_tree1->GetV2();
    double pitt_mean = ffloat->GetParameter(0);
    TH1D *hpitt_free = new TH1D("hpitt_free","Pull Fit",39,-5,8);
    TH1D *hpitt_fixed = new TH1D("hpitt_fixed","Pull Fit",39,-5,8);
    
    for(int ipt=0;ipt<pitt_tree1->GetSelectedRows();ipt++){
      hpitt_free->Fill( (yval_pitt[ipt] - pitt_mean)/yerr_pitt[ipt] );
      hpitt_fixed->Fill( (yval_pitt[ipt])/yerr_pitt[ipt] );
    }
    fDeltaMu_pitt.push_back(ffloat->GetParameter(0));
    fSigmaDeltaMu_pitt.push_back(ffloat->GetParError(0));
    fChiSquare_pitt.push_back(ffloat->GetChisquare());
    fNDF_pitt.push_back(ffloat->GetNDF());
    fProb_pitt.push_back(ffloat->GetProb());
    pad2->cd();
    hpitt_free->Draw();
    hpitt_free->Fit("gaus","Q");
    hpitt_free->SetName("Pull");
    pad2->Update();
    TPaveStats *pspitt_free = (TPaveStats*)hpitt_free->FindObject("stats");
    pspitt_free->SetOptStat(111111);
    pspitt_free->SetY1NDC(0.5);
    pspitt_free->SetX1NDC(0.5);

    c1->SaveAs(Form("compare_%s_vs_%s_by_pitt_freepar.png",set1.Data(),set2.Data()) );
    pad1->cd();
    ger_pitt->Fit("ffix","Q");
    fDeltaMu_pittfix.push_back(ffix->GetParameter(0));
    fSigmaDeltaMu_pittfix.push_back(ffix->GetParError(0));
    fChiSquare_pittfix.push_back(ffix->GetChisquare());
    fNDF_pittfix.push_back(ffix->GetNDF());
    fProb_pittfix.push_back(ffix->GetProb());

    pad2->cd();
    hpitt_fixed->Draw();
    hpitt_fixed->Fit("gaus","Q");
    hpitt_fixed->SetName("Pull");
    pad2->Update();
    TPaveStats *pspitt_fixed = (TPaveStats*)hpitt_fixed->FindObject("stats");
    pspitt_fixed->SetOptStat(111111);
    pspitt_fixed->SetY1NDC(0.5);
    pspitt_fixed->SetX1NDC(0.5);
    
    c1->SaveAs(Form("compare_%s_vs_%s_by_pitt_fixpar.png",set1.Data(),set2.Data()) );

    file1->Close();
    
    iter++;
  }

  cout << " ********* Slug Summary ************ " << endl;
  Int_t icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");

    cout << set1 <<  " vs " << set2 << "&";
    printf("%.1f & %.1f & %.1f & %d  & %.2f\\\\ \n",
	   fDeltaMu_slug[icount],fSigmaDeltaMu_slug[icount],
	   fChiSquare_slug[icount],fNDF_slug[icount],fProb_slug[icount]);
    iter++;
    icount++;
  }
  cout << " ********* fix para " << endl;
  icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");
    
    cout << set1 <<  " vs " << set2 << "&";
    printf("%.0f & %.0f & %.1f & %d & %.2f \\\\ \n",
	   fDeltaMu_slugfix[icount],fSigmaDeltaMu_slugfix[icount],
	   fChiSquare_slugfix[icount],fNDF_slugfix[icount],fProb_slugfix[icount]);
    iter++;
    icount++;
  }

  cout << " ********* Pitt Summary ************ " << endl;
  icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");

    cout << set1 <<  " vs " << set2 << "&";
    printf("%.1f & %.1f & %.1f & %d & %.2f \\\\ \n",
	   fDeltaMu_pitt[icount],fSigmaDeltaMu_pitt[icount],
	   fChiSquare_pitt[icount],fNDF_pitt[icount],fProb_pitt[icount]);
    iter++;
    icount++;
  }
  cout << " ********* fix para " << endl;
  icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");

    cout << set1 <<  " vs " << set2 << "&";
    printf("%.0f & %.0f & %.1f & %d & %.2f\\\\ \n",
	   fDeltaMu_pittfix[icount],fSigmaDeltaMu_pittfix[icount],
	   fChiSquare_pittfix[icount],fNDF_pittfix[icount],fProb_pittfix[icount]);
    iter++;
    icount++;
  }

}

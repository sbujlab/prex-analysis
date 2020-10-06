void GrandAverages(){
  map<TString,TString> fLinks={{"dithering","dit_plus_bmw_grand_average_signed_weighted.root"},
			       {"regall","regall_plus_bmw_grand_average_signed_weighted.root"},
			       {"lagrall","lagrall_plus_bmw_grand_average_signed_weighted.root"}};


  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  TPad *pad1 = new TPad("","",0,0,0.66,1.0);
  TPad *pad2 = new TPad("","",0.66,0,1.0,1.0);
  pad1->Draw();
  pad2->Draw();
  // pad1->SetGridx();
  // pad1->SetGridy();
  vector<Double_t> fMeanSlug;
  vector<Double_t> fErrorSlug;
  vector<Double_t> fChiSquareSlug;
  vector<Int_t> fNDFSlug;
  vector<Double_t> fProbSlug;

  vector<Double_t> fMeanPitt;
  vector<Double_t> fErrorPitt;
  vector<Double_t> fChiSquarePitt;
  vector<Int_t> fNDFPitt;
  vector<Double_t> fProbPitt;

  auto iter= fLinks.begin();
  while(iter!= fLinks.end()){
    TString tag  = (*iter).first;
    TString file_name  = (*iter).second;

    TFile *input = TFile::Open(file_name);
    TTree* slug_tree = (TTree*)input->Get("slug");
    TTree* pitt_tree = (TTree*)input->Get("pitt");

    slug_tree->Draw("sign*Adet*1e9 : Adet.error*1e9 :slug + arm_flag/3.0 ","","goff");
    double *yval_slug = slug_tree->GetV1();
    double *yerr_slug = slug_tree->GetV2();
    int npt = slug_tree->GetSelectedRows();
    TGraphErrors *ger_slug= new TGraphErrors(npt, slug_tree->GetV3(), slug_tree->GetV1(),
					     0,slug_tree->GetV2());

    pad1->cd();
    ger_slug->SetMarkerStyle(20);
    ger_slug->Draw("AP");
    ger_slug->SetTitle( Form("%s Corrected Asymmetry (ppb)  vs Slug; Slug;(ppb)" ,tag.Data()));
    ger_slug->GetYaxis()->SetTitleSize(0.05);
    ger_slug->GetYaxis()->SetTitleOffset(0.75);
    ger_slug->GetXaxis()->SetTitleSize(0.05);
    ger_slug->GetXaxis()->SetTitleOffset(0.75);
    ger_slug->Fit("pol0","Q");
    TF1 *ffit = ger_slug->GetFunction("pol0");
    fMeanSlug.push_back(ffit->GetParameter(0));
    fErrorSlug.push_back(ffit->GetParError(0));
    fChiSquareSlug.push_back(ffit->GetChisquare());
    fNDFSlug.push_back(ffit->GetNDF());
    fProbSlug.push_back(ffit->GetProb());
    
    double slug_mean = ffit->GetParameter(0);
    pad2->cd();
    TH1D *hslug = new TH1D("hslug","",40,-5,8);
    for(int ipt=0;ipt<npt;ipt++)
      hslug->Fill( (yval_slug[ipt] - slug_mean)/yerr_slug[ipt]);
    hslug->Draw();
    hslug->Fit("gaus","Q");
    hslug->SetName("PullFit");
    pad2->Update();
    TPaveStats *psslug = (TPaveStats*)hslug->FindObject("stats");
    psslug->SetOptStat(111111);
    psslug->SetY1NDC(0.5);
    psslug->SetX1NDC(0.5);
    c1->SaveAs(Form("prex_slug_averages_%s_plus_bmw.pdf",tag.Data()));


    pitt_tree->Draw("Adet*1e9 : Adet.error*1e9 :pitt ","","goff");
    double *yval_pitt = pitt_tree->GetV1();
    double *yerr_pitt = pitt_tree->GetV2();
    npt = pitt_tree->GetSelectedRows();
    TGraphErrors *ger_pitt= new TGraphErrors(npt, pitt_tree->GetV3(), pitt_tree->GetV1(),
					     0,pitt_tree->GetV2());

    pad1->cd();
    ger_pitt->SetMarkerStyle(20);
    ger_pitt->Draw("AP");
    ger_pitt->SetTitle( Form("%s Corrected Asymmetry (ppb)  vs Pitt; Pitt;(ppb)" ,tag.Data()));
    ger_pitt->GetYaxis()->SetTitleSize(0.05);
    ger_pitt->GetYaxis()->SetTitleOffset(0.75);
    ger_pitt->GetXaxis()->SetTitleSize(0.05);
    ger_pitt->GetXaxis()->SetTitleOffset(0.75);
    ger_pitt->Fit("pol0","Q");
    ffit = ger_pitt->GetFunction("pol0");
    fMeanPitt.push_back(ffit->GetParameter(0));
    fErrorPitt.push_back(ffit->GetParError(0));
    fChiSquarePitt.push_back(ffit->GetChisquare());
    fNDFPitt.push_back(ffit->GetNDF());
    fProbPitt.push_back(ffit->GetProb());

    double pitt_mean = ffit->GetParameter(0);
    pad2->cd();
    TH1D *hpitt = new TH1D("hpitt","",39,-5,8);
    for(int ipt=0;ipt<npt;ipt++)
      hpitt->Fill( (yval_pitt[ipt] - pitt_mean)/yerr_pitt[ipt]);
    hpitt->Draw();
    hpitt->Fit("gaus","Q");
    hpitt->SetName("PullFit");
    pad2->Update();
    TPaveStats *pspitt = (TPaveStats*)hpitt->FindObject("stats");
    pspitt->SetOptStat(111111);
    pspitt->SetY1NDC(0.5);
    pspitt->SetX1NDC(0.5);
    c1->SaveAs(Form("prex_pitt_averages_%s_plus_bmw.pdf",tag.Data()));


    // NULL asymmetry
    c1->Clear("D");
    pitt_tree->Draw("null_Adet*1e9 : null_Adet.error*1e9 :pitt ","","goff");
    double *nulval_pitt = pitt_tree->GetV1();
    double *nulerr_pitt = pitt_tree->GetV2();
    npt = pitt_tree->GetSelectedRows();
    TGraphErrors *ger_pitt_null= new TGraphErrors(npt, pitt_tree->GetV3(), pitt_tree->GetV1(),
    					     0,pitt_tree->GetV2());

    pad1->cd();
    ger_pitt_null->SetMarkerStyle(20);
    ger_pitt_null->Draw("AP");
    ger_pitt_null->SetTitle( Form("%s Corrected Null (ppb)  vs Pitt; Pitt;(ppb)" ,tag.Data()));
    ger_pitt_null->GetYaxis()->SetTitleSize(0.05);
    ger_pitt_null->GetYaxis()->SetTitleOffset(0.75);
    ger_pitt_null->GetXaxis()->SetTitleSize(0.05);
    ger_pitt_null->GetXaxis()->SetTitleOffset(0.75);
    ger_pitt_null->Fit("pol0","Q");
    
    double nul_mean = ger_pitt_null->GetFunction("pol0")->GetParameter(0);
    double pitt_null = ffit->GetParameter(0);
    pad2->cd();
    TH1D *hpittnul = new TH1D("hpittnul","",39,-5,8);
    for(int ipt=0;ipt<npt;ipt++)
      hpittnul->Fill( (nulval_pitt[ipt] - nul_mean)/nulerr_pitt[ipt]);
    hpittnul->Draw();
    hpittnul->Fit("gaus","Q");
    hpittnul->SetName("PullFit");
    pad2->Update();
    TPaveStats *pspittnul = (TPaveStats*)hpittnul->FindObject("stats");
    pspittnul->SetOptStat(111111);
    pspittnul->SetY1NDC(0.5);
    pspittnul->SetX1NDC(0.5);
    c1->SaveAs(Form("prex_pitt_null_%s_plus_bmw.pdf",tag.Data()));

    input->Close();
    iter++;
  }
  cout << " ********* Slug Summary ************ " << endl;
  Int_t icount=0;
  iter = fLinks.begin();
  while(iter!=fLinks.end()){
    TString tag = (*iter).first;
    printf("%s & %.1f & %.1f & %.1f & %d & %.2f\\\\ \n",
	   tag.Data(), fMeanSlug[icount],fErrorSlug[icount],
	   fChiSquareSlug[icount],fNDFSlug[icount],fProbSlug[icount]);
    iter++;
    icount++;
  }

  cout << " ********* Pitt Summary ************ " << endl;
  icount=0;
  iter = fLinks.begin();
  while(iter!=fLinks.end()){
    TString tag = (*iter).first;
    printf("%s & %.1f & %.1f & %.1f & %d & %.2f\\\\ \n",
	   tag.Data(), fMeanPitt[icount],fErrorPitt[icount],
	   fChiSquarePitt[icount],fNDFPitt[icount],fProbPitt[icount]);
    iter++;
    icount++;
  }
}

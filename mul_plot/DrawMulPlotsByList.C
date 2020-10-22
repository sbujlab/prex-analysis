void DrawMulPlotsByList(){
  TString list_array[]={"50uA","50uA_1arm",
			"70uA_120Hz","70uA_120Hz_1arm",
			"70uA_240Hz","70uA_240Hz_1arm",
			"85uA"};
  Int_t npt=  sizeof(list_array)/sizeof(*list_array);

  gStyle->SetOptStat("ourMeN");
  gStyle->SetStatX(1);
  gStyle->SetStatY(0.93);
  gStyle->SetStatW(0.2);
  gStyle->SetOptFit(1);


  TCanvas *c1 = new TCanvas("c1","c1",1800,600);
  double *fMean_reg = new double[npt];
  double *fError_reg = new double[npt];
  double *fMean_lagr = new double[npt];
  double *fError_lagr = new double[npt];

  double *fMean_fitted_lagr = new double[npt];
  double *fError_fitted_lagr = new double[npt];

  double *fnpt = new double[npt];
  c1->Divide(3,1);
  c1->cd(1);
  gPad->SetLogy();
  c1->cd(2);
  gPad->SetLogy();
  c1->cd(3);
  gPad->SetLogy();
  c1->SaveAs("mul_plot_sorted.pdf[");
  double total_counts = 0;
  for(int ipt=0;ipt<npt;ipt++){
    c1->Clear("D");
    TFile* input = TFile::Open("mul1d_"+list_array[ipt]+".root");
    TH1D* h1dreg = (TH1D*)input->Get("h1dreg");
    TH1D* h1dlagr = (TH1D*)input->Get("h1dlagr");
    TH1D* h1dbcm = (TH1D*)input->Get("h1dbcm");

    fMean_reg[ipt] = h1dreg->GetMean()*1e3;
    fError_reg[ipt] = h1dreg->GetMeanError()*1e3;
    fMean_lagr[ipt] = h1dlagr->GetMean()*1e3;
    fError_lagr[ipt] = h1dlagr->GetMeanError()*1e3;
    fnpt[ipt] = ipt;

    h1dreg->SetTitle(list_array[ipt]+": Regression (ppm); ppm; counts");
    h1dlagr->SetTitle(list_array[ipt]+": Lagrange  (ppm); ppm; counts");
    h1dbcm->SetTitle("Beam Current (uA) ; uA; counts");

    c1->cd(1);
    h1dreg->Draw();
    h1dreg->Fit("gaus");
    h1dreg->GetFunction("gaus")->SetNpx(1000);
    h1dreg->GetFunction("gaus")->SetLineStyle(2);

    c1->cd(2);
    h1dlagr->Draw();
    h1dlagr->Fit("gaus");
    h1dlagr->GetFunction("gaus")->SetNpx(1000);
    h1dlagr->GetFunction("gaus")->SetLineStyle(2);
    fMean_fitted_lagr[ipt] = h1dlagr->GetFunction("gaus")->GetParameter(1) * 1e3;
    fError_fitted_lagr[ipt] = h1dlagr->GetFunction("gaus")->GetParError(1) * 1e3;

    c1->cd(3);
    h1dbcm->SetFillColor(kAzure+1);
    h1dbcm->Draw("b");
    c1->SaveAs("mul_plot_sorted.pdf");
    total_counts += h1dbcm->GetEntries();
    input->Close();
  }
  c1->SaveAs("mul_plot_sorted.pdf]");
  cout << "+++++++++++++++++++" << endl;
  cout << "Total counts: " << (int)total_counts << endl;
  cout << "+++++++++++++++++++" << endl;

  gStyle->SetOptFit(1111);
  TCanvas *c2 = new TCanvas("c2","c2",1200,400);
  c2->cd();
  // TPad *pad1 = new TPad("pad1","pad1",0,0,1,0.33);
  // TPad *pad2 = new TPad("pad2","pad2",0,0.33,1,1);
  // pad1->Draw();
  // pad2->Draw();
  // pad2->SetBottomMargin(0.0);
  // pad1->SetTopMargin(0.0);
  // pad1->SetBottomMargin(0.15);
  c2->SaveAs("mul_grand_fit_sorted.pdf[");
  TGraphErrors *ger_reg = new TGraphErrors(npt,fnpt,fMean_reg,0,fError_reg);
  TGraphErrors *ger_lagr = new TGraphErrors(npt,fnpt,fMean_lagr,0,fError_lagr);
  TGraphErrors *ger_fitted_lagr = new TGraphErrors(npt,fnpt,fMean_fitted_lagr,0,fError_fitted_lagr);
  TH1D *hsig_reg = new TH1D("hsig_reg","hsig_reg",npt,-0.5,npt-0.5);
  TH1D *hsig_lagr = new TH1D("hsig_lagr","hsig_lagr",npt,-0.5,npt-0.5);
  TH1D *hsig_fitted_lagr = new TH1D("hsig_fitted_lagr","hsig_fitted_lagr",npt,-0.5,npt-0.5);
  // pad2->cd();
  ger_reg->SetTitle("Regression: Multiplet Calculated Mean (Upstream Main);;ppb");
  ger_reg->SetMarkerStyle(20);
  ger_reg->Draw("AP");
  ger_reg->Fit("pol0");
  double fPar0_reg = ger_reg->GetFunction("pol0")->GetParameter(0);
  for(int ipt=0;ipt<npt;ipt++){
    hsig_reg->SetBinContent(ipt+1,
			    (fMean_reg[ipt]-fPar0_reg)/fError_reg[ipt]);

    TLatex *data_text = new TLatex(ipt, fMean_reg[ipt], 
				 Form(" %.1f#pm %.1f", fMean_reg[ipt], fError_reg[ipt]));
    data_text->Draw("same");
  }
  
  TH1D *hger_reg = (TH1D*)ger_reg->GetHistogram();
  hger_reg->GetXaxis()->Set(npt,-0.5,npt-0.5);
  // pad1->cd();
  // gStyle->SetOptStat(0);
  // hsig_reg->SetTitle(";;Normalized Residual");
  // hsig_reg->SetFillColor(kGreen+2);
  // hsig_reg->Draw("b");
  for(int ipt=0;ipt<npt;ipt++){
    // hsig_reg->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
    hger_reg->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
  }
  hger_reg->GetYaxis()->SetTitleSize(0.07);
  hger_reg->GetYaxis()->SetTitleOffset(0.5);
  hger_reg->GetYaxis()->SetLabelSize(0.07);
  hger_reg->GetXaxis()->SetLabelSize(0.05);

  c2->SaveAs("mul_grand_fit_sorted.pdf");
  // Lagrange Mul calculated Mean
  // pad2->Clear();
  // pad1->Clear();
  // pad2->cd();
  c2->Clear();
  ger_lagr->SetTitle("Lagrange: Multiplet Calculated Mean (Upstream Main); ;ppb");
  ger_lagr->SetMarkerStyle(20);
  ger_lagr->Draw("AP");
  ger_lagr->Fit("pol0");
  TH1D *hger_lagr = (TH1D*)ger_lagr->GetHistogram();
  hger_lagr->GetXaxis()->Set(npt,-0.5,npt-0.5);
  // double fPar0_lagr = ger_lagr->GetFunction("pol0")->GetParameter(0);
  for(int ipt=0;ipt<npt;ipt++){
    // hsig_lagr->SetBinContent(ipt+1,
    // 			     (fMean_lagr[ipt]-fPar0_lagr)/fError_lagr[ipt]);

    TLatex *data_text = new TLatex(ipt, fMean_lagr[ipt], 
  				   Form(" %.1f#pm %.1f", fMean_lagr[ipt], fError_lagr[ipt]));
    data_text->Draw("same");

  }
  // pad1->cd();
  // hsig_lagr->SetTitle(";;Normalized Residual");
  // hsig_lagr->SetFillColor(kGreen+2);
  // hsig_lagr->Draw("b");

  for(int ipt=0;ipt<npt;ipt++){
    // hsig_lagr->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
    hger_lagr->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
  }
  hger_lagr->GetYaxis()->SetTitleSize(0.07);
  hger_lagr->GetYaxis()->SetTitleOffset(0.5);
  hger_lagr->GetYaxis()->SetLabelSize(0.07);
  hger_lagr->GetXaxis()->SetLabelSize(0.05);


  c2->SaveAs("mul_grand_fit_sorted.pdf");

  // Lagrange Mul Fitted Mean
  // pad2->Clear();
  // pad1->Clear();
  // pad2->cd();
  c2->Clear();
  ger_fitted_lagr->SetTitle("Lagrange: Multiplet Fitted Mean (Upstream Main); ;ppb");
  ger_fitted_lagr->SetMarkerStyle(20);
  ger_fitted_lagr->Draw("AP");
  ger_fitted_lagr->Fit("pol0");
  TH1D *hger_fitted_lagr = (TH1D*)ger_fitted_lagr->GetHistogram();
  hger_fitted_lagr->GetXaxis()->Set(npt,-0.5,npt-0.5);
  // double fPar0_fitted_lagr = ger_fitted_lagr->GetFunction("pol0")->GetParameter(0);
  for(int ipt=0;ipt<npt;ipt++){
    // hsig_fitted_lagr->SetBinContent(ipt+1,
    // 			     (fMean_fitted_lagr[ipt]-fPar0_fitted_lagr)/fError_fitted_lagr[ipt]);

    TLatex *data_text = new TLatex(ipt, fMean_fitted_lagr[ipt], 
  				   Form(" %.1f#pm %.1f", fMean_fitted_lagr[ipt], fError_fitted_lagr[ipt]));
    data_text->Draw("same");

  }
  // pad1->cd();
  // hsig_fitted_lagr->SetTitle(";;Normalized Residual");
  // hsig_fitted_lagr->SetFillColor(kGreen+2);
  // hsig_fitted_lagr->Draw("b");

  for(int ipt=0;ipt<npt;ipt++){
    // hsig_fitted_lagr->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
    hger_fitted_lagr->GetXaxis()->SetBinLabel(ipt+1,list_array[ipt]);
  }
  hger_fitted_lagr->GetYaxis()->SetTitleSize(0.07);
  hger_fitted_lagr->GetYaxis()->SetTitleOffset(0.5);
  hger_fitted_lagr->GetYaxis()->SetLabelSize(0.07);
  hger_fitted_lagr->GetXaxis()->SetLabelSize(0.05);

  c2->SaveAs("mul_grand_fit_sorted.pdf");
  c2->SaveAs("mul_grand_fit_sorted.pdf]");
}

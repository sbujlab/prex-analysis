void Draw2DHistory(){
  TFile* input = TFile::Open("mul2d_run3897-4864.root");
  TCanvas *c1 = new TCanvas("c1","c1",1200,300);
  gStyle->SetOptStat(0);
  c1->SaveAs("prex_2d_plot_run3897-4864.pdf[");
  TH2D* h2dreg = (TH2D*)input->Get("h2dreg");
  TH2D* h2dlagr = (TH2D*)input->Get("h2dlagr");
  TH2D* h2dregdd = (TH2D*)input->Get("h2dregdd");
  TH2D* h2dlagrdd = (TH2D*)input->Get("h2dlagrdd");

  TH2D* h2dreg_slug = (TH2D*)input->Get("h2dreg_slug");
  TH2D* h2dlagr_slug = (TH2D*)input->Get("h2dlagr_slug");
  TH2D* h2dregdd_slug = (TH2D*)input->Get("h2dregdd_slug");
  TH2D* h2dlagrdd_slug = (TH2D*)input->Get("h2dlagrdd_slug");

  TH2D* hist_array[]={h2dreg,h2dlagr,h2dregdd,h2dlagrdd};
  TH2D* hist_slug[]={h2dreg_slug,h2dlagr_slug,h2dregdd_slug,h2dlagrdd_slug};

  TString Title[]={"Regression Corrected us avg (ppm); run; ppm",
		   "Lagrange Corrected us avg (ppm); run; ppm",
		   "Regression Corrected us dd (ppm); run; ppm",
		   "Lagrange Corrected us dd (ppm); run; ppm"};

  TString Title_slug[]={"Regression Corrected us avg (ppm); slug; ppm",
			"Lagrange Corrected us avg (ppm); slug; ppm",
			"Regression Corrected us dd (ppm); slug; ppm",
			"Lagrange Corrected us dd (ppm); slug; ppm"};
  TCandle::SetBoxRange(0.68); // 1 sigma
  TCandle::SetWhiskerRange(0.9974); // 3 sigma
 
  for(int i=0;i<4;i++){
    cout << "================================= " << endl;
    cout << Title[i] << endl;
    cout << "================================= " << endl;
    Int_t nBinsX = hist_array[i]->GetNbinsX();
    for(int ibinx=1;ibinx<=nBinsX;ibinx++){
      TH1D *projectY = hist_array[i]->ProjectionY("",ibinx,ibinx);
      Int_t nBinsY = projectY->GetNbinsX();
      double BinWidth = projectY->GetBinWidth(1);
      Int_t nWindow = 800/BinWidth;
      Bool_t kFound = kFALSE;
      for(int ibiny=1;ibiny<=nBinsY/2-nWindow;ibiny++){
	if(projectY->GetBinContent(ibiny)>0){
	  cout << "-- Run " << hist_array[i]->GetXaxis()->GetBinCenter(ibinx) 
	       << " has outlier ("
	       << projectY->GetBinCenter(ibiny) 
	       << " ppm)" << endl;
	  kFound = kTRUE;
	  break;
	}
      }
      if(kFound)
	continue;

      for(int ibiny=nBinsY;ibiny>=nBinsY/2+nWindow;ibiny--){
	if(projectY->GetBinContent(ibiny)>0){
	  cout << "-- Run " << hist_array[i]->GetXaxis()->GetBinCenter(ibinx) 
	       << " has outlier ("
	       << projectY->GetBinCenter(ibiny) 
	       << " ppm)" << endl;
	  kFound = kTRUE;
	  break;
	}
      }
    }
  }

  for(int i=0;i<4;i++){
    hist_slug[i]->SetTitle(Title_slug[i]);
    hist_slug[i]->SetMarkerSize(0.5);
    hist_slug[i]->Draw("candlex(10111111)"); //zhpawMmb
    c1->SaveAs("prex_2d_plot_run3897-4864.pdf");
  }



  c1->SaveAs("prex_2d_plot_run3897-4864.pdf]");
  input->Close();
}

void DrawMulPlots(){
  gStyle->SetOptStat("ourMeN");
  gStyle->SetStatX(1);
  gStyle->SetStatY(0.93);
  gStyle->SetStatW(0.2);
  gStyle->SetOptFit(1);

  TFile* input = TFile::Open("mulhist_run3897-4864.root");
  TCanvas *c1 = new TCanvas("c1","c1",1300,600);
  c1->Divide(2,1);
  TH1D* h1dreg = (TH1D*)input->Get("h1dreg");
  TH1D* h1dlagr = (TH1D*)input->Get("h1dlagr");
  TH1D* h1dregdd = (TH1D*)input->Get("h1dregdd");
  TH1D* h1dlagrdd = (TH1D*)input->Get("h1dlagrdd");

  TH1D* hist_array[]={h1dreg,h1dlagr,h1dregdd,h1dlagrdd};
  TString Title[]={"Regression Corrected us avg (ppm); ppm; counts",
		   "Lagrange Corrected us avg (ppm); ppm; counts",
		   "Regression Corrected us dd (ppm); ppm; counts",
		   "Lagrange Corrected us dd (ppm); ppm; counts"};
  c1->SaveAs("prex_mul_plot_run3897-4864.pdf[");
  for(int i=0;i<4;i++){
    hist_array[i]->SetTitle(Title[i]);
    c1->cd(1);
    hist_array[i]->Draw();
    hist_array[i]->Fit("gaus");
    hist_array[i]->GetFunction("gaus")->SetNpx(1000);
    hist_array[i]->GetFunction("gaus")->SetLineStyle(2);
    c1->cd(2);
    gPad->SetLogy();
    hist_array[i]->Draw();
    // hist_array[i]->Fit("gaus");
    c1->SaveAs("prex_mul_plot_run3897-4864.pdf");
  }
  c1->SaveAs("prex_mul_plot_run3897-4864.pdf]");
  input->Close();
}

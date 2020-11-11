void DrawMulPlots_pro(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TFile* input = TFile::Open("mulhist_run3897-4864.root");
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogy();
  TH1D* h1dlagr = (TH1D*)input->Get("h1dlagr");

  TH1D* hist_array[]={h1dlagr};
  TString Title[]={";Multiplet asymmetry(ppm); Number of multiplets"};

  c1->SaveAs("mul_plot_run3897-4864_pro.pdf[");
  for(int i=0;i<1;i++){
    c1->cd();
    hist_array[i]->SetTitleFont(22);
    hist_array[i]->SetLabelFont(22);
    hist_array[i]->GetYaxis()->SetTitleFont(22);
    hist_array[i]->GetYaxis()->SetLabelFont(22);

    hist_array[i]->SetTitle(Title[i]);

    hist_array[i]->SetLineWidth(3);
    hist_array[i]->Draw();
    hist_array[i]->Fit("gaus");
    hist_array[i]->GetFunction("gaus")->SetNpx(1000);
    hist_array[i]->GetFunction("gaus")->SetLineStyle(2);
    hist_array[i]->GetFunction("gaus")->SetLineWidth(3);
    TLatex *Entries = new TLatex(650,1e5,Form("%.1f#times 10^{6} entries",h1dlagr->GetEntries()/1e6));
    TLatex *sigma = new TLatex(1000,2e4,Form("#sigma=%.1f ppm",h1dlagr->GetRMS()));
    Entries->SetTextFont(22);
    sigma->SetTextFont(22);

    Entries->Draw("same");
    sigma->Draw("same");
    c1->SaveAs("mul_plot_run3897-4864_pro.pdf");
  }
  c1->SaveAs("mul_plot_run3897-4864_pro.pdf]");
  input->Close();
}

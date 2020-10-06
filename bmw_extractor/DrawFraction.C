void DrawFraction(){
  TFile *plus_input = TFile::Open("regall_plus_bmw_grand_average_signed_weighted.root");
  TFile *bmw_input = TFile::Open("regall_bmw_only_grand_average_signed_weighted.root");

  TTree *plus_tree = (TTree*)plus_input->Get("slug");
  TTree *bmw_tree = (TTree*)bmw_input->Get("slug");

  plus_tree->AddFriend(bmw_tree,"bmw");
  Int_t npt=plus_tree->Draw("Adet.nsamp:bmw.Adet.nsamp:slug:arm_flag","","goff");
  double* v1 = plus_tree->GetV1();
  double* v2 = plus_tree->GetV2();
  double* v3 = plus_tree->GetV3();
  double* v4 = plus_tree->GetV4();
  TH1D *h1d = new TH1D("h1d","BMW data percentage Fraction vs Slug",npt,-0.5,npt-0.5);
  TH1D *h1draw = new TH1D("h1draw","Number of patterns vs Slug",npt,-0.5,npt-0.5);
  for(int ibin=0;ibin<npt;ibin++){
    TString arm_flag ="";
    if(v4[ibin]==1)
      arm_flag = "R";
    if(v4[ibin]==2)
      arm_flag = "L";
    h1d->GetXaxis()->SetBinLabel(ibin+1, Form("%d %s",(int)v3[ibin], arm_flag.Data()));
    h1d->SetBinContent(ibin+1, 100*(v2[ibin]/v1[ibin]));
    h1draw->GetXaxis()->SetBinLabel(ibin+1, Form("%d %s",(int)v3[ibin], arm_flag.Data()));
    h1draw->SetBinContent(ibin+1, v2[ibin]);
  }
  TCanvas* c1 = new TCanvas("c1","c1",1200,400);
  c1->cd();
  h1d->SetFillColor(38);
  h1d->GetYaxis()->SetTitle(" Fraction (%)");
  h1d->GetXaxis()->SetTitle(" Slug ");
  h1d->Draw("bar");

  TCanvas* c2 = new TCanvas("c2","c2",1200,400);
  c2->cd();
  h1draw->SetFillColor(38);
  h1draw->GetYaxis()->SetTitle(" counts");
  h1draw->GetXaxis()->SetTitle(" Slug ");
  h1draw->Draw("bar");
}

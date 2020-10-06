void bmwOnlyTestBySlug(){

  TFile* bmw_input = TFile::Open("regall_plus_bmw_grand_average_signed_weighted.root");
  TFile* prod_input = TFile::Open("regall_prod_only_grand_average_signed_weighted.root");

  // TFile* bmw_input = TFile::Open("lagrall_plus_bmw_grand_average_signed_weighted.root");
  // TFile* prod_input = TFile::Open("lagrall_prod_only_grand_average_signed_weighted.root");

  // TFile* bmw_input = TFile::Open("dit_plus_bmw_grand_average_signed_weighted.root");
  // TFile* prod_input = TFile::Open("dit_prod_only_grand_average_signed_weighted.root");


  TTree *prod_tree = (TTree*)prod_input->Get("slug");
  TTree *bmw_tree = (TTree*)bmw_input->Get("slug");
  prod_tree->AddFriend(bmw_tree,"bmw");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->cd();
  TPad *pad1 = new TPad("","",0,0,0.66,1.0);
  TPad *pad2 = new TPad("","",0.66,0,1.0,1.0);
  pad1->Draw();
  pad2->Draw();

  prod_tree->SetMarkerStyle(8);
  gStyle->SetOptFit(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  pad1->cd();
  prod_tree->Draw("sign*(Adet.mean - bmw.Adet.mean)*1e9: 1e9*sqrt(pow(Adet.error,2) -pow(bmw.Adet.error,2)):slug","Adet.mean - bmw.Adet.mean!=0","goff");


  TGraphErrors* mger = new TGraphErrors(prod_tree->GetSelectedRows(),
					prod_tree->GetV3(), prod_tree->GetV1(),
					0,prod_tree->GetV2());
  mger->SetMarkerStyle(20);
  mger->Draw("AP");
  TF1 *f1 = new TF1("f1","pol0",-1e6,1e6);
  mger->Fit("f1","Q");
  mger->GetXaxis()->SetTitle("Run Number");
  mger->GetYaxis()->SetTitle(" #Delta (ppb)");
  mger->GetYaxis()->SetTitleSize(0.05);
  mger->SetTitle( "Regression : #Delta =A_{prod}-A_{prod+bmw}, #sigma_{#Delta} = #sqrt{ #sigma_{prod}^{2} - #sigma_{prod+bmw}^{2} }");
  double grand_mean = f1->GetParameter(0);

  pad2->cd();
  prod_tree->Draw(Form("(sign*(Adet.mean - bmw.Adet.mean)*1e9 -%f )/ (1e9*sqrt( pow(Adet.error,2) -pow(bmw.Adet.error,2)))", grand_mean),"Adet.mean - bmw.Adet.mean!=0");
  TH1D *h1dptr = (TH1D*)gPad->FindObject("htemp");
  h1dptr->Fit("gaus");
  h1dptr->SetTitle("Pull Fit");
}

void bmwOnlyTest(){
  TFile* input = TFile::Open("bmw_plus_allslugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini");
  TTree *bmw_tree = (TTree*)input->Get("mini_bmw");
  TTree *info_tree = (TTree*)input->Get("mini_info");
  mini_tree->AddFriend(info_tree);
  mini_tree->AddFriend(bmw_tree);
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->cd();
  TPad *pad1 = new TPad("","",0,0,0.66,1.0);
  TPad *pad2 = new TPad("","",0.66,0,1.0,1.0);
  pad1->Draw();
  pad2->Draw();

  mini_tree->SetMarkerStyle(8);
  gStyle->SetOptFit(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  pad1->cd();
  mini_tree->Draw("spin*(reg_asym_us_avg.mean - mini_bmw.reg_asym_us_avg.mean)*1e9: 1e9*sqrt( pow(reg_asym_us_avg.err,2) +pow(mini_bmw.reg_asym_us_avg.err,2)):run+mini/10.0",
		  "arm_flag==0 && kGood==1 && mini_bmw.reg_asym_us_avg.err>0","goff");
    
  TGraphErrors* mger = new TGraphErrors(mini_tree->GetSelectedRows(),
					mini_tree->GetV3(), mini_tree->GetV1(),
					0,mini_tree->GetV2());
  mger->SetMarkerStyle(20);
  mger->Draw("AP");
  TF1 *f1 = new TF1("f1","pol0",-1e6,1e6);
  mger->Fit("f1","Q");
  mger->GetXaxis()->SetTitle("Run Number");
  mger->GetYaxis()->SetTitle(" #Delta (ppb)");
  mger->GetYaxis()->SetTitleSize(0.05);
  mger->SetTitle( "Regression : #Delta =A_{prod}-A_{bmw}, #sigma_{#Delta} = #sqrt{ #sigma_{prod}^{2} + #sigma_{bmw}^{2} }");
  double grand_mean = f1->GetParameter(0);

  pad2->cd();
  mini_tree->Draw(Form("(spin*(reg_asym_us_avg.mean - mini_bmw.reg_asym_us_avg.mean)*1e9 -%f )/ (1e9*sqrt( pow(reg_asym_us_avg.err,2) +pow(mini_bmw.reg_asym_us_avg.err,2)))", grand_mean),
		  "arm_flag==0 && kGood==1 && mini_bmw.reg_asym_us_avg.err>0","");
  TH1D *h1dptr = (TH1D*)gPad->FindObject("htemp");
  h1dptr->Fit("gaus");
  h1dptr->SetTitle("Pull Fit");
}

void DrawCentralValues(){
  
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);
  TCanvas *c1 = new TCanvas("c1","c1",1400,600);
  // c1->SetGridx();
  // c1->SetGridy();
  c1->cd();
  TPad *pad1 = new TPad("","",0,0,0.66,1.0);
  TPad *pad2 = new TPad("","",0.66,0,1.0,1.0);
  pad1->Draw();
  pad2->Draw();
  pad1->SetGridx();
  pad1->SetGridy();
  TF1 *fp0free = new TF1("fp0free","[0]",-100,100);

  TFile *input =  TFile::Open(
  TTree *slug_tree = (TTree*)input->Get("slug");

  pad1->cd();
  TMultiGraph *mg_er = new TMultiGraph();
  TLegend *leg_er = new TLegend(0.99,0.9,0.9,0.7);
  slug_tree->Draw("sign*Adet:Adet.err:slug+arm_flag/6","","goff");
  // pretending it is both-only...
  TGraphErrors* ger_both = new TGraphErrors(slug_tree->GetSelectedRows(),
					    slug_tree->GetV3(),slug_tree->GetV1(),
					    0,slug_tree->GetV2());
    
  ger_both->Fit("fp0free","Q");  
  double* yval = slug_tree->GetV1();
  double* yerr = slug_tree->GetV2();
  int ynpt = slug_tree->GetSelectedRows();
  TH1D *hpull_free = new TH1D("hpull_free","Pull Fit",40,-5,8);
  TH1D *hpull_fixed = new TH1D("hpull_fixed","Pull Fit",40,-5,8);
  double mean = fp0free->GetParameter(0);
  for(int i=0;i<ynpt;i++){
    hpull_free->Fill((yval[i]-mean)/yerr[i]);
    hpull_fixed->Fill(yval[i]/yerr[i]);
  }
  hpull_free->Fit("gaus","Q");
  hpull_fixed->Fit("gaus","Q");
  hpull_free->SetName("Pull");
  hpull_fixed->SetName("Pull");
    
  slug_tree->Draw("sign*Adet:Adet.err:slug+arm_flag/6","arm_flag==1","goff");
  TGraphErrors* ger_right = new TGraphErrors(slug_tree->GetSelectedRows(),
					     slug_tree->GetV3(),slug_tree->GetV1(),
					     0,slug_tree->GetV2());
  slug_tree->Draw("sign*Adet:Adet.err:slug+arm_flag/6","arm_flag==2","goff");
  TGraphErrors* ger_left = new TGraphErrors(slug_tree->GetSelectedRows(),
					    slug_tree->GetV3(),slug_tree->GetV1(),
					    0,slug_tree->GetV2());
    
  ger_both->SetLineColor(kBlack);
  ger_both->SetMarkerColor(kBlack);
  ger_both->SetMarkerStyle(20);
  ger_right->SetLineColor(kBlue);
  ger_right->SetMarkerColor(kBlue);
  ger_right->SetMarkerStyle(20);
  ger_left->SetLineColor(kRed);
  ger_left->SetMarkerColor(kRed);
  ger_left->SetMarkerStyle(20);
    
  mg_er->Add(ger_both,"p");
  mg_er->Add(ger_right,"p");
  mg_er->Add(ger_left,"p");
  leg_er->AddEntry(ger_both,"both-arm","p");
  leg_er->AddEntry(ger_right,"right-arm","p");
  leg_er->AddEntry(ger_left,"left-arm","p");
  pad1->cd();
  mg_er->Draw("A");
  // ger_both->Fit("fp0free","Q");  
  leg_er->Draw("same");
  mg_er->SetTitle(Form(" %s vs %s :Sign Corrected  Apv (ppb) vs Slug; Slug; #Apv (ppb) ",
		       set1.Data(),set2.Data()));
  mg_er->GetXaxis()->SetTitleSize(0.07);
  mg_er->GetXaxis()->SetTitleOffset(0.5);
  mg_er->GetYaxis()->SetTitleSize(0.07);
  mg_er->GetYaxis()->SetTitleOffset(0.5);
  double ymax = mg_er->GetYaxis()->GetXmax();
  double ymin = mg_er->GetYaxis()->GetXmin();
  mg_er->GetYaxis()->SetRangeUser(ymin,ymax+0.3*(ymax-ymin));
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)ger_both->FindObject("stats");
  ps->SetX2NDC(0.9);
  pad2->cd();
  hpull_free->Draw();
  pad2->Update();
  TPaveStats *pspull_free = (TPaveStats*)hpull_free->FindObject("stats");
  pspull_free->SetOptStat(111111);
  pspull_free->SetY1NDC(0.5);
  pspull_free->SetX1NDC(0.5);
  fDeltaMu_slug.push_back(fp0free->GetParameter(0));
  fSigmaDeltaMu_slug.push_back(fp0free->GetParError(0));
  fChiSquare_slug.push_back(fp0free->GetChisquare());
  fProb_slug.push_back(fp0free->GetProb());
  fNDF_slug.push_back(fp0free->GetNDF());
  c1->SaveAs(Form("./plots/central_val_Apv_%s_%s_by_slug.pdf",set1.Data(),set2.Data()));
    
  input->Close();

}

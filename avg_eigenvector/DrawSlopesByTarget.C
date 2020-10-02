void DrawSlopesByTarget(){
  TFile* input = TFile::Open("rootfiles/slug3to94_sorted_evector.root");
  TTree *eig_tree = (TTree*)input->Get("eig");
  TTree *lagr_tree = (TTree*)input->Get("lagr");
  TTree *reg_tree = (TTree*)input->Get("reg");
  TTree *info_tree = (TTree*)input->Get("mini_info");
  eig_tree->AddFriend(info_tree);
  eig_tree->AddFriend(lagr_tree);
  eig_tree->AddFriend(reg_tree);
  gStyle->SetTitleFontSize(0.07);
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(1,3);
  c1->cd(1);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  c1->cd(2);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  c1->cd(3);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);

  eig_tree->SetMarkerStyle(8);
  TString title[] = {"208Pb2", "208Pb10", "208Pb9", "208Pb8"};
  TString slug_cut[] = {" && slug>=12 && slug <=13 ",
			"&& slug>=26 && slug<=28",
			"&& slug>=44 && slug<=45",
			"&& slug>=58 && slug<=59"};
  for(int ic=0;ic<4;ic++){
    c1->Clear("D");
    
    Int_t npt=eig_tree->Draw("run:Entry$","arm_flag==0 && kGood" +slug_cut[ic],"goff");
    double *fRun_ptr = eig_tree->GetV1();
    double *fEntry_ptr = eig_tree->GetV2();
    double *fRun = new double[npt];
    double *fEntry = new double[npt];
    for(int ipt=0;ipt<npt;ipt++){
      fRun[ipt] =  fRun_ptr[ipt];
      fEntry[ipt] =  fEntry_ptr[ipt];
    }

    int npt_cut = eig_tree->Draw("Entry$","arm_flag==0 && kGood==2" +slug_cut[ic],"goff");
    double *fEntry_cut = eig_tree->GetV1();
    double cut_up = fEntry_cut[npt_cut-1];
    double cut_low = fEntry_cut[0];
    
    c1->cd(1);
    eig_tree->Draw("reg.us_avg_evMon2*1e3:Entry$",
		    "arm_flag==0 && kGood"+slug_cut[ic],"goff");
    TGraph *g_reg_slope = new TGraph(eig_tree->GetSelectedRows(),eig_tree->GetV2(),eig_tree->GetV1());
    g_reg_slope->SetMarkerStyle(20);
    g_reg_slope->SetMarkerColor(kRed);

    eig_tree->Draw("lagr.us_avg_evMon2*1e3:Entry$",
		   "arm_flag==0 && kGood" +slug_cut[ic],"goff");
    TGraph *g_lagr_slope = new TGraph(eig_tree->GetSelectedRows(),eig_tree->GetV2(),eig_tree->GetV1());
    g_lagr_slope->SetMarkerStyle(20);
    g_lagr_slope->SetMarkerColor(kBlue);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(g_reg_slope);
    mg->Add(g_lagr_slope);
    mg->Draw("AP");
    TH1D* hmg = (TH1D*)mg->GetHistogram();
    hmg->SetTitle(" us_avg_evMon2 Slopes(ppm/um) vs Run ; Run Number;Slope(ppm/um)");
    
    double x_min = fEntry[0];
    double x_max = fEntry[npt-1];
    double y_max = hmg->GetYaxis()->GetXmax();
    double y_min = hmg->GetYaxis()->GetXmin();
    hmg->GetXaxis()->Set(npt,x_min,x_max);
    hmg->GetXaxis()->SetTickLength(0.0);
    for(int ibin=1;ibin<=npt;ibin++)
      hmg->GetXaxis()->SetBinLabel(ibin, "");

    hmg->GetXaxis()->SetLabelSize(0.05);
    hmg->GetXaxis()->SetTitleSize(0.03);
    hmg->GetYaxis()->SetLabelSize(0.1);
    hmg->GetYaxis()->SetTitleSize(0.07);
    hmg->GetYaxis()->SetTitleOffset(0.8);

    Int_t last_run = -1;
    for(int ipt=0;ipt<npt;ipt++){
      if(fRun[ipt]!=last_run){
    	last_run = fRun[ipt];
    	hmg->GetXaxis()->SetBinLabel(ipt+1, Form("%d",(int)fRun[ipt]));
    	TLine *line = new TLine(x_min+ipt,y_min,x_min+ipt,y_max);
    	line->SetLineWidth(1);
    	line->SetLineStyle(7);
    	line->SetLineColorAlpha(14,0.5);
    	line->Draw("same");
      }
    }
    TLegend *leg = new TLegend(0.85,0.7,0.99,0.9);
    leg->AddEntry(g_reg_slope,"Regression","p");
    leg->AddEntry(g_lagr_slope,"Lagrange","p");
    leg->Draw("same");
    // TPaveText *pave_text = new TPaveText(fEntry[0],y_max,
    // 					 fEntry[npt/3],y_max+(y_max-y_min)*0.13);
    // pave_text->AddText(title[ic]);
    // pave_text->Draw();
    
    TBox *box1 = new TBox(cut_low, y_min,cut_up, y_max);
    box1->SetFillColorAlpha(kGray,0.5);
    box1->Draw();
    
    // Second one
    c1->cd(2);

    eig_tree->Draw("reg.us_avg_evMon0*1e3:Entry$",
		   "arm_flag==0 && kGood"+slug_cut[ic],"goff");
    TGraph *g_reg_slope2 = new TGraph(eig_tree->GetSelectedRows(),eig_tree->GetV2(),eig_tree->GetV1());
    g_reg_slope2->SetMarkerStyle(20);
    g_reg_slope2->SetMarkerColor(kRed);

    eig_tree->Draw("lagr.us_avg_evMon0*1e3:Entry$",
		   "arm_flag==0 && kGood" +slug_cut[ic],"goff");
    TGraph *g_lagr_slope2 = new TGraph(eig_tree->GetSelectedRows(),eig_tree->GetV2(),eig_tree->GetV1());
    g_lagr_slope2->SetMarkerStyle(20);
    g_lagr_slope2->SetMarkerColor(kBlue);

    TMultiGraph *mg2 = new TMultiGraph();
    mg2->Add(g_reg_slope2);
    mg2->Add(g_lagr_slope2);
    mg2->Draw("AP");
    TH1D* hmg2 = (TH1D*)mg2->GetHistogram();
    hmg2->SetTitle(" us_avg_evMon0 Slopes(ppm/um) vs Run ; Run Number;Slope(ppm/um)");
    
    x_min = fEntry[0];
    x_max = fEntry[npt-1];
    y_max = hmg2->GetYaxis()->GetXmax();
    y_min = hmg2->GetYaxis()->GetXmin();
    hmg2->GetXaxis()->Set(npt,x_min,x_max);
    hmg2->GetXaxis()->SetTickLength(0.0);
    for(int ibin=1;ibin<=npt;ibin++)
      hmg2->GetXaxis()->SetBinLabel(ibin, "");

    hmg2->GetXaxis()->SetLabelSize(0.05);
    hmg2->GetXaxis()->SetTitleSize(0.03);
    hmg2->GetYaxis()->SetLabelSize(0.1);
    hmg2->GetYaxis()->SetTitleSize(0.07);
    hmg2->GetYaxis()->SetTitleOffset(0.8);

    last_run = -1;
    for(int ipt=0;ipt<npt;ipt++){
      if(fRun[ipt]!=last_run){
    	last_run = fRun[ipt];
    	hmg2->GetXaxis()->SetBinLabel(ipt+1, Form("%d",(int)fRun[ipt]));
    	TLine *line = new TLine(x_min+ipt,y_min,x_min+ipt,y_max);
    	line->SetLineWidth(1);
    	line->SetLineStyle(7);
    	line->SetLineColorAlpha(14,0.5);
    	line->Draw("same");
      }
    }
    
    TBox *box2 = new TBox(cut_low, y_min,cut_up, y_max);
    box2->SetFillColorAlpha(kGray,0.5);
    box2->Draw();

    c1->cd(3);
    
    eig_tree->Draw("reg.us_avg_evMon1*1e3:Entry$",
		   "arm_flag==0 && kGood"+slug_cut[ic],"goff");
    TGraph *g_reg_slope3 = new TGraph(eig_tree->GetSelectedRows(),eig_tree->GetV2(),eig_tree->GetV1());
    g_reg_slope3->SetMarkerStyle(20);
    g_reg_slope3->SetMarkerColor(kRed);

    eig_tree->Draw("lagr.us_avg_evMon1*1e3:Entry$",
		   "arm_flag==0 && kGood" +slug_cut[ic],"goff");
    TGraph *g_lagr_slope3 = new TGraph(eig_tree->GetSelectedRows(),eig_tree->GetV2(),eig_tree->GetV1());
    g_lagr_slope3->SetMarkerStyle(20);
    g_lagr_slope3->SetMarkerColor(kBlue);

    TMultiGraph *mg3 = new TMultiGraph();
    mg3->Add(g_reg_slope3);
    mg3->Add(g_lagr_slope3);
    mg3->Draw("AP");
    TH1D* hmg3 = (TH1D*)mg3->GetHistogram();
    hmg3->SetTitle(" us_avg_evMon1 Slopes(ppm/um) vs Run ; Run Number;Slope(ppm/um)");
    
    x_min = fEntry[0];
    x_max = fEntry[npt-1];
    y_max = hmg3->GetYaxis()->GetXmax();
    y_min = hmg3->GetYaxis()->GetXmin();
    hmg3->GetXaxis()->Set(npt,x_min,x_max);
    hmg3->GetXaxis()->SetTickLength(0.0);
    hmg3->GetYaxis()->SetLabelSize(0.1);
    hmg3->GetYaxis()->SetTitleSize(0.07);
    hmg3->GetYaxis()->SetTitleOffset(0.8);

    for(int ibin=1;ibin<=npt;ibin++)
      hmg3->GetXaxis()->SetBinLabel(ibin, "");

    hmg3->GetXaxis()->SetLabelSize(0.05);
    hmg3->GetXaxis()->SetTitleSize(0.03);
    last_run = -1;
    for(int ipt=0;ipt<npt;ipt++){
      if(fRun[ipt]!=last_run){
    	last_run = fRun[ipt];
    	hmg3->GetXaxis()->SetBinLabel(ipt+1, Form("%d",(int)fRun[ipt]));
    	TLine *line = new TLine(x_min+ipt,y_min,x_min+ipt,y_max);
    	line->SetLineWidth(1);
    	line->SetLineStyle(7);
    	line->SetLineColorAlpha(14,0.5);
    	line->Draw("same");
      }
    }
    
    TBox *box3 = new TBox(cut_low, y_min,cut_up, y_max);
    box3->SetFillColorAlpha(kGray,0.5);
    box3->Draw();
    
    c1->SaveAs(title[ic]+"_slopes.pdf");
  }
  
}

void CompareSlopes6bpm(Int_t slug=94){
  gStyle->SetOptFit(1);
  TFile *input = TFile::Open(Form("rootfiles/MergedLagrange_slug%d.root",slug));
  TTree *mini_tree = (TTree*)input->Get("mini_eig6bpm");
  mini_tree->AddFriend("mini_lagr6bpm");
  mini_tree->AddFriend("mini");
  mini_tree->AddFriend("dit=slope",
		       Form("rootfiles/dit_eigslopes_6bpm_slug%d.root",slug));
  Int_t nBPM=6;
  
  vector<TString> DVlist={"asym_us_avg","asym_usl","asym_usr"};
  Int_t nDV = DVlist.size();

  for(int idv=0;idv<nDV;idv++){
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Print(Form("plots/slug%d_%s_compare_6bpm.pdf[",slug,DVlist[idv].Data()));
    // mini_tree->Draw(Form("%s.rms/ppm:Entry$",DVlist[idv].Data()),"","goff");
    // double *raw_val = mini_tree->GetV1();
    // double *rawx_val = mini_tree->GetV2();
    // TGraph *g1_raw_rms = new TGraph(mini_tree->GetSelectedRows(),
    // 				    rawx_val,raw_val);
    // g1_raw_rms->SetMarkerStyle(20);
    // g1_raw_rms->Draw("ALP");
    // g1_raw_rms->SetTitle("Raw RMS (ppm); Minirun# ; RMS(ppm)");
    // c1->Print(Form("plots/compare_6bpm_slug%d.pdf",slug));
    
    TMultiGraph *mgrms = new TMultiGraph(); // rms
    TLegend *legrms = new TLegend(0.9,0.7,1.0,0.9);
    mini_tree->Draw(Form("reg_%s.rms/ppm:Entry$",DVlist[idv].Data()),"","goff");
    double *y_val = mini_tree->GetV1();
    double *x_val = mini_tree->GetV2();
    TGraph *g1_reg_rms = new TGraph(mini_tree->GetSelectedRows(),
				    x_val,y_val);
    g1_reg_rms->SetMarkerStyle(20);
    g1_reg_rms->SetMarkerColor(kRed);
    g1_reg_rms->SetLineColor(kRed);

    mini_tree->Draw(Form("lagr_%s.rms/ppm:Entry$",DVlist[idv].Data()),"","goff");
    double *y_val2 = mini_tree->GetV1();
    double *x_val2 = mini_tree->GetV2();
    TGraph *g1_lag_rms = new TGraph(mini_tree->GetSelectedRows(),
				    x_val2,y_val2);
    g1_lag_rms->SetMarkerStyle(20);
    g1_lag_rms->SetMarkerColor(kBlue);
    g1_lag_rms->SetLineColor(kBlue);
    mgrms->Add(g1_reg_rms,"lp");
    mgrms->Add(g1_lag_rms,"lp");
    legrms->AddEntry(g1_reg_rms,"Eig-Reg");
    legrms->AddEntry(g1_lag_rms,"Lagrange");
    mgrms->Draw("A");
    mgrms->SetTitle("Corrected RMS (ppm): " + DVlist[idv] +";Minirun# ; RMS(ppm)");
    legrms->Draw("same");
    c1->Print(Form("plots/slug%d_%s_compare_6bpm.pdf",slug,DVlist[idv].Data()));

    mini_tree->Draw(Form("sqrt(lagr_%s.rms**2 - reg_%s.rms**2)/ppm:Entry$",DVlist[idv].Data(),DVlist[idv].Data()),"","goff");
    double *diff_y = mini_tree->GetV1();
    double *diff_x = mini_tree->GetV2();
    TGraph *g1_diff_rms = new TGraph(mini_tree->GetSelectedRows(),
				     diff_x,diff_y);
    g1_diff_rms->SetMarkerStyle(20);
    g1_diff_rms->Draw("ALP");
    g1_diff_rms->SetTitle(DVlist[idv] + ": #sqrt{#sigma^{2}_{lagr} -#sigma^{2}_{reg}} (ppm); Minirun# ; RMS(ppm)");
    c1->Print(Form("plots/slug%d_%s_compare_6bpm.pdf",slug,DVlist[idv].Data()));

    c1->cd();
    c1->Clear();
    TPad *pad1  = new TPad("pad1","pad1",0.0,0.0,0.7,1.0);
    TPad *pad2  = new TPad("pad2","pad2",0.7,0.0,1.0,1.0);
    pad1->Draw();
    pad2->Draw();
 
    TString det_name = DVlist[idv];
    det_name.ReplaceAll("asym_","");
    for(int i=0;i<nBPM;i++){
      gStyle->SetStatH(0.2);
      gStyle->SetStatW(0.3);
      c1->Clear("D");
      pad1->cd();
      mini_tree->Draw(Form("sign%d*mini_eig6bpm.%s_evMon%d/(ppm/um):Entry$",i,det_name.Data(),i),
		      "","goff");
      double *y_val = mini_tree->GetV1();
      double *x_val = mini_tree->GetV2();
      TGraph *g1_reg_slope = new TGraph(mini_tree->GetSelectedRows(),
					x_val,y_val);
      g1_reg_slope->SetMarkerStyle(20);
      g1_reg_slope->SetMarkerColor(kRed);
      g1_reg_slope->SetLineColor(kRed);
    
      mini_tree->Draw(Form("sign%d*dit.%s_evMon%d/(ppm/um):Entry$",i,det_name.Data(),i),
		      "","goff");
      double *y_val2 = mini_tree->GetV1();
      double *x_val2 = mini_tree->GetV2();
      TGraph *g1_lag_slope = new TGraph(mini_tree->GetSelectedRows(),
					x_val2,y_val2);
      g1_lag_slope->SetMarkerStyle(20);
      g1_lag_slope->SetMarkerColor(kBlue);
      g1_lag_slope->SetLineColor(kBlue);
    
      TMultiGraph *mgslope = new TMultiGraph(); // rms
      TLegend *legslope = new TLegend(0.9,0.7,1.0,0.9);
    
      mgslope->Add(g1_reg_slope,"lp");
      mgslope->Add(g1_lag_slope,"lp");
      legslope->AddEntry(g1_reg_slope,"Eig-Reg");
      legslope->AddEntry(g1_lag_slope,"Lagrange");
      mgslope->Draw("A");
      mgslope->SetTitle(Form("Slope : %s_vs_evMon%d (ppm/um) ; minirun # ; Slope (ppm/um)",
			     det_name.Data(),i));
      legslope->Draw("same");
      pad2->cd();
      mini_tree->Draw(Form("dit.%s_evMon%d/mini_eig6bpm.%s_evMon%d",
			   det_name.Data(),i,det_name.Data(),i),"","");
      TH1D *h1d = (TH1D*)gPad->FindObject("htemp");
      h1d->SetTitle("ratio (dit_slope /reg_slope) in eigen-basis");

      c1->Print(Form("plots/slug%d_%s_compare_6bpm.pdf",slug,DVlist[idv].Data())); 
    }
    c1->Print(Form("plots/slug%d_%s_compare_6bpm.pdf]",slug,DVlist[idv].Data())); 
  }

 

  // TMultiGraph *mgavg = new TMultiGraph(); // avg
  // TLegend *legavg = new TLegend(0.9,0.7,1.0,0.9);
  // mini_tree->Draw("reg_asym_us_avg/ppb:reg_asym_us_avg.err/ppb:Entry$","","goff");
  // double *y_mean = mini_tree->GetV1();
  // double *y_err = mini_tree->GetV2();
  // double *x_mean = mini_tree->GetV3();
  // TGraphErrors *g1_reg_avg = new TGraphErrors(mini_tree->GetSelectedRows(),
  // 					      x_mean,y_mean,0,y_err);
  
  // g1_reg_avg->SetMarkerStyle(20);
  // g1_reg_avg->SetMarkerColor(kRed);
  // g1_reg_avg->SetLineColor(kRed);
  // mini_tree->Draw("lagr_asym_us_avg/ppb:lagr_asym_us_avg.err/ppb:Entry$+0.4","","goff");
  // double *y_mean2 = mini_tree->GetV1();
  // double *y_err2 = mini_tree->GetV2();
  // double *x_mean2 = mini_tree->GetV3();
  // TGraphErrors *g1_lag_avg = new TGraphErrors(mini_tree->GetSelectedRows(),
  // 					      x_mean2,y_mean2,0,y_err2);

  // g1_lag_avg->SetMarkerStyle(20);
  // g1_lag_avg->SetMarkerColor(kBlue);
  // g1_lag_avg->SetLineColor(kBlue);
  // mgavg->Add(g1_reg_avg,"p");
  // mgavg->Add(g1_lag_avg,"p");
  // legavg->AddEntry(g1_reg_avg,"Eig-Reg");
  // legavg->AddEntry(g1_lag_avg,"Lagrange");
  // mgavg->Draw("A");
  // double ymax = mgavg->GetYaxis()->GetXmax();
  // double ymin = mgavg->GetYaxis()->GetXmin();
  // mgavg->GetYaxis()->SetRangeUser(ymin,ymax+0.5*(ymax-ymin));
  // mgavg->SetTitle("Corrected Mean (ppb); Minirun# ; Mean(ppb)");
  // legavg->Draw("same");

  // g1_reg_avg->Fit("pol0");
  // gPad->Update();
  // TPaveStats *ps1  =(TPaveStats*)g1_reg_avg->FindObject("stats");
  // ps1->SetName("reg");
  // ps1->SetTextColor(kRed);
  // ps1->SetX1NDC(0.1);
  // ps1->SetX2NDC(0.3);
  // ps1->SetY1NDC(0.8);
  // ps1->SetY2NDC(0.9);
  // g1_lag_avg->Fit("pol0");
  // g1_lag_avg->GetFunction("pol0")->SetLineColor(kBlue);
  // g1_lag_avg->GetFunction("pol0")->SetLineStyle(9);
  // gPad->Update();
  // TPaveStats *ps2  =(TPaveStats*)g1_lag_avg->FindObject("stats");
  // ps2->SetName("lag");
  // ps2->SetTextColor(kBlue);
  // ps2->SetX1NDC(0.1);
  // ps2->SetX2NDC(0.3);
  // ps2->SetY1NDC(0.8);
  // ps2->SetY2NDC(0.7);

  // c1->Print(Form("plots/compare_6bpm_slug%d.pdf",slug));

  
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Print(Form("plots/slug%d_diff_evMon_6bpm.pdf[",slug));
  c2->cd();
  for(int i=0;i<nBPM;i++){
    TString channel = Form("diff_evMon%d",i);
    mini_tree->Draw(Form("sign%d*%s/nm:%s.err/nm:Entry$",i,channel.Data(),channel.Data())
		    ,"","goff");
    double *y_mean = mini_tree->GetV1();
    double *y_err = mini_tree->GetV2();
    double *x_mean = mini_tree->GetV3();
    TGraphErrors *g1_avg = new TGraphErrors(mini_tree->GetSelectedRows(),
						x_mean,y_mean,0,y_err);
  
    g1_avg->SetMarkerStyle(20);
    g1_avg->Draw("AP");
    g1_avg->SetTitle(channel+" Mean (nm) ; minirun# ; diff (nm)");
    g1_avg->Fit("pol0","Q");
    c2->Print(Form("plots/slug%d_diff_evMon_6bpm.pdf",slug));

    mini_tree->Draw(Form("%s.rms/um:Entry$",channel.Data())
		    ,"","goff");
    double *y_rms = mini_tree->GetV1();
    double *x_mean2 = mini_tree->GetV2();
    TGraph *g1_rms = new TGraph(mini_tree->GetSelectedRows(),
				x_mean2,y_mean);
  
    g1_rms->SetMarkerStyle(20);
    g1_rms->Draw("AP");
    g1_rms->SetTitle(channel+" RMS (um) ; minirun# ; RMS (um)");
    g1_rms->Fit("pol0","Q");
    c2->Print(Form("plots/slug%d_diff_evMon_6bpm.pdf",slug));

  }
  c2->Print(Form("plots/slug%d_diff_evMon_6bpm.pdf]",slug));

  // for(int i=0;i<nBPM;i++){
  //   gStyle->SetStatH(0.2);
  //   gStyle->SetStatW(0.3);
  //   c1->Clear("D");
  //   pad1->cd();
  //   mini_tree->Draw(Form("mini_eig.us_avg_evMon%d*diff_evMon%d/ppb:Entry$",i,i),"","goff");
  //   double *y_val = mini_tree->GetV1();
  //   double *x_val = mini_tree->GetV2();
  //   TGraph *g1_reg_slope = new TGraph(mini_tree->GetSelectedRows(),
  // 				      x_val,y_val);
  //   g1_reg_slope->SetMarkerStyle(20);
  //   g1_reg_slope->SetMarkerColor(kRed);
  //   g1_reg_slope->SetLineColor(kRed);
    
  //   mini_tree->Draw(Form("dit.us_avg_evMon%d*diff_evMon%d/ppb:Entry$",i,i),"","goff");
  //   double *y_val2 = mini_tree->GetV1();
  //   double *x_val2 = mini_tree->GetV2();
  //   TGraph *g1_lag_slope = new TGraph(mini_tree->GetSelectedRows(),
  // 				      x_val2,y_val2);
  //   g1_lag_slope->SetMarkerStyle(20);
  //   g1_lag_slope->SetMarkerColor(kBlue);
  //   g1_lag_slope->SetLineColor(kBlue);
    
  //   TMultiGraph *mgslope = new TMultiGraph(); // rms
  //   TLegend *legslope = new TLegend(0.9,0.7,1.0,0.9);
    
  //   mgslope->Add(g1_reg_slope,"lp");
  //   mgslope->Add(g1_lag_slope,"lp");
  //   legslope->AddEntry(g1_reg_slope,"Eig-Reg");
  //   legslope->AddEntry(g1_lag_slope,"Lagrange");
  //   mgslope->Draw("A");
  //   mgslope->SetTitle(Form("Correction on us_avg by evMon%d (ppb) ; minirun # ; Correction (ppb)",i));
  //   legslope->Draw("same");
  //   pad2->cd();
  //   mini_tree->Draw(Form("(dit.us_avg_evMon%d-mini_eig.us_avg_evMon%d)*diff_evMon%d/ppb",i,i,i),"","");
  //   TH1D *h1d = (TH1D*)gPad->FindObject("htemp");
  //   h1d->SetTitle("Correction Difference (lagr - eigReg) (ppb) ");
  //   c1->Print(Form("plots/compare_6bpm_slug%d.pdf",slug));
  // }

}

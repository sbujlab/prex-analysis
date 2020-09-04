void CompareSlopes_allrun(TString det_name="us_avg",TString cut_text="arm_flag==0"){
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();
  TPad *pad1 = new TPad("pad1","",0,0.5, 0.66, 1.0);
  TPad *pad2 = new TPad("pad2","",0,0, 0.66,0.5);
  TPad *pad3 = new TPad("pad3","",0.66,0, 1.0, 0.5);
  TPad *pad4 = new TPad("pad4","",0.66,0.5, 1.0,1.0);
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  c1->SaveAs(Form("plots/slopes_%s_allrun.pdf[",det_name.Data()));
  Int_t nIV =12;
  TFile* input = TFile::Open("rootfiles/all_slugs.root");
  TTree *eig_tree = (TTree*)input->Get("eig");
  TTree *reg_tree = (TTree*)input->Get("reg");
  TTree *lagr_tree = (TTree*)input->Get("lagr");
  eig_tree->AddFriend(reg_tree);
  eig_tree->AddFriend(lagr_tree);
  vector<double> fDiff;
  vector<double> fDiffFrac;
  vector<double> fDiffWidth;
  vector<double> fSlope;
  TH1D *h1dptr;
  for(int iiv=0; iiv<nIV;iiv++){
    pad1->cd();
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    eig_tree->Draw(Form(" (reg.%s_evMon%d)*1e3 : run ", det_name.Data(), iiv),
		   cut_text);
    double *reg_slope = eig_tree->GetV1();
    double *reg_run = eig_tree->GetV2();
    TGraph *g_reg_slope  = new TGraph(eig_tree->GetSelectedRows(), reg_run,reg_slope);
    g_reg_slope->SetMarkerStyle(22);
    g_reg_slope->SetMarkerSize(0.5);
    g_reg_slope->SetMarkerColor(kRed);
    eig_tree->Draw(Form(" (lagr.%s_evMon%d)*1e3 : run ", det_name.Data(), iiv),
		   cut_text);
    double *lagr_slope = eig_tree->GetV1();
    double *lagr_run = eig_tree->GetV2();
    TGraph *g_lagr_slope  = new TGraph(eig_tree->GetSelectedRows(), lagr_run,lagr_slope);
    g_lagr_slope->SetMarkerStyle(23);
    g_lagr_slope->SetMarkerSize(0.5);
    g_lagr_slope->SetMarkerColor(kBlue);

    mg->SetTitle( Form(" %s_evMon%d Slope(ppm/um)", det_name.Data(),iiv));
    mg->Add(g_reg_slope,"p");
    mg->Add(g_lagr_slope,"p");
    leg->AddEntry(g_reg_slope,"regression","p");
    leg->AddEntry(g_lagr_slope,"Lagrange","p");
    mg->Draw("A");
    double yaxis_max = mg->GetYaxis()->GetXmax();
    double yaxis_min = mg->GetYaxis()->GetXmin();
    mg->GetYaxis()->SetRangeUser(yaxis_min, yaxis_max + 0.3*( yaxis_max-yaxis_min));
    gPad->Modified();
    gPad->Update();
    mg->GetYaxis()->SetTitle("Slope (ppm/um)");
    mg->GetXaxis()->SetTitle("run number");
    leg->Draw("same");
    
    pad4->cd();
    eig_tree->Draw(Form(" (reg.%s_evMon%d)*1e3 :(lagr.%s_evMon%d)*1e3  ", det_name.Data(), iiv,det_name.Data(), iiv),
		   cut_text);
    
    pad2->cd();
    eig_tree->Draw(Form(" (lagr.%s_evMon%d - reg.%s_evMon%d)*1e3 : run ", det_name.Data(), iiv,det_name.Data(), iiv),
		   cut_text,"goff");
    double *diff_slope = eig_tree->GetV1();
    double *diff_run  = eig_tree->GetV2();
    TGraph *g_diff_slope = new TGraph(eig_tree->GetSelectedRows(),diff_run,diff_slope);
    g_diff_slope->SetMarkerStyle(7);
    g_diff_slope->SetMarkerSize(0.5);
    g_diff_slope->Draw("AP");
    g_diff_slope->SetTitle( Form(" Difference in slopes : Lagrange - Regression ; run number;Slope(ppm/um) "));
    

    pad3->cd();
    eig_tree->Draw(Form("lagr.%s_evMon%d*1e3 ", det_name.Data(), iiv),
		   cut_text);
    h1dptr = (TH1D*)gPad->FindObject("htemp");
    fSlope.push_back( h1dptr->GetMean() );

    eig_tree->Draw(Form(" (lagr.%s_evMon%d - reg.%s_evMon%d)*1e3 ", det_name.Data(), iiv,det_name.Data(), iiv),
		   cut_text);
    h1dptr = (TH1D*)gPad->FindObject("htemp");
    fDiff.push_back( h1dptr->GetMean() );
    fDiffWidth.push_back( h1dptr->GetRMS() );
    c1->SaveAs(Form("plots/slopes_%s_allrun.pdf",det_name.Data()));
  }
  int nlines = fDiff.size();
  for(int i=0;i<nlines;i++){
    printf("evMon%d & %.1f \t  & %.1f \t & %.1f \\\\ \n",
	   i, fSlope[i], fDiff[i],fDiffWidth[i]);
  }

  c1->SaveAs(Form("plots/slopes_%s_allrun.pdf]",det_name.Data()));
}

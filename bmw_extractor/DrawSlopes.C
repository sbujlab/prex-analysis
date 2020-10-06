void DrawSlopes(TString det_name,TString cut_text, Int_t slug_start, Int_t slug_end);
void DrawSlopes(){
  DrawSlopes("us_avg","arm_flag==0",3,94);   
  DrawSlopes("us_dd","arm_flag==0",3,94);
  DrawSlopes("us_avg","arm_flag==0",1,2);   
  DrawSlopes("us_dd","arm_flag==0",1,2);
}

void DrawSlopes(TString det_name,TString cut_text, Int_t slug_start, Int_t slug_end){
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->cd();
  TPad *pad1 = new TPad("pad1","",0,0.5, 0.66, 1.0);
  TPad *pad2 = new TPad("pad2","",0,0, 0.66,0.5);
  TPad *pad1b = new TPad("pad1b","", 0.66, 0.5, 1.,1);
  TPad *pad2b = new TPad("pad2b","", 0.66, 0.0, 1,0.5);

  pad1->Draw();  pad1b->Draw();
  pad2->Draw();  pad2b->Draw();
  pad1->SetRightMargin(0.02);
  pad2->SetRightMargin(0.02);
  TString plot_filename = Form("plots/slopes_evectors_%s_slug%d-%d.pdf",det_name.Data(),slug_start,slug_end);
  c1->SaveAs(plot_filename+"[");
  Int_t nIV ;
  if( slug_start>=3)
    nIV=12;
  if( slug_end<=2)
    nIV=10;

  TChain *eig_tree = new TChain("eig");
  TChain *reg_tree = new TChain("reg");
  TChain *lagr_tree = new TChain("lagr");
  for(int islug=slug_start;islug<=slug_end;islug++){
    eig_tree->AddFile(Form("./rootfiles/slug%d_sorted_eigenvector_allbpm.root",islug));
    reg_tree->AddFile(Form("./rootfiles/slug%d_sorted_eigenvector_allbpm.root",islug));
    lagr_tree->AddFile(Form("./rootfiles/slug%d_sorted_eigenvector_allbpm.root",islug));
  }
  eig_tree->AddFriend(reg_tree);
  eig_tree->AddFriend(lagr_tree);
  int npt  = eig_tree->Draw("slug", cut_text,"goff");
  double *slug_ptr = eig_tree->GetV1();
  double *fCounter  = new double[npt];
  vector<Double_t> fSlug;
  for(int ipt=0;ipt<npt;ipt++){
    fSlug.push_back( slug_ptr[ipt] );
    fCounter[ipt] = ipt;
  }
  TH1D *h1dptr;
  for(int iiv=0; iiv<nIV;iiv++){
    pad1b->cd();
    eig_tree->Draw(Form(" (reg.%s_evMon%d)*1e3  ", det_name.Data(), iiv),cut_text);
    double *reg_slope = eig_tree->GetV1();
    TGraph *g_reg_slope  = new TGraph(eig_tree->GetSelectedRows(), fCounter,reg_slope);
    g_reg_slope->SetMarkerStyle(6);
    g_reg_slope->SetMarkerColor(kRed);
    g_reg_slope->SetLineColor(kRed);
    
    eig_tree->Draw(Form(" (lagr.%s_evMon%d)*1e3 ", det_name.Data(), iiv),cut_text);
    h1dptr = (TH1D*)gPad->FindObject("htemp");
    h1dptr->SetTitle(" Lagrange Slope(ppm/um)");

    pad1->cd();
    double *lagr_slope = eig_tree->GetV1();
    TGraph *g_lagr_slope  = new TGraph(eig_tree->GetSelectedRows(), fCounter,lagr_slope);
    g_lagr_slope->SetMarkerStyle(6);
    g_lagr_slope->SetMarkerColor(kBlue);
    g_lagr_slope->SetLineColor(kBlue);
    
    TLegend *leg = new TLegend(0.99,0.99,0.8,0.8);
    leg->AddEntry(g_lagr_slope,"Lagrange","lp");
    leg->AddEntry(g_reg_slope,"Regression","lp");
    // TH1F *hlagr_slope  = g_lagr_slope->GetHistogram();
    // hlagr_slope->GetXaxis()->Set(npt,-0.5,npt-0.5);
    // for(int ipt=1;ipt<=npt;ipt++){
    //   hlagr_slope->GetXaxis()->SetBinLabel(ipt, Form( "%d.%d",(int)fRun[ipt],(int)fMini[ipt]));
    // }
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(g_reg_slope);
    mg->Add(g_lagr_slope);
    mg->Draw("AP");
    mg->SetTitle(Form(" Slope: %s_evMon%d (ppm/um); Slug ; Slope(ppm/um)",
		      det_name.Data(),iiv));

    double y_max =  mg->GetYaxis()->GetXmax();
    double y_min =  mg->GetYaxis()->GetXmin();
    mg->GetYaxis()->SetRangeUser(y_min, y_max + 0.3*(y_max-y_min));
    mg->GetXaxis()->SetTickLength(0.0);
    leg->Draw("same");
    TH1F *hmg_slope  = mg->GetHistogram();
    Int_t nBin = npt;
    hmg_slope->GetXaxis()->Set(nBin,-0.5,npt-0.5);
    Int_t current_slug=-1;
    for(int ibin=0;ibin<nBin;ibin++){
      if(fSlug[ibin]!=current_slug){
	current_slug = fSlug[ibin];
	hmg_slope->GetXaxis()->SetBinLabel(ibin+1, Form( "%d",(int)fSlug[ibin]));
	TLine *line = new TLine(ibin+1,y_min, ibin+1,y_max);
	line->SetLineWidth(1);
	line->SetLineStyle(3);
	line->SetLineColorAlpha(14,0.5);
	line->Draw("same");
      }
    }

    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.6);
    hmg_slope->GetXaxis()->SetLabelSize(0.04);

    pad2b->cd();
    eig_tree->Draw(Form(" (reg.%s_evMon%d - lagr.%s_evMon%d)*1e3 ", det_name.Data(), iiv,det_name.Data(), iiv), cut_text);
    h1dptr = (TH1D*)gPad->FindObject("htemp");

    double *diff_slope = eig_tree->GetV1();
    TGraph *g_diff_slope = new TGraph(eig_tree->GetSelectedRows(),fCounter,diff_slope);
    g_diff_slope->SetMarkerStyle(7);      

    pad2->cd();
    g_diff_slope->Draw("APL");
    g_diff_slope->SetTitle(" Difference in slopes : regression - Lagrange ; Slug;Slope(ppm/um) ");
    g_diff_slope->GetXaxis()->SetTickLength(0.0);
    g_diff_slope->SetMarkerStyle(6);
    TH1F *hdiff_slope  = g_diff_slope->GetHistogram();
    hdiff_slope->GetXaxis()->Set(npt,-0.5,npt-0.5);
    y_max = g_diff_slope->GetYaxis()->GetXmax();
    y_min = g_diff_slope->GetYaxis()->GetXmin();
    current_slug=-1;
    for(int ibin=0;ibin<npt;ibin++){
      if(fSlug[ibin]!=current_slug){
	current_slug = fSlug[ibin];
	hdiff_slope->GetXaxis()->SetBinLabel(ibin+1, Form( "%d",(int)fSlug[ibin]));
	TLine *line = new TLine(ibin+1,y_min, ibin+1,y_max);
	line->SetLineWidth(1);
	line->SetLineStyle(3);
	line->SetLineColorAlpha(14,0.5);
	line->Draw("same");
      }
    }

    g_diff_slope->GetYaxis()->SetTitleSize(0.07);
    g_diff_slope->GetYaxis()->SetTitleOffset(0.6);
    hdiff_slope->GetXaxis()->SetLabelSize(0.04);
    c1->SaveAs(plot_filename);
  }

  c1->SaveAs(plot_filename+"]");
}

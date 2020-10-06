void DrawRMS(){
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  TCanvas *c1 = new TCanvas("c1","c1",1600,600);
  c1->cd();
  c1->SetRightMargin(0.02);
  TString plot_filename = "rms.pdf";
  
  TFile *input = TFile::Open("bmw_plus_allslugs.root");
  
  TTree *mini_tree = (TTree*)input->Get("mini");
  TTree *info_tree = (TTree*)input->Get("mini_info");
  mini_tree->AddFriend(info_tree);

  TString cut_text = "arm_flag==0 && kGood==1";
  int npt  = mini_tree->Draw("slug", cut_text,"goff");
  double *slug_ptr = mini_tree->GetV1();
  double *run_ptr = mini_tree->GetV2();
  double *fCounter  = new double[npt];
  vector<Double_t> fSlug;
  for(int ipt=0;ipt<npt;ipt++){
    fSlug.push_back( slug_ptr[ipt] );
    fCounter[ipt] = ipt;
  }
  TH1D *h1dptr;

  mini_tree->Draw("reg_asym_us_avg.rms*1e6:Entry$",cut_text);
  double *reg_rms = mini_tree->GetV1();
  TGraph *g_reg_rms  = new TGraph(mini_tree->GetSelectedRows(), fCounter,reg_rms);
  g_reg_rms->SetMarkerStyle(6);
    
  g_reg_rms->Draw("AP");
  g_reg_rms->SetTitle(" reg_asym_us_avg RMS (ppm) vs slug; Slug ; RMS(ppm)");


  double y_max =  g_reg_rms->GetYaxis()->GetXmax();
  double y_min =  g_reg_rms->GetYaxis()->GetXmin();
  g_reg_rms->GetYaxis()->SetRangeUser(y_min, y_max + 0.3*(y_max-y_min));
  g_reg_rms->GetXaxis()->SetTickLength(0.0);
  TH1F *hreg_rms  = g_reg_rms->GetHistogram();
  Int_t nBin = npt;
  hreg_rms->GetXaxis()->Set(nBin,-0.5,npt-0.5);
  Int_t current_slug=-1;
  for(int ibin=0;ibin<nBin;ibin++){
    if(fSlug[ibin]!=current_slug){
      current_slug = fSlug[ibin];
      hreg_rms->GetXaxis()->SetBinLabel(ibin+1, Form( "%d",(int)fSlug[ibin]));
      TLine *line = new TLine(ibin+1,y_min, ibin+1,y_max);
      line->SetLineWidth(1);
      line->SetLineStyle(3);
      line->SetLineColorAlpha(14,0.5);
      line->Draw("same");
    }
  }

  g_reg_rms->GetYaxis()->SetTitleSize(0.07);
  g_reg_rms->GetYaxis()->SetTitleOffset(0.6);
  hreg_rms->GetXaxis()->SetLabelSize(0.04);
  c1->SaveAs(plot_filename);
}

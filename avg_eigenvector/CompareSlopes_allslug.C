void CompareSlopes_allslug(TString det_name="us_avg",TString cut_text="arm_flag==0"){
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->cd();
  TPad *pad1 = new TPad("pad1","",0,0, 0.66, 0.33);
  TPad *pad2 = new TPad("pad2","",0,0.33, 0.66,0.66);
  TPad *pad3 = new TPad("pad3","",0,0.66, 0.66, 1.0);  

  TPad *pad1b = new TPad("pad1b","", 0.66, 0.0, 1.0,0.33);
  TPad *pad2b = new TPad("pad2b","", 0.66, 0.33, 1.0,0.66);
  TPad *pad3b = new TPad("pad3b","", 0.66, 0.66, 1.0,1.0);

  pad1->Draw();  pad1b->Draw();
  pad2->Draw();  pad2b->Draw();
  pad3->Draw();  pad3b->Draw();

  c1->SaveAs(Form("plots/slopes_%s_slug.pdf[",det_name.Data()));
  Int_t nIV =12;
  TFile* input = TFile::Open("rootfiles/all_slugs.root");
  TTree *eig_tree = (TTree*)input->Get("eig");
  TTree *reg_tree = (TTree*)input->Get("reg");
  TTree *lagr_tree = (TTree*)input->Get("lagr");
  eig_tree->AddFriend(reg_tree);
  eig_tree->AddFriend(lagr_tree);
  vector<Double_t> fRegSlopeBySlug;
  vector<Double_t> fLagrSlopeBySlug;
  vector<Double_t> fDiffSlopeMeanBySlug;
  vector<Double_t> fDiffSlopeRMSBySlug;
  vector<Double_t> fSlugIdx;
  for(int islug=1;islug<=94;islug++){
    TString slug_cut = Form("&& slug==%d",islug);
    int npt  = eig_tree->Draw("run:mini", cut_text + slug_cut,"goff");
    if(npt == 0)
      continue;
    fSlugIdx.push_back(islug);
    double *run_ptr = eig_tree->GetV1();
    double *mini_ptr = eig_tree->GetV2();
    vector<double> fRun;
    vector<double> fMini;
    double *fCounter  = new double[npt];
    for(int ipt=0;ipt<npt;ipt++){
      fRun.push_back( run_ptr[ipt] );
      fMini.push_back( mini_ptr[ipt] );
      fCounter[ipt] = ipt;
    }
    TH1D *h1dptr;
    for(int iiv=0; iiv<nIV;iiv++){
      pad1b->cd();
      eig_tree->Draw(Form(" (reg.%s_evMon%d)*1e3  ", det_name.Data(), iiv),
      		     cut_text+slug_cut);
      h1dptr = (TH1D*)gPad->FindObject("htemp");
      fRegSlopeBySlug.push_back(h1dptr->GetMean());
      double *reg_slope = eig_tree->GetV1();
      TGraph *g_reg_slope  = new TGraph(eig_tree->GetSelectedRows(), fCounter,reg_slope);
      g_reg_slope->SetMarkerStyle(22);
      g_reg_slope->SetMarkerColor(kRed);
      g_reg_slope->SetLineColor(kRed);
      pad1->cd();
      g_reg_slope->Draw("APL");
      g_reg_slope->SetTitle( Form("Slug %d, %s_evMon%d: Regression Slope(ppm/um)", islug,det_name.Data(),iiv));
      TH1F *hreg_slope  = g_reg_slope->GetHistogram();
      hreg_slope->GetXaxis()->Set(npt,-0.5,npt-0.5);
      for(int ipt=1;ipt<=npt;ipt++){
	hreg_slope->GetXaxis()->SetBinLabel(ipt, Form( "%d.%d",(int)fRun[ipt],(int)fMini[ipt]));
      }
      
      pad2b->cd();
      eig_tree->Draw(Form(" (lagr.%s_evMon%d)*1e3 ", det_name.Data(), iiv),
      		     cut_text+slug_cut);
      h1dptr = (TH1D*)gPad->FindObject("htemp");
      fLagrSlopeBySlug.push_back(h1dptr->GetMean());

      double *lagr_slope = eig_tree->GetV1();
      TGraph *g_lagr_slope  = new TGraph(eig_tree->GetSelectedRows(), fCounter,lagr_slope);
      g_lagr_slope->SetMarkerStyle(23);
      g_lagr_slope->SetMarkerColor(kBlue);
      g_lagr_slope->SetLineColor(kBlue);
      pad2->cd();
      g_lagr_slope->Draw("APL");
      g_lagr_slope->SetTitle( Form("Slug %d, %s_evMon%d: Lagrange Slope(ppm/um)", islug,det_name.Data(),iiv));
      TH1F *hlagr_slope  = g_lagr_slope->GetHistogram();
      hlagr_slope->GetXaxis()->Set(npt,-0.5,npt-0.5);
      for(int ipt=1;ipt<=npt;ipt++){
	hlagr_slope->GetXaxis()->SetBinLabel(ipt, Form( "%d.%d",(int)fRun[ipt],(int)fMini[ipt]));
      }

      pad3b->cd();
      eig_tree->Draw(Form(" (reg.%s_evMon%d - lagr.%s_evMon%d)*1e3 ", det_name.Data(), iiv,det_name.Data(), iiv),
		     cut_text+slug_cut);
      h1dptr = (TH1D*)gPad->FindObject("htemp");
      fDiffSlopeMeanBySlug.push_back(h1dptr->GetMean());
      fDiffSlopeRMSBySlug.push_back(h1dptr->GetRMS());
      
      double *diff_slope = eig_tree->GetV1();
      TGraph *g_diff_slope = new TGraph(eig_tree->GetSelectedRows(),fCounter,diff_slope);
      g_diff_slope->SetMarkerStyle(7);      

      pad3->cd();
      g_diff_slope->Draw("APL");
      g_diff_slope->SetTitle( Form(" Slug %d, Difference in slopes : regression - Lagrange ; run number;Slope(ppm/um) ",islug));
      TH1F *hdiff_slope  = g_diff_slope->GetHistogram();
      hdiff_slope->GetXaxis()->Set(npt,-0.5,npt-0.5);
      for(int ipt=1;ipt<=npt;ipt++){
	hdiff_slope->GetXaxis()->SetBinLabel(ipt, Form( "%d.%d",(int)fRun[ipt],(int)fMini[ipt]));
      }
      

      c1->SaveAs(Form("plots/slopes_%s_slug.pdf",det_name.Data()));
    }
  }
  
  c1->SaveAs(Form("plots/slopes_%s_slug.pdf]",det_name.Data()));
  TCanvas* c2 = new TCanvas("c2","c2",1200,1000);
  c2->Divide(1,3);
  Int_t nSlug = fSlugIdx.size();
  double *fslug = new double[nSlug];
  for(int i=0;i<nSlug;i++)
    fslug[i] = fSlugIdx[i];
  c2->SaveAs(Form("compare_slope_%s_summary.pdf[",det_name.Data()));
  for(int iiv=0;iiv<nIV;iiv++){
    double *reg_slope  = new double[nSlug];
    double *lagr_slope  = new double[nSlug];
    double *diff_slope_mean  = new double[nSlug];
    double *diff_slope_rms  = new double[nSlug];
    for(int islug=0;islug<nSlug;islug++){
      diff_slope_mean[islug]  = fDiffSlopeMeanBySlug[islug*nIV+iiv];
      diff_slope_rms[islug]  = fDiffSlopeRMSBySlug[islug*nIV+iiv];
      reg_slope[islug] = fRegSlopeBySlug[islug*nIV+iiv];
      lagr_slope[islug] = fLagrSlopeBySlug[islug*nIV+iiv];
    }
    c2->cd(1);
    TMultiGraph *mg  = new TMultiGraph();
    TGraph *g_reg_slope = new TGraph(nSlug,fslug,reg_slope);
    g_reg_slope->SetMarkerStyle(20);
    g_reg_slope->SetMarkerColor(kRed);
    g_reg_slope->SetLineColor(kRed);
    TGraph *g_lagr_slope = new TGraph(nSlug,fslug,lagr_slope);
    g_lagr_slope->SetMarkerStyle(20);
    g_lagr_slope->SetMarkerColor(kBlue);
    g_lagr_slope->SetLineColor(kBlue);
    mg->Add(g_reg_slope,"lp");
    mg->Add(g_lagr_slope,"lp");
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(g_reg_slope,"Regression","p");
    leg->AddEntry(g_lagr_slope,"Lagrange","p");
    mg->Draw("A");
    mg->SetTitle( Form("%s evMon%d Slopes; slug; Slope(ppm/um)",det_name.Data(), iiv));
    double ymax = mg->GetYaxis()->GetXmax();
    double ymin = mg->GetYaxis()->GetXmin();
    mg->GetYaxis()->SetRangeUser(ymin, ymax+0.3*(ymax-ymin));
    mg->GetXaxis()->SetRangeUser(0,100);
    gPad->Modified();
    gPad->Update();
    leg->Draw("same");
    
    c2->cd(2);
    TGraph *g_diff_slope_mean = new TGraph(nSlug,fslug,diff_slope_mean);
    g_diff_slope_mean->SetMarkerStyle(20);
    g_diff_slope_mean->Draw("APL");
    g_diff_slope_mean->SetTitle("Slope Differences : Mean vs Slugs; Slug; Mean(ppm/um)");
    g_diff_slope_mean->GetXaxis()->SetRangeUser(0,100);
    c2->cd(3);
    TGraph *g_diff_slope_rms = new TGraph(nSlug,fslug,diff_slope_rms);
    g_diff_slope_rms->SetMarkerStyle(20);
    g_diff_slope_rms->Draw("APL");
    g_diff_slope_rms->SetTitle("Slope Differences : RMS vs Slugs; Slug; RMS(ppm/um)");
    g_diff_slope_rms->GetXaxis()->SetRangeUser(0,100);
    c2->SaveAs(Form("compare_slope_%s_summary.pdf",det_name.Data()));
  }
  c2->SaveAs(Form("compare_slope_%s_summary.pdf]",det_name.Data()));

  int nslug = fSlugIdx.size();
  FILE *output_log = fopen("output.log","w");
  for(int i=0;i<nslug;i++){
    fprintf(output_log,"%d \t", i);
    for(int iev=0;iev<nIV;iev++){
      fprintf(output_log, "%f \t", fDiffSlopeMeanBySlug[i*nIV+iev]);
      fprintf(output_log, "%f \t", fDiffSlopeRMSBySlug[i*nIV+iev]);
      double quad_sum = sqrt( pow(fDiffSlopeMeanBySlug[i*nIV+iev],2)+pow(fDiffSlopeRMSBySlug[i*nIV+iev],2));
      fprintf(output_log, "%f \t", quad_sum);
    }
    fprintf(output_log,"\n");
  }
  fclose(output_log);
}

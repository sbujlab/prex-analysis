void CompareSlopes(Int_t kSwitch, Int_t run){
  gStyle->SetOptFit(1);
  TFile *input = TFile::Open(Form("rootfiles/mc_lagrange_%d.000.root",run));
  TTree *mini_tree;
  Int_t nBPM=0;
  TString run_title=Form("Run%d: ",run);
  TString filename_tag="";
  TString eig_treename="";
  TString lagr_treename="";
  TString dit_prefix = "lagr";
  if(kSwitch==0){
    nBPM = 5;
    filename_tag="5bpm";
    eig_treename ="mini_eig";
    lagr_treename ="mini_dit";
    dit_prefix = "dit";
  }else if (kSwitch==1){
    nBPM = 12;
    filename_tag="allbpm";
    eig_treename ="mini_eigall";
    lagr_treename ="mini_lagrall";
  }
  mini_tree = (TTree*)input->Get(eig_treename);
  mini_tree->AddFriend(lagr_treename);
  mini_tree->AddFriend("mini");
  mini_tree->AddFriend("dit=slope",
		       Form("rootfiles/dit_eigslopes_%s_run%d.root",filename_tag.Data(),run));

  vector<TString> DVlist={"asym_us_avg","asym_usl","asym_usr"};
  vector<TString> det_name={"us_avg","usl","usr"};
  Int_t nDV = DVlist.size();
  for(int idv=0;idv<nDV;idv++){
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Print(Form("plots/run%d_%s_compare_%s.pdf[",run,
		   DVlist[idv].Data(),filename_tag.Data()));
    TMultiGraph *mgrms = new TMultiGraph(); // rms
    TLegend *legrms = new TLegend(0.9,0.7,1.0,0.9);
    mini_tree->Draw(Form("eig_%s.rms/ppm:Entry$:run:mini",DVlist[idv].Data()),"","goff");
    Int_t nmini = mini_tree->GetSelectedRows();
    double *y_val = mini_tree->GetV1();
    double *x_val = new double[nmini];
    
    vector<Int_t> fRun;
    vector<Int_t> fMini;
    Double_t *yval = mini_tree->GetV3();
    Double_t *xval = mini_tree->GetV4();

    for(int i =0;i<nmini;i++){
      x_val[i] = i;
      fRun.push_back(yval[i]);
      fMini.push_back(xval[i]);
    }

    TGraph *g1_reg_rms = new TGraph(mini_tree->GetSelectedRows(),
				    x_val,y_val);
    g1_reg_rms->SetMarkerStyle(20);
    g1_reg_rms->SetMarkerColor(kRed);
    g1_reg_rms->SetLineColor(kRed);

    mini_tree->Draw(Form("%s_%s.rms/ppm:Entry$",
			 dit_prefix.Data(),DVlist[idv].Data()),"","goff");
    double *y_val2 = mini_tree->GetV1();
    TGraph *g1_lag_rms = new TGraph(mini_tree->GetSelectedRows(),
				    x_val,y_val2);
    g1_lag_rms->SetMarkerStyle(20);
    g1_lag_rms->SetMarkerColor(kBlue);
    g1_lag_rms->SetLineColor(kBlue);
    mgrms->Add(g1_reg_rms,"lp");
    mgrms->Add(g1_lag_rms,"lp");
    legrms->AddEntry(g1_reg_rms,"Regression");
    
    if(kSwitch==0)
      legrms->AddEntry(g1_lag_rms,"Dithering");
    else
      legrms->AddEntry(g1_lag_rms,"Lagrange");
    
    mgrms->Draw("A");
    // Labeling run numbers; 
    TH1D *hmgrms = (TH1D*)mgrms->GetHistogram();
    hmgrms->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
    Int_t last_run_label=0;
    double ymax  = hmgrms->GetYaxis()->GetXmax();
    double ymin  = hmgrms->GetYaxis()->GetXmin();
    for(int i=0;i<nmini;i++){
      if(fRun[i]!=last_run_label){
	last_run_label = fRun[i];
	hmgrms->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
	TLine *line_buff= new TLine(i,ymin,i,ymax);
	line_buff->SetLineWidth(1);
	line_buff->SetLineStyle(4);
	line_buff->SetLineColor(14);
	line_buff->Draw("same");
      }else
	hmgrms->GetXaxis()->SetBinLabel(i+1,"");
    }
    hmgrms->SetTitle(run_title+"Corrected RMS (ppm): " + DVlist[idv] +";Run Number; RMS(ppm)");
    hmgrms->GetXaxis()->SetTitleOffset(1.25);
    legrms->Draw("same");
    c1->Print(Form("plots/run%d_%s_compare_%s.pdf[",run,
		   DVlist[idv].Data(),filename_tag.Data()));
    /// 
    /// 
    /// Quadrature diff 
    ///
    /// 
    mini_tree->Draw(Form("sqrt(%s_%s.rms**2 - eig_%s.rms**2)/ppm:Entry$",
			 dit_prefix.Data(),DVlist[idv].Data(),DVlist[idv].Data()),"","goff");
    double *diff_y = mini_tree->GetV1();
    TGraph *g1_diff_rms = new TGraph(mini_tree->GetSelectedRows(),
				     x_val,diff_y);
    g1_diff_rms->SetMarkerStyle(20);
    g1_diff_rms->Draw("ALP");
    TH1D *hg1_diff_rms = (TH1D*)g1_diff_rms->GetHistogram();
    hg1_diff_rms->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
    last_run_label=0;
    ymax = hg1_diff_rms->GetYaxis()->GetXmax();
    ymin = hg1_diff_rms->GetYaxis()->GetXmin();
    for(int i=0;i<nmini;i++){
      if(fRun[i]!=last_run_label){
	last_run_label = fRun[i];
	hg1_diff_rms->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
	TLine *line_buff = new TLine(i,ymin,i,ymax);
	line_buff->SetLineWidth(1);
	line_buff->SetLineStyle(4);
	line_buff->SetLineColor(14);
	line_buff->Draw("same");
      }else
	hg1_diff_rms->GetXaxis()->SetBinLabel(i+1,"");
    }
    hg1_diff_rms->SetTitle(run_title+DVlist[idv] + Form(": #sqrt{#sigma^{2}_{%s} -#sigma^{2}_{reg}} (ppm);Run Number; RMS(ppm)",dit_prefix.Data()));
    hg1_diff_rms->GetXaxis()->SetTitleOffset(1.25);
    c1->Print(Form("plots/run%d_%s_compare_%s.pdf",run,
		   DVlist[idv].Data(),filename_tag.Data()));
    
    //
    // Eigenvector Corrections
    //

    TMultiGraph *mg_corr = new TMultiGraph();
    TLegend *leg_corr = new TLegend(0.9,0.6,0.99,0.9);
    Int_t color_offset=0;
    for(int iiv=0;iiv<nBPM;iiv++){
      mini_tree->Draw(Form("fabs( dit.%s_evMon%d * diff_evMon%d.rms )/ppm:Entry$",
			   det_name[idv].Data(),iiv,iiv),"","goff");
      double *y_corr = mini_tree->GetV1();
      double *x_corr = mini_tree->GetV2();
      TGraph *g_ev_corr = new TGraph(mini_tree->GetSelectedRows(),
				     x_corr,y_corr);
      g_ev_corr->SetMarkerStyle(20);
      if(iiv+1==5 || iiv+1==9)
	color_offset ++;
      g_ev_corr->SetMarkerColor(iiv+1+color_offset);
      g_ev_corr->SetLineColor(iiv+1+color_offset);
      mg_corr->Add(g_ev_corr,"lp");
      leg_corr->AddEntry(g_ev_corr,Form("Eigenvector %d",iiv),"lp");
    }
    mg_corr->Draw("A");

    TH1D *hmg_corr = (TH1D*) mg_corr->GetHistogram();
    hmg_corr->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
    Int_t last_run_corr = 0;
    double ymax_corr = hmg_corr->GetYaxis()->GetXmax();
    double ymin_corr = hmg_corr->GetYaxis()->GetXmin();
    for(int i=0;i<nmini;i++){
      if(fRun[i]!=last_run_corr){
	last_run_corr = fRun[i];
	hmg_corr->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
	TLine *line_buff= new TLine(i,ymin_corr,i,ymax_corr);
	line_buff->SetLineWidth(1);
	line_buff->SetLineStyle(4);
	line_buff->SetLineColor(14);
	line_buff->Draw("same");
      }else
	hmg_corr->GetXaxis()->SetBinLabel(i+1,"");
    }
    if(kSwitch==0)
      hmg_corr->SetTitle("Dithering Corrections Width (ppm) by Eigenvectors; Run Number ; (PPM)");
    else 
      hmg_corr->SetTitle("Lagrange Corrections Width (ppm) by Eigenvectors; Run Number ; (PPM)");
    
    leg_corr->Draw("same");
    c1->Print(Form("plots/run%d_%s_compare_%s.pdf",run,
		   DVlist[idv].Data(),filename_tag.Data()));
    
    /// 
    /// Slopes Comparisons
    ///  
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
      mini_tree->Draw(Form("sign%d*%s.%s_evMon%d/(ppm/um):Entry$",
			   i,eig_treename.Data(),det_name.Data(),i),
		      "","goff");
      double *y_val = mini_tree->GetV1();
      TGraph *g1_reg_slope = new TGraph(mini_tree->GetSelectedRows(),
					x_val,y_val);
      g1_reg_slope->SetMarkerStyle(20);
      g1_reg_slope->SetMarkerColor(kRed);
      g1_reg_slope->SetLineColor(kRed);
    
      mini_tree->Draw(Form("sign%d*dit.%s_evMon%d/(ppm/um):Entry$",
			   i,det_name.Data(),i),
		      "","goff");
      double *y_val2 = mini_tree->GetV1();
      TGraph *g1_lag_slope = new TGraph(mini_tree->GetSelectedRows(),
					x_val,y_val2);
      g1_lag_slope->SetMarkerStyle(20);
      g1_lag_slope->SetMarkerColor(kBlue);
      g1_lag_slope->SetLineColor(kBlue);
    
      TMultiGraph *mgslope = new TMultiGraph(); // rms
      TLegend *legslope = new TLegend(0.9,0.7,1.0,0.9);
    
      mgslope->Add(g1_reg_slope,"lp");
      mgslope->Add(g1_lag_slope,"lp");
      legslope->AddEntry(g1_reg_slope,"Regression");
      if(kSwitch==0)
	legslope->AddEntry(g1_lag_slope,"Dithering");
      else
	legslope->AddEntry(g1_lag_slope,"Lagrange");
      mgslope->Draw("A");
      TH1D *hmgslope = (TH1D*)mgslope->GetHistogram();
      hmgslope->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
      Int_t last_run_label = 0;
      double ymax = hmgslope->GetYaxis()->GetXmax();
      double ymin = hmgslope->GetYaxis()->GetXmin();
      for(int i=0;i<nmini;i++){
	if(fRun[i]!=last_run_label){
	  last_run_label = fRun[i];
	  hmgslope->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
	  TLine *line_buff = new TLine(i,ymin,i,ymax);
	  line_buff->SetLineWidth(1);
	  line_buff->SetLineStyle(4);
	  line_buff->SetLineColor(14);
	  line_buff->Draw("same");
	}else
	  hmgslope->GetXaxis()->SetBinLabel(i+1,"");
      }
      hmgslope->SetTitle(run_title+Form("Slope : %s_vs_evMon%d (ppm/um);Run Number; Slope (ppm/um)",
					det_name.Data(),i));
      hmgslope->GetXaxis()->SetTitleOffset(1.25);
      legslope->Draw("same");
      pad2->cd();
      mini_tree->Draw(Form("( fabs(dit.%s_evMon%d)-fabs(%s.%s_evMon%d) ) /(ppm/um)",
			   det_name.Data(),i,
			   eig_treename.Data(),det_name.Data(),i),"","");
      TH1D *h1d = (TH1D*)gPad->FindObject("htemp");
      if(h1d!=NULL){
	h1d->SetTitle("|dit_slope| - |reg_slope| (ppm/um) ");
	h1d->SetTitleSize(0.05);
      }
      c1->Print(Form("plots/run%d_%s_compare_%s.pdf",run,
		     DVlist[idv].Data(),filename_tag.Data()));
    }
    c1->Print(Form("plots/run%d_%s_compare_%s.pdf]",run,
		   DVlist[idv].Data(),filename_tag.Data()));
  }

  ///////////////////////////////
  /////////////////////////////// diff_evMon
  ///////////////////////////////
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Print(Form("plots/run%d_diff_evMon_%s.pdf[",run,filename_tag.Data()));
  c2->cd();
  mini_tree->Draw("run:mini","","goff");
  vector<Int_t> fRun;
  vector<Int_t> fMini;
  Double_t *yval = mini_tree->GetV1();
  Double_t *xval = mini_tree->GetV2();
  Int_t nmini = mini_tree->GetSelectedRows();
  for(int i =0;i<nmini;i++){
    fRun.push_back(yval[i]);
    fMini.push_back(xval[i]);
  }

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
    
    TH1D *hg1_avg = (TH1D*)g1_avg->GetHistogram();
    hg1_avg->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
    Int_t last_run_label =0;
    double ymax = hg1_avg->GetYaxis()->GetXmax();
    double ymin = hg1_avg->GetYaxis()->GetXmin();
    for(int i=0;i<nmini;i++){
      if(fRun[i]!=last_run_label){
	last_run_label = fRun[i];
	hg1_avg->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
	TLine *line_buff = new TLine(i,ymin,i,ymax);
	line_buff->SetLineWidth(1);
	line_buff->SetLineStyle(4);
	line_buff->SetLineColor(14);
	line_buff->Draw("same");
      }else
	hg1_avg->GetXaxis()->SetBinLabel(i+1,"");
    }
    hg1_avg->SetTitle(run_title+ channel+" Mean (nm) ;Run Number; diff (nm)");
    hg1_avg->GetXaxis()->SetTitleOffset(1.25);
    g1_avg->Fit("pol0","Q");
    c2->Print(Form("plots/run%d_diff_evMon_%s.pdf",run,filename_tag.Data()));

    mini_tree->Draw(Form("%s.rms/um:Entry$",channel.Data())
		    ,"","goff");
    double *y_rms = mini_tree->GetV1();
    double *x_mean2 = mini_tree->GetV2();
    TGraph *g1_rms = new TGraph(mini_tree->GetSelectedRows(),
				x_mean2,y_mean);
  
    g1_rms->SetMarkerStyle(20);
    g1_rms->Draw("ALP");
    TH1D *hg1_rms = (TH1D*)g1_rms->GetHistogram();
    hg1_rms->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
    last_run_label=0;
    ymax = hg1_rms->GetYaxis()->GetXmax();
    ymin = hg1_rms->GetYaxis()->GetXmin();
    for(int i=0;i<nmini;i++){
      if(fRun[i]!=last_run_label){
	last_run_label = fRun[i];
	hg1_rms->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
	TLine *line_buff = new TLine(i,ymin,i,ymax);
	line_buff->SetLineWidth(1);
	line_buff->SetLineStyle(4);
	line_buff->SetLineColor(14);
	line_buff->Draw("same");
      }else
	hg1_rms->GetXaxis()->SetBinLabel(i+1,"");
    }
    hg1_rms->SetTitle(run_title+channel+" RMS (um) ;Run Number; RMS (um)");
    hg1_rms->GetXaxis()->SetTitleOffset(1.25);
    g1_rms->Fit("pol0","Q");
    c2->Print(Form("plots/run%d_diff_evMon_%s.pdf",run,filename_tag.Data()));

  }
  c2->Print(Form("plots/run%d_diff_evMon_%s.pdf]",run,filename_tag.Data()));

}

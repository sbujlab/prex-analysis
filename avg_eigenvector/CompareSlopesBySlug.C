void CompareSlopesBySlug(Int_t slug=94){
  gStyle->SetOptFit(1);
  TFile *input = TFile::Open(Form("rootfiles/slug%d_sorted_eigenvector_allbpm.root",slug));
  TTree *eig_tree;
  Int_t nBPM=12;
  if(slug<=2)
    nBPM = 10;
  TString slug_title=Form("Slug%d: ",slug);
  TString filename_tag="allbpm";
  
  eig_tree = (TTree*)input->Get("eig");
  eig_tree->AddFriend("reg");
  eig_tree->AddFriend("lagr");

  eig_tree->Draw("run:mini","","goff");
  vector<Int_t> fRun;
  vector<Int_t> fMini;
  Double_t *yval = eig_tree->GetV1();
  Double_t *xval = eig_tree->GetV2();
  Int_t nmini = eig_tree->GetSelectedRows();
  for(int i =0;i<nmini;i++){
    fRun.push_back(yval[i]);
    fMini.push_back(xval[i]);
  }

  vector<TString> DVlist={"asym_us_avg","asym_usl","asym_usr","asym_us_dd"};
  vector<TString> det_name={"us_avg","usl","usr","us_dd"};
  Int_t nDV = DVlist.size();
  for(int idv=0;idv<nDV;idv++){
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Print(Form("plots/slug%d_%s_compare_%s.pdf[",slug,
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
    TString arm_cut="";
    if( DVlist[idv].Contains("l"))
      arm_cut = "arm_flag!=1";
    if( DVlist[idv].Contains("r"))
      arm_cut = "arm_flag!=2";
    if( DVlist[idv].Contains("avg") || DVlist[idv].Contains("dd") )
      arm_cut = "arm_flag==0";

    for(int i=0;i<nBPM;i++){
      gStyle->SetStatH(0.2);
      gStyle->SetStatW(0.3);
      c1->Clear("D");
      pad1->cd();
      eig_tree->Draw(Form("reg.%s_evMon%d*1e3:Entry$",
			   det_name.Data(),i),
		      arm_cut,"goff");
      double *y_val = eig_tree->GetV1();
      double *x_val = eig_tree->GetV2();
      TGraph *g1_reg_slope = new TGraph(eig_tree->GetSelectedRows(),
					x_val,y_val);
      g1_reg_slope->SetMarkerStyle(20);
      g1_reg_slope->SetMarkerColor(kRed);
      g1_reg_slope->SetLineColor(kRed);
    
      eig_tree->Draw(Form("lagr.%s_evMon%d*1e3:Entry$",
			   det_name.Data(),i),
		      arm_cut,"goff");
      double *y_val2 = eig_tree->GetV1();
      TGraph *g1_lag_slope = new TGraph(eig_tree->GetSelectedRows(),
					x_val,y_val2);
      g1_lag_slope->SetMarkerStyle(20);
      g1_lag_slope->SetMarkerColor(kBlue);
      g1_lag_slope->SetLineColor(kBlue);
    
      TMultiGraph *mgslope = new TMultiGraph(); // rms
      TLegend *legslope = new TLegend(0.9,0.7,1.0,0.9);
    
      mgslope->Add(g1_reg_slope,"lp");
      mgslope->Add(g1_lag_slope,"lp");
      legslope->AddEntry(g1_reg_slope,"Regression");
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
      hmgslope->SetTitle(slug_title+Form("Slope : %s_vs_evMon%d (ppm/um);Run Number; Slope (ppm/um)",
					 det_name.Data(),i));
      hmgslope->GetXaxis()->SetTitleOffset(1.25);
      legslope->Draw("same");
      pad2->cd();
      eig_tree->Draw(Form("( fabs(lagr.%s_evMon%d)-fabs(reg.%s_evMon%d) ) *1e3",
			   det_name.Data(),i,det_name.Data(),i),arm_cut,"");
      TH1D *h1d = (TH1D*)gPad->FindObject("htemp");
      if(h1d!=NULL){
	h1d->SetTitle("|dit_slope| - |reg_slope| (ppm/um) ");
	h1d->SetTitleSize(0.05);
      }
      c1->Print(Form("plots/slug%d_%s_compare_%s.pdf",slug,
		     DVlist[idv].Data(),filename_tag.Data()));
    }
    c1->Print(Form("plots/slug%d_%s_compare_%s.pdf]",slug,
		   DVlist[idv].Data(),filename_tag.Data()));
  }

  ///////////////////////////////
  /////////////////////////////// diff_evMon
  ///////////////////////////////
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Print(Form("plots/slug%d_diff_evMon_%s.pdf[",slug,filename_tag.Data()));
  c2->cd();

  for(int i=0;i<nBPM;i++){
    TString channel = Form("diff_evMon%d",i);
    eig_tree->Draw(Form("%s*1e6:%s.err*1e6:Entry$",channel.Data(),channel.Data())
		    ,"","goff");
    double *y_mean = eig_tree->GetV1();
    double *y_err = eig_tree->GetV2();
    double *x_mean = eig_tree->GetV3();
    TGraphErrors *g1_avg = new TGraphErrors(eig_tree->GetSelectedRows(),
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
    hg1_avg->SetTitle(slug_title+ channel+" Mean (nm) ;Run Number; diff (nm)");
    hg1_avg->GetXaxis()->SetTitleOffset(1.25);
    g1_avg->Fit("pol0","Q");
    c2->Print(Form("plots/slug%d_diff_evMon_%s.pdf",slug,filename_tag.Data()));

    eig_tree->Draw(Form("%s.rms*1e3:Entry$",channel.Data())
		    ,"","goff");
    double *y_rms = eig_tree->GetV1();
    double *x_mean2 = eig_tree->GetV2();
    TGraph *g1_rms = new TGraph(eig_tree->GetSelectedRows(),
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
    hg1_rms->SetTitle(slug_title+channel+" RMS (um) ;Run Number; RMS (um)");
    hg1_rms->GetXaxis()->SetTitleOffset(1.25);
    g1_rms->Fit("pol0","Q");
    c2->Print(Form("plots/slug%d_diff_evMon_%s.pdf",slug,filename_tag.Data()));

  }
  c2->Print(Form("plots/slug%d_diff_evMon_%s.pdf]",slug,filename_tag.Data()));

}

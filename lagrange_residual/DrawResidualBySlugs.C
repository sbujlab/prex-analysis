/*
  author: Tao Ye
  Last update: Aug 2019
*/

void DrawResidualBySlugs(){
  gStyle->SetOptStat(111111);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.35);
  gStyle->SetTitleFontSize(0.07);


  TString det_name[]={"usl","usr","us_avg","us_dd"};
  TString redundant_cut[7]={ "(run>=4735)","(run<4439 || run>4450)",
			     "(0)","(0)",
			     "(run<4745)","(run>=4439 && run<=4450)",
			     "(0)"};

  for(int idet =0;idet<4;idet++){
    TString arm_cut ="";
    if(det_name[idet].Contains("l"))
      arm_cut = "arm_flag != 1";
    if(det_name[idet].Contains("r"))
      arm_cut = "arm_flag != 2";
    if(det_name[idet].Contains("avg") || det_name[idet].Contains("dd"))
      arm_cut = "arm_flag==0";

    TCanvas *c1 = new TCanvas("c1","c1",1200,1000);
    c1->Divide(1,2);
    c1->cd(1);
    TPad *pad1 = new TPad("pad1","pad1",0,0,1,0.4);
    pad1->SetTopMargin(0.0);
    TPad *pad2 = new TPad("pad2","pad2",0,0.4,1,1);
    pad2->SetTopMargin(0.2);
    pad2->SetBottomMargin(0.0);
    pad1->Draw();
    pad2->Draw();
    pad1->SetGridx();
    pad1->SetGridy();
    c1->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();

    c1->Print(Form("./plots/compare_residual_%s_all_slugs.pdf[",det_name[idet].Data()));
    for(int icoil=1;icoil<=7;icoil++){
      double reg_val[94];
      double reg_val_error[94];
      double reg_rms[94];
      double reg_red_rms[94];

      double lagr_val[94];
      double lagr_val_error[94];
      
      double lagr_red_val[94];
      double lagr_red_val_error[94];
      
      double lagr_rms[94];
      double lagr_red_rms[94];
      
      double slug_value[94];
      double slug_red_value[94];
      double slug_reg_value[94];
      
      double abs_val_diff[94];
      double slug_abs_value[94];
      double abs_val_diff_below0[94];
      double slug_abs_value_below0[94];
      
      int counts= 0;
      int red_counts = 0;
      int reg_counts = 0;
      int abs_counts = 0;
      int abs_below0_counts = 0;
      TH1D *h1dptr;
      for(int islug=1;islug<=94;islug++){
	if(islug==7)
	  continue;
	TFile* reg_file = TFile::Open(Form("rootfiles/slug%d_regression_residual.root",islug));
	TFile* lagr_file = TFile::Open(Form("rootfiles/slug%d_lagrange_residual.root",islug));
	TTree *reg_res = (TTree*) reg_file->Get("res");
	TTree *lagr_res = (TTree*) lagr_file->Get("res");
	TString slug_cut = Form("&& slug==%d",islug);
	TString coil_flag = Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil);
	TString cut_text = arm_cut + slug_cut +coil_flag;
	int npt1=reg_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			       cut_text,"goff");
	if(npt1!=0){
	  reg_counts ++;
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  reg_val[reg_counts-1] = h1dptr->GetMean();
	  reg_rms[reg_counts-1] = h1dptr->GetRMS();
	  double nsamples = h1dptr->GetEntries();
	  reg_val_error[reg_counts-1] = reg_rms[reg_counts-1]/ sqrt(nsamples);
	  slug_reg_value[reg_counts-1] = islug;
	}

	npt1=lagr_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			    cut_text+"&&!"+redundant_cut[icoil-1],"goff");
	if(npt1!=0){
	  counts++;
	  slug_value[counts-1] = islug;
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  lagr_val[counts-1] = h1dptr->GetMean();
	  lagr_rms[counts-1] = h1dptr->GetRMS();
	  double nsamples =  h1dptr->GetEntries();
	  lagr_val_error[counts-1 ]= lagr_rms[counts-1]/sqrt(nsamples);
	}
	int npt2=lagr_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			 cut_text+"&&"+redundant_cut[icoil-1],"goff");
		  
	if(npt2!=0){
	  red_counts++;
	  slug_red_value[red_counts-1] = islug;
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  lagr_red_val[red_counts-1] = h1dptr->GetMean();
	  lagr_red_rms[red_counts-1] = h1dptr->GetRMS();
	  double nsamples = h1dptr->GetEntries();
	  lagr_red_val_error[red_counts-1] = lagr_red_rms[red_counts-1]/sqrt(nsamples);
	}

	if(npt2!=0 || npt1!=0){
	  abs_counts++;
	  lagr_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			 cut_text,"goff");
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  double val1=h1dptr->GetMean();
	  reg_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			cut_text,"goff");
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  double val2=h1dptr->GetMean();
	  abs_val_diff[abs_counts-1] = fabs(val2)- fabs(val1);
	  slug_abs_value[abs_counts-1] = islug;
	  if( abs_val_diff[abs_counts-1] <0){
	    abs_below0_counts++;
	    abs_val_diff_below0[abs_below0_counts-1] = fabs(val2)- fabs(val1);
	    slug_abs_value_below0[abs_below0_counts-1] = islug;
	  }
	    
	}
	
	reg_file->Close();
	lagr_file->Close();
      } // end slug loop
      TMultiGraph *mg_val = new TMultiGraph();
      TMultiGraph *mg_rms = new TMultiGraph();
      TLegend *leg_val = new TLegend(0.9,0.7,1.0,0.9);
      TLegend *leg_rms = new TLegend(0.9,0.7,1.0,0.9);
      TLegend *leg_abs = new TLegend(0.9,0.7,1.0,0.9);
      TMultiGraph *mg_abs_diff = new TMultiGraph();
      TGraph *g_abs_diff = new TGraph(abs_counts, slug_abs_value, abs_val_diff);
      TGraph *g_abs_diff_below0 = new TGraph(abs_below0_counts, slug_abs_value_below0, abs_val_diff_below0);
      g_abs_diff->SetMarkerStyle(20);
      g_abs_diff_below0->SetMarkerStyle(23);
      g_abs_diff_below0->SetMarkerSize(1.2);
      g_abs_diff_below0->SetMarkerColor(kMagenta);
      mg_abs_diff->Add(g_abs_diff,"lp");
      mg_abs_diff->Add(g_abs_diff_below0,"p");
      leg_abs->AddEntry(g_abs_diff_below0, "Below Zero","p");
      TGraphErrors *g_reg_val = new TGraphErrors(reg_counts, slug_reg_value, reg_val, 0, reg_val_error);
      TGraphErrors *g_lagr_val = new TGraphErrors(counts, slug_value, lagr_val, 0 , lagr_val_error);
      TGraphErrors *g_lagr_red_val = new TGraphErrors(red_counts, slug_red_value, lagr_red_val, 0 , lagr_red_val_error);

      pad1->cd();
      mg_abs_diff->Draw("A");
      mg_abs_diff->GetXaxis()->SetRangeUser(0,100);
      mg_abs_diff->GetXaxis()->SetTitleSize(0.05);
      mg_abs_diff->GetYaxis()->SetTitleSize(0.05);
      mg_abs_diff->GetXaxis()->SetLabelSize(0.08);
      mg_abs_diff->GetYaxis()->SetLabelSize(0.08);
      mg_abs_diff->SetTitle("|Reg| -|Lagr|; Slug ; (ppm/um)");
      leg_abs->Draw("same");
      
      g_reg_val->SetMarkerStyle(20);
      g_reg_val->SetMarkerColor(kRed);
      g_reg_val->SetLineColor(kRed);
      g_lagr_val->SetMarkerStyle(20);
      g_lagr_val->SetMarkerColor(kBlue);
      g_lagr_val->SetLineColor(kBlue);
      g_lagr_red_val->SetMarkerStyle(4);
      g_lagr_red_val->SetMarkerColor(kBlue);
      g_lagr_red_val->SetLineColor(kBlue);
      mg_val->Add(g_reg_val,"lp");
      mg_val->Add(g_lagr_val,"lp");
      mg_val->Add(g_lagr_red_val,"lp");
      leg_val->AddEntry(g_reg_val, " Regression ","lp");
      leg_val->AddEntry(g_lagr_val, " Lagrange ","lp");
      leg_val->AddEntry(g_lagr_red_val, " Redundant ","lp");
      
      pad2->cd();
      mg_val->Draw("A");
      leg_val->Draw("same");
      mg_val->SetTitle(Form("%s vs coil %d residual sensitivity(ppm/count);slug;Residual Sensitivity(ppm/count)",det_name[idet].Data(),icoil));
      mg_val->GetXaxis()->SetRangeUser(0,100);
      mg_val->GetXaxis()->SetLabelSize(0.05);
      mg_val->GetYaxis()->SetLabelSize(0.05);
      mg_val->GetYaxis()->SetTitleSize(0.05);
      
      TGraph *g_reg_rms = new TGraph(reg_counts, slug_reg_value, reg_rms);
      TGraph *g_lagr_rms = new TGraph(counts, slug_value, lagr_rms);
      TGraph *g_lagr_red_rms = new TGraph(red_counts, slug_red_value, lagr_red_rms);
      g_reg_rms->SetMarkerStyle(20);
      g_reg_rms->SetMarkerColor(kRed);
      g_reg_rms->SetLineColor(kRed);
      g_lagr_rms->SetMarkerStyle(20);
      g_lagr_rms->SetMarkerColor(kBlue);
      g_lagr_rms->SetLineColor(kBlue);
      g_lagr_red_rms->SetMarkerStyle(4);
      g_lagr_red_rms->SetMarkerColor(kBlue);
      g_lagr_red_rms->SetLineColor(kBlue);
      mg_rms->Add(g_reg_rms,"lp");
      mg_rms->Add(g_lagr_rms,"lp");
      mg_rms->Add(g_lagr_red_rms,"lp");
      leg_rms->AddEntry(g_reg_rms, " Regression ","lp");
      leg_rms->AddEntry(g_lagr_rms, " Lagrange ","lp");
      leg_rms->AddEntry(g_lagr_red_rms, " Redundant ","lp");
      c1->cd(2);
      mg_rms->Draw("A");
      mg_rms->GetXaxis()->SetRangeUser(0,100);
      leg_rms->Draw("same");
      mg_rms->SetTitle(Form("%s vs coil %d residual Width (ppm/count);slug;Residual Width(ppm/count)",det_name[idet].Data(),icoil));

      c1->Print(Form("./plots/compare_residual_%s_all_slugs.pdf",det_name[idet].Data()));
    } // end of coil loop
    c1->Print(Form("./plots/compare_residual_%s_all_slugs.pdf]",det_name[idet].Data()));

  } // end of detector loop
}

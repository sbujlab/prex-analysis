/*
  author: Tao Ye
  Last update: Aug 2019
*/

void DrawResidualBySlugs_v2(){
  gStyle->SetOptStat(111111);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.35);
  gStyle->SetTitleFontSize(0.07);


  TString det_name[]={"usl","usr"};
  TString redundant_cut[7]={ "(run>=4735)","(run<4439 || run>4450)",
			     "(0)","(0)",
			     "(run<4745)","(run>=4439 && run<=4450)",
			     "(0)"};
  int ndet = sizeof(det_name)/sizeof(*det_name);
  for(int idet =0;idet<ndet;idet++){
    TString arm_cut ="";
    if(det_name[idet].Contains("l"))
      arm_cut = "arm_flag != 1";
    if(det_name[idet].Contains("r"))
      arm_cut = "arm_flag != 2";
    if(det_name[idet].Contains("avg") || det_name[idet].Contains("dd"))
      arm_cut = "arm_flag==0";

    TCanvas *c1 = new TCanvas("c1","c1",1200,1000);
    c1->Divide(1,2);
    c1->Print(Form("./plots/compare_residual_%s_all_slugs_v2.pdf[",det_name[idet].Data()));
    for(int icoil=1;icoil<=7;icoil++){
      double reg_val[94];
      double reg_val_error[94];

      double lagr_val[94];
      double lagr_val_error[94];
      double lagr_red_val[94];
      double lagr_red_val_error[94];


      double lagr_rms[94];
      double lagr_red_rms[94];
      double reg_rms[94];
      double reg_red_rms[94];
      
      double slug_value[94];
      double slug_red_value[94];
      double slug_reg_value[94];
      
      double raw_val[94];
      double raw_rms[94];
      double raw_val_error[94];
      double slug_raw_value[94];

      int counts= 0;
      int red_counts = 0;
      int reg_counts = 0;
      int raw_counts = 0;
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
	
	int npt_raw =reg_res->Draw(Form("%s_coil%d_sens*1e6",det_name[idet].Data(),icoil),
				   cut_text,"goff");
	if(npt_raw!=0){
	  raw_counts ++;
	  if(h1dptr!=NULL){
	    h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	    raw_val[raw_counts-1] = h1dptr->GetMean();
	    raw_rms[raw_counts-1] = h1dptr->GetRMS();
	    double nsamples = h1dptr->GetEntries();
	    raw_val_error[raw_counts-1] = raw_rms[raw_counts-1]/ sqrt(nsamples);
	    slug_raw_value[raw_counts-1] = islug;
	  }
	}
	
	int npt1=reg_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			       cut_text,"goff");
	if(npt1!=0){
	  reg_counts ++;
	  if(h1dptr!=NULL){
	    h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	    reg_val[reg_counts-1] = h1dptr->GetMean();
	    reg_rms[reg_counts-1] = h1dptr->GetRMS();
	    double nsamples = h1dptr->GetEntries();
	    reg_val_error[reg_counts-1] = reg_rms[reg_counts-1]/ sqrt(nsamples);
	    slug_reg_value[reg_counts-1] = islug;
	  }
	}

	npt1=lagr_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			    cut_text+"&&!"+redundant_cut[icoil-1],"goff");
	if(npt1!=0){
	  counts++;
	  slug_value[counts-1] = islug;
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  if(h1dptr!=NULL){
	    lagr_val[counts-1] = h1dptr->GetMean();
	    lagr_rms[counts-1] = h1dptr->GetRMS();
	    double nsamples =  h1dptr->GetEntries();
	    lagr_val_error[counts-1 ]= lagr_rms[counts-1]/sqrt(nsamples);
	  }
	}
	int npt2=lagr_res->Draw(Form("%s_coil%d_res*1e6",det_name[idet].Data(),icoil),
			 cut_text+"&&"+redundant_cut[icoil-1],"goff");
		  
	if(npt2!=0){
	  red_counts++;
	  slug_red_value[red_counts-1] = islug;
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  if(h1dptr!=NULL){
	    lagr_red_val[red_counts-1] = h1dptr->GetMean();
	    lagr_red_rms[red_counts-1] = h1dptr->GetRMS();
	    double nsamples = h1dptr->GetEntries();
	    lagr_red_val_error[red_counts-1] = lagr_red_rms[red_counts-1]/sqrt(nsamples);
	  }
	}
	
	reg_file->Close();
	lagr_file->Close();
      } // end slug loop
      c1->cd(1);
      TMultiGraph *mg_val = new TMultiGraph();

      TLegend *leg_val = new TLegend(0.9,0.7,1.0,0.9);

      TGraphErrors *g_raw_val = new TGraphErrors(raw_counts, slug_raw_value, raw_val, 0, raw_val_error);

      TGraphErrors *g_reg_val = new TGraphErrors(reg_counts, slug_reg_value, reg_val, 0, reg_val_error);
      TGraphErrors *g_lagr_val = new TGraphErrors(counts, slug_value, lagr_val, 0 , lagr_val_error);
      TGraphErrors *g_lagr_red_val = new TGraphErrors(red_counts, slug_red_value, lagr_red_val, 0 , lagr_red_val_error);

      g_raw_val->SetMarkerStyle(20);
      g_raw_val->SetMarkerColor(kBlack);
      g_raw_val->SetLineColor(kBlack);

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
      mg_val->SetTitle(Form("%s vs coil %d residual sensitivity(ppm/count);Slug Number; Sensitivity(ppm/count)",det_name[idet].Data(),icoil));
      mg_val->Draw("A");
      leg_val->Draw("same");
      mg_val->GetXaxis()->SetLimits(0,100);
      mg_val->GetXaxis()->SetLabelSize(0.04);
      mg_val->GetXaxis()->SetTitleSize(0.05);

      mg_val->GetYaxis()->SetLabelSize(0.05);
      mg_val->GetYaxis()->SetTitleSize(0.05);

      c1->cd(2);
      g_raw_val->Draw("APL");
      g_raw_val->SetTitle(Form("%s vs coil %d Raw sensitivity(ppm/count);Slug Number; Sensitivity(ppm/count)",det_name[idet].Data(),icoil));
      g_raw_val->GetXaxis()->SetRangeUser(0,100);
      g_raw_val->GetXaxis()->SetLabelSize(0.04);
      g_raw_val->GetXaxis()->SetTitleSize(0.05);

      g_raw_val->GetYaxis()->SetLabelSize(0.05);
      g_raw_val->GetYaxis()->SetTitleSize(0.05);

      c1->Print(Form("./plots/compare_residual_%s_all_slugs_v2.pdf",det_name[idet].Data()));
    } // end of coil loop
    c1->Print(Form("./plots/compare_residual_%s_all_slugs_v2.pdf]",det_name[idet].Data()));

  } // end of detector loop
}

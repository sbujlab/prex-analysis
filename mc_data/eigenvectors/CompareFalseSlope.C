void CompareFalseSlope(Int_t run){
  TFile *input = TFile::Open(Form("rootfiles/mc_lagrange_%d.000.root",run));
  TTree *mini_tree = (TTree*)input->Get("mini");

  TFile *input_normal = TFile::Open(Form("rootfiles/dit_eigslopes_allbpm_run%d.root",run));
  TFile *input_false = TFile::Open(Form("rootfiles/dit_eigslopes_allbpm_false_run%d.root",run));
  TTree *slope_true = (TTree*)input_normal->Get("slope");
  TTree *slope_false = (TTree*)input_false->Get("slope");
  slope_true->AddFriend(mini_tree);
  slope_false->AddFriend(mini_tree);
  mini_tree->AddFriend(slope_false,"false");
  mini_tree->AddFriend(slope_true,"true");
  vector<TString> DVlist={"asym_us_avg","asym_usl","asym_usr"};
  vector<TString> det_name={"us_avg","usl","usr"};
  TString run_title = Form("Run %d:   ", run);
  Int_t nBPM=12;
  Int_t nDV = DVlist.size();
  for(int idv=0;idv<nDV;idv++){
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Print(Form("plots/run%d_%s_compare_false.pdf[",run, DVlist[idv].Data()));
    
    TPad *pad1  = new TPad("pad1","pad1",0.0,0.0,0.7,1.0);
    TPad *pad2  = new TPad("pad2","pad2",0.7,0.0,1.0,1.0);
    pad1->Draw();
    pad2->Draw();
    mini_tree->Draw("run:mini","","goff");
    Int_t nmini = mini_tree->GetSelectedRows();
    double *yval = mini_tree->GetV1();
    double *xval = mini_tree->GetV2();
    double *x_val = new double[nmini];
    vector<Int_t> fRun;
    vector<Int_t> fMini;
    for(int i =0;i<nmini;i++){
      x_val[i] = i;
      fRun.push_back(yval[i]);
      fMini.push_back(xval[i]);
    }

    TString det_name = DVlist[idv];
    det_name.ReplaceAll("asym_","");
    for(int i=0;i<nBPM;i++){
      gStyle->SetStatH(0.2);
      gStyle->SetStatW(0.3);
      c1->Clear("D");
      pad1->cd();
      slope_true->Draw(Form("sign%d*%s_evMon%d/(ppm/um):Entry$",i,det_name.Data(),i),
			 "","goff");
      double *y_val = slope_true->GetV1();
      TGraph *g1_norm_slope = new TGraph(slope_true->GetSelectedRows(),
					 x_val,y_val);
      g1_norm_slope->SetMarkerStyle(20);
      g1_norm_slope->SetMarkerColor(kRed);
      g1_norm_slope->SetLineColor(kRed);
    
      slope_false->Draw(Form("sign%d*%s_evMon%d/(ppm/um):Entry$",i,det_name.Data(),i),
			"","goff");
      double *y_val2 = slope_false->GetV1();
      TGraph *g1_false_slope = new TGraph(slope_false->GetSelectedRows(),
					  x_val,y_val2);
      g1_false_slope->SetMarkerStyle(20);
      g1_false_slope->SetMarkerColor(kBlue);
      g1_false_slope->SetLineColor(kBlue);
    
      TMultiGraph *mgslope = new TMultiGraph(); // rms
      TLegend *legslope = new TLegend(0.9,0.7,1.0,0.9);
      mgslope->Add(g1_norm_slope,"lp");
      mgslope->Add(g1_false_slope,"lp");
      legslope->AddEntry(g1_norm_slope,"Lagrange Normal");
      legslope->AddEntry(g1_false_slope,"Lagrange False");
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
      mini_tree->Draw(Form("( fabs(true.%s_evMon%d)-fabs(false.%s_evMon%d) ) /(ppm/um)",
			   det_name.Data(),i,det_name.Data(),i),"","");
      TH1D *h1d = (TH1D*)gPad->FindObject("htemp");
      if(h1d!=NULL){
	h1d->SetTitle("|true_slope| - |false_slope| (ppm/um) ");
	h1d->SetTitleSize(0.05);
      }
      c1->Print(Form("plots/run%d_%s_compare_false.pdf",run,DVlist[idv].Data()));
    }


    c1->Print(Form("plots/run%d_%s_compare_false.pdf]",run,DVlist[idv].Data()));
  }
}

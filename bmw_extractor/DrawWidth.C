void DrawWidth(){
  TFile* input = TFile::Open("bmw_plus_allslugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini");
  TTree *info_tree = (TTree*)input->Get("mini_info");
  mini_tree->AddFriend(info_tree);

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  c1->cd(2);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  mini_tree->SetMarkerStyle(8);
  TString forced_cut = "( (run==3450&&mini==1) || ( run==3803&&mini==0) || (run==3840&&mini==7)|| (run==3853&&mini==2) || (run==4174&&mini==0) || (run==4174&&mini==1) || (run==4339&&mini==4) || run==3342 || run==3370 )";
  
  TString title[] = {"208Pb2", "208Pb10", "208Pb9", "208Pb8",
		     "208Pb5","208Pb7","208Pb6"};
  TString slug_cut[] = {" && slug>=12 && slug <=13 ",
			"&& slug>=26 && slug<=28",
			"&& slug>=44 && slug<=45",
			"&& slug>=58 && slug<=59",
			" && slug>=68 && slug<=71 && run!=4614",
			" && slug>=86 && slug<=88",
			"&& slug>=93 && slug<=94"};
  Int_t ncut = sizeof(slug_cut)/sizeof(*slug_cut);
  for(int ic=0;ic<ncut;ic++){
    c1->Clear("D");
    c1->cd(1);
    mini_tree->Draw("reg_asym_us_avg.rms*1e6:Entry$:yield_bcm_target",
		    "arm_flag==0 && kGood"+slug_cut[ic],"COLZ");
    TH2F* h2dptr = (TH2F*)gPad->FindObject("htemp");
    h2dptr->SetTitle(" Regression(all-bpm) Corrected asym_us_avg RMS(ppm) vs Run ; Run Number;RMS(ppm);Beam Current(uA)");
    
    int npt_cut = mini_tree->Draw("Entry$",
				  "arm_flag==0 && kGood==2"+slug_cut[ic],"goff");
    double *fEntry_cut = mini_tree->GetV1();
    double cut_up;
    double cut_low;
    if(npt_cut!=0){
      cut_up = fEntry_cut[npt_cut-1];
      cut_low = fEntry_cut[0];
    }

    Int_t npt=mini_tree->Draw("run:Entry$","arm_flag==0 && kGood" +slug_cut[ic],"goff");
    double *fRun = mini_tree->GetV1();
    double *fEntry = mini_tree->GetV2();
    double x_min = fEntry[0];
    double x_max = fEntry[npt-1];
    double y_max = h2dptr->GetYaxis()->GetXmax();
    double y_min = h2dptr->GetYaxis()->GetXmin();
    int nbins = x_max - x_min;
    h2dptr->GetXaxis()->Set(nbins,x_min-0.5,x_max-0.5);
    h2dptr->GetXaxis()->SetTickLength(0.0);
    for(int ibin=1;ibin<=nbins;ibin++)
      h2dptr->GetXaxis()->SetBinLabel(ibin, "");

    h2dptr->GetXaxis()->SetLabelSize(0.05);
    h2dptr->GetXaxis()->SetTitleSize(0.03);
    Int_t last_run = -1;
    for(int ipt=0;ipt<npt;ipt++){
      if(fRun[ipt]!=last_run){
	last_run = fRun[ipt];
	h2dptr->GetXaxis()->SetBinLabel(fEntry[ipt]-x_min+1, Form("%d",(int)fRun[ipt]));
	TLine *line = new TLine(fEntry[ipt],y_min,fEntry[ipt],y_max);
	line->SetLineWidth(1);
	line->SetLineStyle(7);
	line->SetLineColorAlpha(14,0.5);
	line->Draw("same");
      }
    }
    c1->cd(1);
    if(npt_cut!=0){
      TBox *box1 = new TBox(cut_low, y_min,cut_up, y_max);
      box1->SetFillColorAlpha(kGray,0.5);
      box1->Draw();
    }

    TPaveText *pave_text = new TPaveText(fEntry[0],y_max-(y_max-y_min)*0.2,
					 fEntry[npt/3],y_max);
    pave_text->AddText(title[ic]);
    pave_text->Draw();

    c1->cd(2);
    mini_tree->Draw("yield_usl*1e3:Entry$","arm_flag==0 && kGood"+slug_cut[ic],"goff");
    TGraph* gusl = new TGraph(mini_tree->GetSelectedRows(),
			      mini_tree->GetV2(),mini_tree->GetV1());
    
    mini_tree->Draw("yield_usr*1e3:Entry$:run","arm_flag==0 && kGood"+slug_cut[ic],"goff");
    double *fEntry_det = mini_tree->GetV2();
    double *fRun_det = mini_tree->GetV3();
    TGraph* gusr = new TGraph(mini_tree->GetSelectedRows(),
			      mini_tree->GetV2(),mini_tree->GetV1());
    gusl->SetMarkerStyle(20);
    gusl->SetMarkerColor(kGreen);
    gusr->SetMarkerStyle(20);
    gusr->SetMarkerColor(kMagenta);
    TMultiGraph *mgdet = new TMultiGraph();
    mgdet->Add(gusl,"P");
    mgdet->Add(gusr,"P");
    mgdet->Draw("AP");
    mgdet->SetTitle("Detector Yield; Run Number ; Normalized yield(mV/uA)");
    TLegend *legdet = new TLegend(0.85,0.7,0.99,0.9);
    legdet->AddEntry(gusl,"USL","p");
    legdet->AddEntry(gusr,"USR","p");
    legdet->Draw("same");

    TH1D *hmgdet = (TH1D*)mgdet->GetHistogram();
    hmgdet->GetXaxis()->Set(nbins,x_min-0.5,x_max-0.5);
    y_max = mgdet->GetYaxis()->GetXmax();
    y_min = mgdet->GetYaxis()->GetXmin();
    mgdet->GetXaxis()->SetTickLength(0.0);
    for(int ipt=0;ipt<npt;ipt++){
      if(fRun_det[ipt]!=last_run){
	last_run = fRun_det[ipt];
	hmgdet->GetXaxis()->SetBinLabel(fEntry_det[ipt]-x_min+1, Form("%d",(int)fRun_det[ipt]));
	TLine *line = new TLine(fEntry[ipt],y_min,fEntry[ipt],y_max);
	line->SetLineWidth(1);
	line->SetLineStyle(7);
	line->SetLineColorAlpha(14,0.5);
	line->Draw("same");
      }
    }
    mgdet->GetXaxis()->SetLabelSize(0.05);
    mgdet->GetXaxis()->SetTitleSize(0.03);
    mgdet->GetYaxis()->SetLabelSize(0.05);
    mgdet->GetYaxis()->SetTitleSize(0.05);
    
    if(npt_cut!=0){
      TBox *box1 = new TBox(cut_low, y_min,cut_up, y_max);
      box1->SetFillColorAlpha(kGray,0.5);
      box1->Draw();
    }

    c1->SaveAs(title[ic]+".pdf");
  }
  
}

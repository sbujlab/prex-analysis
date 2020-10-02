void DrawWidth(){
  TFile* input = TFile::Open("bmw_plus_allslugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini");
  TTree *info_tree = (TTree*)input->Get("mini_info");
  mini_tree->AddFriend(info_tree);

  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->SetRightMargin(0.15);
  c1->SetLeftMargin(0.1);
  mini_tree->SetMarkerStyle(8);
  TString forced_cut = "( (run==3450&&mini==1) || ( run==3803&&mini==0) || (run==3840&&mini==7)|| (run==3853&&mini==2) || (run==4174&&mini==0) || (run==4174&&mini==1) || (run==4339&&mini==4) || run==3342 || run==3370 )";
  TString title[] = {"208Pb2", "208Pb10", "208Pb9", "208Pb8"};
  TString slug_cut[] = {" && slug>=12 && slug <=13 ",
			"&& slug>=26 && slug<=28",
			"&& slug>=44 && slug<=45",
			"&& slug>=58 && slug<=59"};
  for(int ic=0;ic<4;ic++){
    c1->cd();
    c1->Clear();
    mini_tree->Draw("reg_asym_us_avg.rms*1e6:Entry$:yield_bcm_target",
		    "arm_flag==0 && kGood"+slug_cut[ic],"COLZ");
    TH2F* h2dptr = (TH2F*)gPad->FindObject("htemp");
    h2dptr->SetTitle(" Regression(all-bpm) Corrected asym_us_avg RMS(ppm) vs Run ; Run Number;RMS(ppm);Beam Current(uA)");
    
    int npt_cut = mini_tree->Draw("Entry$",
				  "arm_flag==0 && kGood==2"+slug_cut[ic],"goff");
    double *fEntry_cut = mini_tree->GetV1();
    double cut_up = fEntry_cut[npt_cut-1];
    double cut_low = fEntry_cut[0];


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
    c1->cd();
        
    TBox *box1 = new TBox(cut_low, y_min,cut_up, y_max);
    box1->SetFillColorAlpha(kGray,0.5);
    box1->Draw();

    TPaveText *pave_text = new TPaveText(fEntry[0],y_max-(y_max-y_min)*0.2,
					 fEntry[npt/3],y_max);
    pave_text->AddText(title[ic]);
    pave_text->Draw();

    c1->SaveAs(title[ic]+".pdf");
  }
  
}

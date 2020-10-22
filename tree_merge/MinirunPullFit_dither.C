void MinirunPullFit_dither(){
  TFile *input = TFile::Open("./rootfiles/MergedSum_all_slugs.root");
  TTree *mini_tree = (TTree*)input->Get("burst_mulc_dit");
  mini_tree->AddFriend("burst_mulc_dit_combo");
  mini_tree->AddFriend("mini_info","./rootfiles/runinfo_all_slugs.root");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,500);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  mini_tree->Draw("spin*dit_asym_us_avg*1e9: dit_asym_us_avg.err*1e9:run+mini/10.0",
		  "arm_flag==0 && kGood==1","goff");
  TGraphErrors* mger_both = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_both->SetMarkerStyle(7);
  mini_tree->Draw("spin*dit_asym_usr*1e9: dit_asym_usr.err*1e9:run+mini/10.0",
		  "arm_flag==1 && kGood==1","goff");
  TGraphErrors* mger_usr = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_usr->SetMarkerStyle(7);

  mini_tree->Draw("spin*dit_asym_usl*1e9: dit_asym_usl.err*1e9:run+mini/10.0",
		  "arm_flag==2 && kGood==1","goff");
  TGraphErrors* mger_usl = new TGraphErrors(mini_tree->GetSelectedRows(),
					    mini_tree->GetV3(),mini_tree->GetV1(),
					    0,mini_tree->GetV2());
  mger_usl->SetMarkerStyle(7);

  TMultiGraph *mgall = new TMultiGraph();
  mgall->Add(mger_both,"p");
  mgall->Add(mger_usr,"p");
  mgall->Add(mger_usl,"p");
  
  mgall->SetTitle("PREX Dithering Corrected Asymmetry(ppb) ; run; (ppb)"); 
  gStyle->SetOptFit(1);
  mgall->Draw("AP");
  mgall->Fit("pol0");
  double mean_fit = mgall->GetFunction("pol0")->GetParameter(0);


  c2->cd();
  double max_val = mgall->GetYaxis()->GetXmax();
  double min_val = mgall->GetYaxis()->GetXmin();
  TH1D* hpull = new TH1D("hpull","Pull Fit", 150,-6,6);
  // gStyle->SetOptStat(0);
  // gStyle->SetOptFit(0);
  mini_tree->Draw(Form( "(spin*dit_asym_us_avg*1e9 - %f)/(dit_asym_us_avg.err*1e9)>>hpull",mean_fit),
		  "arm_flag==0 && kGood==1","goff");
  mini_tree->Draw(Form( "(spin*dit_asym_usr*1e9 - %f)/(dit_asym_usr.err*1e9) >>+hpull",mean_fit),
		  "arm_flag==1 && kGood==1","goff");
  mini_tree->Draw(Form( "(spin*dit_asym_usl*1e9 - %f)/(dit_asym_usl.err*1e9) >>+hpull",mean_fit),
		  "arm_flag==2 && kGood==1","goff");
  hpull->Draw("PE");
  hpull->SetMarkerStyle(7);
  hpull->SetLineColor(kBlack);
  hpull->SetLineWidth(1);
  hpull->Fit("gaus");
  hpull->GetFunction("gaus")->SetNpx(200);
  hpull->SetTitle(";Normalized difference from average; Number of mini-runs");
  hpull->SetTitleFont(22);
  hpull->SetLabelFont(22);
  hpull->GetYaxis()->SetTitleFont(22);
  hpull->GetYaxis()->SetLabelFont(22);

  c1->SaveAs("prex_dither_minirun_pullfit.pdf");
  c2->SaveAs("prex_dither_minirun_pull_distribution.pdf");
}

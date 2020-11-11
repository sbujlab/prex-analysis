void MinirunPullFit_Lagrange(){
  TFile *input = TFile::Open("./rootfiles/MergedLagrange_all_slugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini_lagrall");
  mini_tree->AddFriend("mini_info","./rootfiles/runinfo_all_slugs.root");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,500);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  mini_tree->Draw("spin*lagr_asym_us_avg*1e9: lagr_asym_us_avg.err*1e9:run+mini/10.0",
		  "arm_flag==0 && kGood==1","goff");
  TGraphErrors* mger_both = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_both->SetMarkerStyle(7);
  mini_tree->Draw("spin*lagr_asym_usr*1e9: lagr_asym_usr.err*1e9:run+mini/10.0",
		  "arm_flag==1 && kGood==1","goff");
  TGraphErrors* mger_usr = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_usr->SetMarkerStyle(7);

  mini_tree->Draw("spin*lagr_asym_usl*1e9: lagr_asym_usl.err*1e9:run+mini/10.0",
		  "arm_flag==2 && kGood==1","goff");
  TGraphErrors* mger_usl = new TGraphErrors(mini_tree->GetSelectedRows(),
					    mini_tree->GetV3(),mini_tree->GetV1(),
					    0,mini_tree->GetV2());
  mger_usl->SetMarkerStyle(7);

  TMultiGraph *mgall = new TMultiGraph();
  mgall->Add(mger_both,"p");
  mgall->Add(mger_usr,"p");
  mgall->Add(mger_usl,"p");
  
  mgall->SetTitle("PREX Lagrange Corrected Asymmetry(ppb) ; run; (ppb)"); 
  gStyle->SetOptFit(1);
  mgall->Draw("AP");
  mgall->Fit("pol0");
  double mean_fit = mgall->GetFunction("pol0")->GetParameter(0);


  c2->cd();
  double max_val = mgall->GetYaxis()->GetXmax();
  double min_val = mgall->GetYaxis()->GetXmin();
  TH1D* hpull = new TH1D("hpull","Pull Fit", 150,-6,6);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  mini_tree->Draw(Form( "(spin*lagr_asym_us_avg*1e9 - %f)/(lagr_asym_us_avg.err*1e9)>>hpull",mean_fit),
		  "arm_flag==0 && kGood==1","goff");
  mini_tree->Draw(Form( "(spin*lagr_asym_usr*1e9 - %f)/(lagr_asym_usr.err*1e9) >>+hpull",mean_fit),
		  "arm_flag==1 && kGood==1","goff");
  mini_tree->Draw(Form( "(spin*lagr_asym_usl*1e9 - %f)/(lagr_asym_usl.err*1e9) >>+hpull",mean_fit),
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
  TLatex *entries = new TLatex(2,170,Form("Entries:%d",(int)hpull->GetEntries()));

  TLatex *chisq = new TLatex(2,140,Form("#chi^{2}/ndf =%.1f/%d",  
					hpull->GetFunction("gaus")->GetChisquare(),
					hpull->GetFunction("gaus")->GetNDF()));
  TLatex *prob = new TLatex(2,125,Form("Prob. =%.1f",  
				       hpull->GetFunction("gaus")->GetProb()));

  TLatex *sigma = new TLatex(2,155,Form("#sigma=%.2f#pm %.2f",  
					hpull->GetFunction("gaus")->GetParameter(2),
					hpull->GetFunction("gaus")->GetParError(2)));
  entries->SetTextFont(22);
  entries->SetTextSize(0.04);
  entries->Draw("same");
  sigma->SetTextFont(22);
  sigma->SetTextSize(0.04);
  sigma->Draw("same");
  chisq->SetTextFont(22);
  chisq->SetTextSize(0.04);
  chisq->Draw("same");
  // prob->SetTextFont(22);
  // prob->SetTextSize(0.04);
  // prob->Draw("same");

  c1->SaveAs("prex_lagrange_minirun_pullfit.pdf");
  c2->SaveAs("prex_lagrange_minirun_pull_distribution.pdf");
}

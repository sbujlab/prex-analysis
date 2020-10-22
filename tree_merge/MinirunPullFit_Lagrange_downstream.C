void MinirunPullFit_Lagrange_downstream(){
  TFile *input = TFile::Open("./rootfiles/MergedLagrange_all_slugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini_lagrall");
  mini_tree->AddFriend("mini_info","./rootfiles/runinfo_all_slugs.root");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,500);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  mini_tree->Draw("spin*lagr_asym_ds_avg*1e9: lagr_asym_ds_avg.err*1e9:run+mini/10.0",
		  "arm_flag==0 && kGood==1","goff");
  TGraphErrors* mger_both = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_both->SetMarkerStyle(7);
  mini_tree->Draw("spin*lagr_asym_dsr*1e9: lagr_asym_dsr.err*1e9:run+mini/10.0",
		  "arm_flag==1 && kGood==1","goff");
  TGraphErrors* mger_dsr = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_dsr->SetMarkerStyle(7);

  mini_tree->Draw("spin*lagr_asym_dsl*1e9: lagr_asym_dsl.err*1e9:run+mini/10.0",
		  "arm_flag==2 && kGood==1","goff");
  TGraphErrors* mger_dsl = new TGraphErrors(mini_tree->GetSelectedRows(),
					    mini_tree->GetV3(),mini_tree->GetV1(),
					    0,mini_tree->GetV2());
  mger_dsl->SetMarkerStyle(7);

  TMultiGraph *mgall = new TMultiGraph();
  mgall->Add(mger_both,"p");
  mgall->Add(mger_dsr,"p");
  mgall->Add(mger_dsl,"p");
  
  mgall->SetTitle("PREX Lagrange Corrected Asymmetry(ppb) - Downstream Main; run; (ppb)"); 
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
  mini_tree->Draw(Form( "(spin*lagr_asym_ds_avg*1e9 - %f)/(lagr_asym_ds_avg.err*1e9)>>hpull",mean_fit),
		  "arm_flag==0 && kGood==1","goff");
  mini_tree->Draw(Form( "(spin*lagr_asym_dsr*1e9 - %f)/(lagr_asym_dsr.err*1e9) >>+hpull",mean_fit),
		  "arm_flag==1 && kGood==1","goff");
  mini_tree->Draw(Form( "(spin*lagr_asym_dsl*1e9 - %f)/(lagr_asym_dsl.err*1e9) >>+hpull",mean_fit),
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

  c1->SaveAs("prex_lagrange_downstream_minirun_pullfit.pdf");
  c2->SaveAs("prex_lagrange_downstream_minirun_pull_distribution.pdf");
}

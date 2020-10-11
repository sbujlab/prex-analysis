void MinirunPullFit_RawReweighted(){
  TFile *input = TFile::Open("./rootfiles/MergedLagrange_all_slugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini");
  mini_tree->AddFriend("mini_lagrall");
  mini_tree->AddFriend("mini_regall");
  mini_tree->AddFriend("mini_info","./rootfiles/runinfo_all_slugs.root");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,500);
  c1->cd();
  // TPad *pad1 = new TPad("pad1","pad1",0,0, 0.66, 1.0);
  // TPad *pad2 = new TPad("pad2","pad2",0.66,0, 1.0, 1.0);
  // pad1->Draw();
  // pad2->Draw();
  // pad1->cd();
  mini_tree->Draw("spin*asym_us_avg*1e9: lagr_asym_us_avg.err*1e9:run+mini/10.0",
		  "arm_flag==0 && kGood==1","goff");
  TGraphErrors* mger_both = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_both->SetMarkerStyle(7);
  mini_tree->Draw("spin*asym_usr*1e9: lagr_asym_usr.err*1e9:run+mini/10.0",
		  "arm_flag==1 && kGood==1","goff");
  TGraphErrors* mger_usr = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_usr->SetMarkerStyle(7);

  mini_tree->Draw("spin*asym_usl*1e9: lagr_asym_usl.err*1e9:run+mini/10.0",
		  "arm_flag==2 && kGood==1","goff");
  TGraphErrors* mger_usl = new TGraphErrors(mini_tree->GetSelectedRows(),
					    mini_tree->GetV3(),mini_tree->GetV1(),
					    0,mini_tree->GetV2());
  mger_usl->SetMarkerStyle(7);

  TMultiGraph *mgall = new TMultiGraph();
  mgall->Add(mger_both,"p");
  mgall->Add(mger_usr,"p");
  mgall->Add(mger_usl,"p");
  
  mgall->SetTitle("PREX Lagrange Raw Asymmetry(ppb) ; run; (ppb)"); 
  gStyle->SetOptFit(1);
  mgall->Draw("AP");
  mgall->Fit("pol0");
  double mean_fit = mgall->GetFunction("pol0")->GetParameter(0);



  // pad2->cd();
  // double max_val = mgall->GetYaxis()->GetXmax();
  // double min_val = mgall->GetYaxis()->GetXmin();
  // TH1D* hpull = new TH1D("hpull","Pull Fit", 200,-5,10);
  // gStyle->SetOptStat("ourMeN");
  // gStyle->SetStatH(0.2);
  // gStyle->SetStatW(0.25);
  // mini_tree->Draw(Form( "(spin*lagr_asym_us_avg*1e9 - %f)/(lagr_asym_us_avg.err*1e9)>>hpull",mean_fit),
  // 		  "arm_flag==0 && kGood==1","goff");
  // mini_tree->Draw(Form( "(spin*lagr_asym_usr*1e9 - %f)/(lagr_asym_usr.err*1e9) >>+hpull",mean_fit),
  // 		  "arm_flag==1 && kGood==1","goff");
  // mini_tree->Draw(Form( "(spin*lagr_asym_usl*1e9 - %f)/(lagr_asym_usl.err*1e9) >>+hpull",mean_fit),
  // 		  "arm_flag==2 && kGood==1","goff");
  // hpull->Draw();
  // hpull->Fit("gaus");
  c1->SaveAs("prex_lagrange_weighted_raw_minirun_pullfit.pdf");
}

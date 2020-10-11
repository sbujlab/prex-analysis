void DoubleCheckReg5bpm(){
  TFile *input = TFile::Open("./rootfiles/PostpanMerged_all_slugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini");
  mini_tree->AddFriend("mini_info","./rootfiles/runinfo_all_slugs.root");
  
  mini_tree->Draw("spin*reg_asym_us_avg*1e9: reg_asym_us_avg.err*1e9:run+mini/10.0",
		  "arm_flag==0 && kGood>0","goff");
  TGraphErrors* mger_both = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_both->SetMarkerStyle(20);
  mini_tree->Draw("spin*reg_asym_usr*1e9: reg_asym_usr.err*1e9:run+mini/10.0",
		  "arm_flag==1 && kGood>0","goff");
  TGraphErrors* mger_usr = new TGraphErrors(mini_tree->GetSelectedRows(),
					     mini_tree->GetV3(),mini_tree->GetV1(),
					     0,mini_tree->GetV2());
  mger_usr->SetMarkerStyle(20);

  mini_tree->Draw("spin*reg_asym_usl*1e9: reg_asym_usl.err*1e9:run+mini/10.0",
		  "arm_flag==2 && kGood>0","goff");
  TGraphErrors* mger_usl = new TGraphErrors(mini_tree->GetSelectedRows(),
					    mini_tree->GetV3(),mini_tree->GetV1(),
					    0,mini_tree->GetV2());
  mger_usl->SetMarkerStyle(20);

  TMultiGraph *mgall = new TMultiGraph();
  mgall->Add(mger_both,"p");
  mgall->Add(mger_usr,"p");
  mgall->Add(mger_usl,"p");
  
  mgall->SetTitle("PREX Postpan Asymmetry(ppb) ; run; (ppb)"); 
  gStyle->SetOptFit(1);
  mgall->Draw("AP");
  mgall->Fit("pol0");
}

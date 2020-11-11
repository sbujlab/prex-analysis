void MinirunPullFit_deltaA(){
  TFile *input = TFile::Open("./rootfiles/MergedLagrange_all_slugs.root");
  TTree *mini_tree = (TTree*)input->Get("mini_lagrall");
  mini_tree->AddFriend("mini_regall");
  mini_tree->AddFriend("mini_info","./rootfiles/runinfo_all_slugs.root");
  
  // TCanvas *c1 = new TCanvas("c1","c1",1200,500);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  mini_tree->Draw("spin*(lagr_asym_us_avg-reg_asym_us_avg)*1e9",
		  "arm_flag==0 && kGood==1 && slug>86");
  // c1->cd();
  // mini_tree->Draw("spin*(lagr_asym_us_avg-reg_asym_us_avg)*1e9:run",
  // 		  "arm_flag==0 && kGood==1 && slug==32","*");
  
}

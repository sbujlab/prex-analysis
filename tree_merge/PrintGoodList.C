void PrintGoodList(){
  TFile *input = TFile::Open("rootfiles/runinfo_all_slugs.root");
  TTree *mini_info = (TTree*)input->Get("mini_info");
  
  TFile *postpan = TFile::Open("rootfiles/PostpanMerged_all_slugs.root");
  TTree *mini_tree = (TTree*)postpan->Get("mini");
  
  cout << " All good run " <<  mini_info->GetEntries("kGood>0") << endl;
  cout << " All good run with filter " <<  mini_info->GetEntries("kGood==1") << endl;
  Int_t npt = mini_info->GetEntries();
  Int_t run, minirun, kGood, spin,arm_flag;
  mini_info->SetBranchAddress("run",&run);
  mini_info->SetBranchAddress("mini",&minirun);
  mini_info->SetBranchAddress("kGood", &kGood);
  mini_info->SetBranchAddress("spin", &spin);
  mini_info->SetBranchAddress("arm_flag",&arm_flag);

  Double_t reg_asym_usl,reg_asym_usr,reg_asym_us_avg;
  mini_tree->GetBranch("reg_asym_usl")->GetLeaf("mean")->SetAddress(&reg_asym_usl);
  mini_tree->GetBranch("reg_asym_usr")->GetLeaf("mean")->SetAddress(&reg_asym_usr);
  mini_tree->GetBranch("reg_asym_us_avg")->GetLeaf("mean")->SetAddress(&reg_asym_us_avg);

  Double_t reg_asym_usl_err,reg_asym_usr_err,reg_asym_us_avg_err;
  mini_tree->GetBranch("reg_asym_usl")->GetLeaf("err")->SetAddress(&reg_asym_usl_err);
  mini_tree->GetBranch("reg_asym_usr")->GetLeaf("err")->SetAddress(&reg_asym_usr_err);
  mini_tree->GetBranch("reg_asym_us_avg")->GetLeaf("err")->SetAddress(&reg_asym_us_avg_err);
  
  FILE *output = fopen("prex_mini_avg_signed_tao.txt","w");
  for(int i=0;i<npt;i++){
    mini_info->GetEntry(i);
    mini_tree->GetEntry(i);
    if(kGood>0){
      fprintf(output,"%d,%d", run,minirun);
      if(arm_flag==0)
      	fprintf(output,",%.2f,%.2f\n",
		spin*reg_asym_us_avg*1e9,
		reg_asym_us_avg_err*1e9);
      if(arm_flag==1)
      	fprintf(output,",%.2f,%.2f\n",
		spin*reg_asym_usr*1e9,
		reg_asym_usr_err*1e9);
      if(arm_flag==2)
      	fprintf(output,",%.2f,%.2f\n",
		spin*reg_asym_usl*1e9,
		reg_asym_usl_err*1e9);
    }
  }
  fclose(output);
  input->Close();
}

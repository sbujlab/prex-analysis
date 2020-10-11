void PrintGoodList(){
  TFile *input = TFile::Open("allslug_runinfo.root");
  TTree *mini_info = (TTree*)input->Get("mini_info");
  cout << " All good run " <<  mini_info->GetEntries("kGood>0") << endl;
  cout << " All good run with filter " <<  mini_info->GetEntries("kGood==1") << endl;
  Int_t npt = mini_info->GetEntries();
  Int_t run, minirun, kGood;
  mini_info->SetBranchAddress("run",&run);
  mini_info->SetBranchAddress("mini",&minirun);
  mini_info->SetBranchAddress("kGood", &kGood);
  FILE *output = fopen("prex_runlist_tao.txt","w");
  for(int i=0;i<npt;i++){
    mini_info->GetEntry(i);
    if(kGood>0){
      fprintf(output,"%d, %d\n", run,minirun);
    }
  }
  fclose(output);
  input->Close();
}

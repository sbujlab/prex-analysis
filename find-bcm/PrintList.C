void PrintList(){
  TFile *input = TFile::Open("prex_bcm_status.root");
  TTree *T = (TTree*) input->Get("T");

  Int_t nEntries = T->GetEntries();
  TString device_list[]={"bcm_an_us","bcm_an_ds",
			 "bcm_an_ds3","bcm_an_ds10",
			 "bcm_dg_us","bcm_dg_ds",
			 "unser"};
  Int_t slug,run_number,bcm_index;
  Int_t bcm_flag[7];
  T->SetBranchAddress("run",&run_number);
  T->SetBranchAddress("slug",&slug);
  T->SetBranchAddress("bcm_index",&bcm_index);
  FILE *output = fopen("normalizing_bcm.txt","w");
  for(int ievt=0;ievt<nEntries;ievt++){
    T->GetEntry(ievt);
    fprintf(output,"%d,%d,%s\n",
	    run_number,
	    bcm_index,
	    device_list[bcm_index].Data());
  }
  fclose(output);
}

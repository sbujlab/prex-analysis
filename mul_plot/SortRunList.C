void SortRunList(){
  TString list_array[]={"50uA",
			"50uA_1arm",
			"70uA_120Hz",
			"70uA_120Hz_1arm",
			"70uA_240Hz",
			"70uA_240Hz_1arm",
			"85uA"};

  TString fCuts[]={"yield_bcm_target<53&&kGood==1 && arm_flag==0",
		   "yield_bcm_target<53&&kGood==1 && arm_flag!=0",
		   "yield_bcm_target>55&&yield_bcm_target<73&& kGood==1 && arm_flag==0 && run<3877",
		   "yield_bcm_target>55&&yield_bcm_target<73&& kGood==1 && arm_flag!=0 && run<3877",
		   "yield_bcm_target>55&&yield_bcm_target<75&& kGood==1 && arm_flag==0 && run>=3877",
		   "yield_bcm_target>55&&yield_bcm_target<75&& kGood==1 && arm_flag!=0 && run>=3877",
		   "yield_bcm_target>75&& kGood==1 && arm_flag==0 && run>3877"};
  			
  Int_t nSel=  sizeof(list_array)/sizeof(*list_array);
  
  TFile *merged_file = TFile::Open("../tree_merge/rootfiles/MergedSum_all_slugs.root");
  TTree *burst_tree = (TTree*)merged_file->Get("burst");
  burst_tree->AddFriend("mini_info","../tree_merge/rootfiles/runinfo_all_slugs.root");
  
  for(int i=0;i<nSel;i++){
    FILE *output = fopen(("./list/"+list_array[i]).Data(),"w");
    Int_t npt = burst_tree->Draw("run",fCuts[i],"goff");
    double *fRun = burst_tree->GetV1();
    double run_number = -1;

    for(int ipt=0;ipt<npt;ipt++){
      if(run_number==fRun[ipt])
	continue;
      
      run_number = fRun[ipt];
      fprintf(output,"%d\n",(int)run_number);
    }
    fclose(output);
  }

}

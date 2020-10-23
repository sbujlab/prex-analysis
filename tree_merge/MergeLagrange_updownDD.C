#include "TaSumStat.cc"
void MergeLagrange_updownDD(Int_t slug);
void MergeLagrange_updownDD(Int_t slug){
  TString label = Form("slug%d",slug);
  TString filename  = label+".list";

  if(slug==9999){
    filename = "slug1-94.list";
    label = "all_slugs";
  }

  TString path = "./prex-runlist/simple_list/";
  FILE *runlist = fopen((path+filename).Data(),"r");

  vector<TString> tree_name_array={"mini","mini_lagrall","mini_regall"};
  vector<TTree*> output_tree_array;
  Int_t nTree = tree_name_array.size();
  vector<Int_t> fRunNumberArray;
  TString lagrange_path= "/lustre/expphy/volatile/halla/parity/prex-respin2/LagrangeOutput/rootfiles/";
  TaSumStat fSumStatBuilder;
  if(runlist!=NULL){
    TFile *output = TFile::Open("./rootfiles/MergedLagrange_updownDD_"+label+".root","RECREATE");
    while(!feof(runlist)){
      Int_t run_number=0;
      fscanf(runlist,"%d/n",&run_number);
      if(run_number==0)
	continue;
      fRunNumberArray.push_back(run_number);
      Int_t seg_number=0;
      TFile *this_file;
      TString rootfile_name = lagrange_path+Form("prexRespin2_lagrange_eigen_updownDD_%d.%03d.root",
						run_number,seg_number);
      if(gSystem->AccessPathName(rootfile_name)==0){
	this_file = TFile::Open(rootfile_name);
	cout << " -- Found ROOT file ";
	cout << this_file->GetName() << endl;
	vector<TString>::iterator it_tree  = tree_name_array.begin();

	while(it_tree!=tree_name_array.end()){
	  TTree *aTree_ptr = (TTree*)this_file->Get(*it_tree);
	  if(aTree_ptr==NULL){
	    cout << " ** TTree " 
		 << *it_tree
		 << " is not found. " << endl;
	    continue;
	  }
	  fSumStatBuilder.collect_branchlist_from_postpan(aTree_ptr);
	  it_tree++;
	} //end of tree loop
	// seg_number++;
	// rootfile_name = postpan_path+Form("prexPrompt_%d_%03d_regress_postpan.root",
	// 				  run_number,seg_number);
	this_file->Close();
      } // end of if file exists
    } // end of the loop over run in a list 
    fclose(runlist);

    fSumStatBuilder.construct_branches(output);

    // ===== Second-pass
    auto iter_run  = fRunNumberArray.begin();
    while(iter_run!=fRunNumberArray.end()){
      Int_t run_number = *iter_run;
      Int_t seg_number=0;
      TString rootfile_name = lagrange_path+Form("prexRespin2_lagrange_eigen_%d.%03d.root",
						run_number,seg_number);
      TFile *this_file;
      if(gSystem->AccessPathName(rootfile_name)==0){
    	this_file = TFile::Open(rootfile_name);
    	cout << " -- Reading ROOT file ";
    	cout << this_file->GetName() << endl;
    	for(int it=0;it<nTree;it++){
    	  TTree *aTree = (TTree*)this_file->Get(tree_name_array[it]);
    	  if(aTree==NULL){  // just in case a tree output is missing
    	    fSumStatBuilder.load_null_postpan_stat_by_name(tree_name_array[it]);
    	    cout << " ** TTree " 
    		 << tree_name_array[it] 
    		 << " is not found and stat is filled with zeros " << endl;
    	    continue;
    	  }else
    	    fSumStatBuilder.load_postpan_stat_ptr(aTree);
    	  fSumStatBuilder.set_run_number(run_number);
    	  Int_t burst_counter=0;
    	  Int_t nEntries = aTree->GetEntries();
    	  for(int ievt=0;ievt<nEntries;ievt++){
    	    aTree->GetEntry(ievt);
    	    burst_counter = ievt;
    	    fSumStatBuilder.set_burst_counter(burst_counter);
    	    fSumStatBuilder.set_minirun_counter(burst_counter);
    	    fSumStatBuilder.write_sum_postpan_stat_by_name(tree_name_array[it]);
    	    fSumStatBuilder.fill_tree_by_name( tree_name_array[it]);
    	  }
    	} // end of tree loop;
    	this_file->Close();
      } // end of file exists
      iter_run++;
    } // end of run number loop
    output->cd();
    cout << " -- Writing Tree " << endl;
    fSumStatBuilder.write_trees_to_output(output);
    cout << "-- Closing " << output->GetName()<< endl;
    output->Close();
    cout << "-- Saved output " << output->GetName()<< endl;
  } // end of if runlist exists
}


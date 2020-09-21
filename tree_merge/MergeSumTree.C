#include "TaSumStat.cc"
void MergeSumTree(Int_t slug);

void MergeSumTree(){
  for(int i=1;i<=94;i++)
    MergeSumTree(i);
}

void MergeSumTree(Int_t slug){
  TString label = Form("slug%d",slug);
  TString filename  = label+".list";
  TString path = "./prex-runlist/simple_list/";
  FILE *runlist ;

  if(slug==9999){
    filename = "all.list";
    label = "all_slugs";
  }
  runlist = fopen((path+filename).Data(),"r");
  
  if(runlist==NULL){
    cout << " -- Error: runlist " 
	 << filename
	 << " is not found! " << endl;
    return;
  }

  TFile *output = TFile::Open("./rootfiles/MergedSum_"+label+".root","RECREATE");
  TaSumStat fSumStatBuilder;
  TString fTreeNameArray[] = {"muls","burst","burst_mulc",
			      "burst_mulc_dit","burst_mulc_dit_combo",
			      "burst_mulc_lrb_burst","burst_mulc_lrb_alldet_burst",
			      "burst_mulc_lrb_all_burst",
			      "burst_mulc_lrb","burst_mulc_lrb_alldet",
			      "burst_mulc_lrb_all"};
  Int_t nTree = sizeof(fTreeNameArray)/sizeof(*fTreeNameArray);
  vector<Int_t> fRunNumberArray;
  // ===== First-Pass
  while(!feof(runlist)){
    Int_t run_number=0;
    Int_t seg_number=0;
    fscanf(runlist,"%d/n",&run_number);
    if(run_number==0)
      continue;
    fRunNumberArray.push_back(run_number);
    TString qw_path = TString(gSystem->Getenv("QW_ROOTFILES") );
    TString rootfile_name = Form("prexPrompt_pass2_%d.%03d.root",
				 run_number,seg_number);
    TFile *this_file;

    if(gSystem->AccessPathName(qw_path+rootfile_name)==0){
      this_file = TFile::Open(qw_path+rootfile_name);
      cout << " -- Found ROOT file ";
      cout << this_file->GetName() << endl;
      for(int it=0;it<nTree;it++){
	TTree *aTree = (TTree*)this_file->Get(fTreeNameArray[it]);
	if(aTree==NULL){
	  cout << " ** TTree " 
	       << fTreeNameArray[it] 
	       << " is not found. " << endl;
	  continue;
	}
	fSumStatBuilder.collect_branchlist_from_input(aTree);
      } // end of tree loop;

      this_file->Close();

    } // end of if file exists
  } // end of runlist loop

  fSumStatBuilder.construct_branches(output);
  // ===== Second-pass
  auto iter_run  = fRunNumberArray.begin();
  while(iter_run!=fRunNumberArray.end()){
    Int_t run_number = *iter_run;
    Int_t seg_number=0;
    TString qw_path = TString(gSystem->Getenv("QW_ROOTFILES") );
    TString rootfile_name = Form("prexPrompt_pass2_%d.%03d.root",
				 run_number,seg_number);
    TFile *this_file;
    if(gSystem->AccessPathName(qw_path+rootfile_name)==0){
      this_file = TFile::Open(qw_path+rootfile_name);
      cout << " -- Reading ROOT file ";
      cout << this_file->GetName() << endl;
      for(int it=0;it<nTree;it++){
	TTree *aTree = (TTree*)this_file->Get(fTreeNameArray[it]);
	if(aTree==NULL){  // just in case a tree output is missing
	  fSumStatBuilder.load_null_stat_by_name(fTreeNameArray[it]);
	  cout << " ** TTree " 
	       << fTreeNameArray[it] 
	       << " is not found and stat is filled with zeros " << endl;
	  continue;
	}else
	  fSumStatBuilder.load_japan_stat_ptr(aTree);
	fSumStatBuilder.set_run_number(run_number);
	Int_t burst_counter=0;
	Int_t nEntries = aTree->GetEntries();

	for(int ievt=0;ievt<nEntries;ievt++){
	  aTree->GetEntry(ievt);
	  burst_counter = ievt;
	  fSumStatBuilder.set_burst_counter(burst_counter);
	  fSumStatBuilder.set_minirun_counter(burst_counter);
	  fSumStatBuilder.write_sum_stat_by_name(fTreeNameArray[it]);
	  fSumStatBuilder.fill_tree_by_name( fTreeNameArray[it]);
	}
      } // end of tree loop;
      this_file->Close();
    } // end of file exists
    iter_run++;
  } // end of run number loop

  output->cd();
  cout << " -- Writing Tree " << endl;
  fSumStatBuilder.write_trees_to_output(output);
  cout << " -- Writing Tree : Done " << endl;
  cout << " -- Closing " << output->GetName() << endl;
  output->Close();
}


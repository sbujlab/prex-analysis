/*
author: Cameron Clarke
last update: February 2021

*/

void Add_RCDB_from_Agg_All(TString agginfile, TString eiginfile, TString eigoutfile, TString treename = "mini"){ //(int start, int end)
  // FIXME replace with read-from-agg tree
  TFile agg_input(agginfile);

  TFile eig_input(eiginfile);
  TFile output(eigoutfile,"recreate");
  //TFile eig_input(Form("dataRootfiles/Full_mini_eigen_reg_allbpms_AT.root"));
  //TFile output(Form("rootfiles/rcdb_AT_eigenvectors.root"),"recreate");

  TTree *agg_tree;
  TTree *mini_tree;
  agg_input.GetObject("agg", agg_tree);
  eig_input.GetObject(treename, mini_tree);

  TChain* out_mini_tree = (TChain*)mini_tree -> CloneTree(0);

  Int_t run, burst_counter;
  mini_tree->SetBranchAddress("run",&run);
  mini_tree->SetBranchAddress("mini",&burst_counter);

  // Add new branches
  // 
  // - Loop over all "rcdb_" branches in agg tree
  // - - agg->SetBranchAddress(name,placeholder_agg);
  // - - make clone/copies of target trees
  // - - eig_reg_tree_part_avg->Branch(name,placeholder_new);
  // - Build index on run, minirun for all trees
  // - Loop over all entries in mini/eig/reg/lagr
  // - - Find corresponding indexed entry in agg
  // - - Set new = old values
  // - Write tree
  //
  // Adding new branches?
  //eig_reg_tree_part_avg->Branch();

  std::string index1 = "run_number";
  std::string index2 = "minirun_n";
  agg_tree->BuildIndex(index1.c_str(),index2.c_str());
  mini_tree->BuildIndex("run","mini");
  TTreeIndex *mini_indexedTree = (TTreeIndex *) mini_tree->GetTreeIndex();
  TTreeIndex *agg_indexedTree = (TTreeIndex *) agg_tree->GetTreeIndex();

  TObjArray *var_list=agg_tree->GetListOfBranches();
  TIter var_iter(var_list);
  std::vector<Double_t> old_vals;
  std::vector<Double_t> new_vals;
  std::vector<TString> names;
  Int_t nRCDB = 0;
  
  while (auto *var= var_iter.Next()){
    TString name(var->GetName());
    if (name.Contains("rcdb_") || name.Contains("_main_det_")) {
      names.push_back(name);
      //double old_tmpval = 0.0;
      ////double new_tmpval = 0.0;
      old_vals.push_back(0.0);
      new_vals.push_back(0.0);
      nRCDB++;
    }
  }

  for (int m = 0 ; m < nRCDB ; m++) {
    // Add new branch to out_mini_tree
    agg_tree->SetBranchAddress(names[m],&old_vals[m]);
    out_mini_tree->Branch(names[m],&new_vals[m]);
  }

  Long64_t nEntries = mini_tree->GetEntries();
  Int_t last_run = -1;

  for(int ievt=0;ievt<nEntries;ievt++){
    mini_tree           -> GetEntry(ievt);
    agg_tree->GetEntryWithIndex((Double_t)run,(Double_t)burst_counter);
    // For each branch in rcdb branches list set new = old
    //////// lagr_tree->GetEntry(ievt);
    for (int o = 0 ; o < nRCDB ; o++) {
      new_vals.at(o) = old_vals.at(o);
      //Printf("RCDB data %d = %f",o,old_vals.at(o));
    }
    out_mini_tree           -> Fill();
  }

  out_mini_tree           -> Write(treename,TObject::kOverwrite);

  TKey *key;
  TIter nextkey(eig_input.GetListOfKeys(),kIterBackward);
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname); 
    if (!cl) continue;
    if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)eig_input.Get(key->GetName());
      if ((TString)key->GetName() == treename || ((TString)key->GetName()).Contains("sum")) {
        // Skip the tree we've been working with
        // Also ignore the run-wise summary trees
        continue;
      }
      // Avoid writing the data of a TTree more than once.
      // Note this assume that older cycles are (as expected) older
      // snapshots of the TTree meta data.
      if (!output.FindObject(key->GetName())) {
        output.cd();
        TTree *newT = T->CloneTree(-1,"fast");
        newT->Write();
      }
    }
  }

  eig_input.Close();
  agg_input.Close();
  output.Close();
}

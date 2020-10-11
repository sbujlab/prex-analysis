#include "TaSumStat.hh"
TaSumStat::TaSumStat(){

  aUnit.ppm = 1e-6;
  aUnit.ppb = 1e-9;
  aUnit.mm = 1.0;
  aUnit.um = 1e-3;
  aUnit.nm = 1e-6;

}
JAPAN_STAT TaSumStat::init_japan_stat(){
  JAPAN_STAT this_stat;
  this_stat.hw_sum=0.0;
  this_stat.hw_sum_m2=0.0;
  this_stat.hw_sum_err=0.0;
  this_stat.num_samples=0.0;
  return this_stat;
}

POSTPAN_STAT TaSumStat::init_postpan_stat(){
  POSTPAN_STAT this_stat;
  this_stat.mean=0.0;
  this_stat.err=0.0;
  this_stat.rms=0.0;
  this_stat.num_samples=0.0;
  return this_stat;
}

JAPAN_STAT TaSumStat::invalid_japan_stat(){
  JAPAN_STAT this_stat;
  this_stat.hw_sum=0.0;
  this_stat.hw_sum_m2=-1;
  this_stat.hw_sum_err=-1;
  this_stat.num_samples=0.0;
  return this_stat;
}
void TaSumStat::merge_japan_stat(JAPAN_STAT &dest_stat,JAPAN_STAT in_stat){
  // for merging short minirun to previous one

  if(in_stat.num_samples==0){ // to avoid div-0 problem
    cout << " **  "
	 << __PRETTY_FUNCTION__
	 << " number of samples is zero, will be skipped " << endl;
    return;
  }
  double mean_1 = dest_stat.hw_sum;
  double m2_1  = dest_stat.hw_sum_m2;
  double nsamp_1 = dest_stat.num_samples;

  double mean_2 = in_stat.hw_sum;
  double m2_2  = in_stat.hw_sum_m2;
  double nsamp_2 = in_stat.num_samples;
  double delta_mean = mean_2 - mean_1;

  m2_1 += m2_2;
  m2_1 += nsamp_1*nsamp_2*pow(delta_mean,2)/(nsamp_1+nsamp_2);

  dest_stat.hw_sum += nsamp_2*delta_mean/(nsamp_1+nsamp_2);
  dest_stat.num_samples+=nsamp_2;
  dest_stat.hw_sum_m2 = m2_1;
  dest_stat.hw_sum_err = TMath::Sqrt(dest_stat.hw_sum_m2)/dest_stat.num_samples;
} 

void TaSumStat::write_sum_stat(SUM_STAT &dest_stat,JAPAN_STAT in_stat){
  dest_stat.mean = in_stat.hw_sum;

  if(in_stat.num_samples==0){// to avoid div-0 problem
    dest_stat.err = -1;
    dest_stat.rms = -1;
  }else{
    dest_stat.err = in_stat.hw_sum_err;
    dest_stat.rms = TMath::Sqrt(in_stat.hw_sum_m2/in_stat.num_samples);
  }
  dest_stat.num_samples  = in_stat.num_samples;    
}


void TaSumStat::write_sum_stat(SUM_STAT &dest_stat,POSTPAN_STAT in_stat){
  dest_stat.mean = in_stat.mean;

  if(pow(in_stat.rms/in_stat.err,2)<=0 || in_stat.err<0){// to avoid div-0 problem
    dest_stat.err = -1;
    dest_stat.rms = -1;
    dest_stat.num_samples = 0.0;
  }else{
    dest_stat.err = in_stat.err;
    dest_stat.rms = in_stat.rms;
    dest_stat.num_samples = pow(in_stat.rms/in_stat.err,2);
  }
}


void TaSumStat::write_sum_stat_by_name(TString tree_name){
  Int_t nch = fJStatMap[tree_name].size();
    for( int ich=0;ich<nch;ich++){
      write_sum_stat( fSumStatMap[tree_name][ich],
		      fJStatMap[tree_name][ich]);
    } // end of branch loop;
}
void TaSumStat::write_sum_postpan_stat_by_name(TString tree_name){
  Int_t nch = fJStatMap[tree_name].size();
    for( int ich=0;ich<nch;ich++){
      write_sum_stat( fSumStatMap[tree_name][ich],
		      fPStatMap[tree_name][ich]);
    } // end of branch loop;
}

void TaSumStat::fill_tree_by_name(TString tree_name){
  fTreeMap[tree_name]->Fill();
}

void TaSumStat::construct_branches(TFile* output){
  output->cd();
  TString leaflist = "mean/D:err:rms:num_samples";
  TTree* aTree;
  Int_t nTree = fTreeNameList.size();
  for(int it=0;it<nTree;it++){
    TString fTreeName = fTreeNameList[it];
    fTreeMap[fTreeName] =  new TTree(fTreeName,"");
    fTreeArray.push_back( fTreeMap[fTreeName] );
  }

  for(int it=0;it<nTree;it++){
    aTree = fTreeArray[it];
    aTree->Branch("unit",&aUnit,"ppm/D:ppb:um:mm:nm");
    aTree->Branch("run",&run_number,"run/I");
    aTree->Branch("BurstCounter",&burst_counter,"BurstCounter/I");
    aTree->Branch("minirun",&minirun_counter,"minirun/I");
    TString fTreeName = aTree->GetName();
    vector<TString> channel_list = fBranchNameListMap[fTreeName];
    Int_t nch  = channel_list.size();
    JAPAN_STAT fzero = init_japan_stat(); // seems necessary;
    POSTPAN_STAT fzero_pp = init_postpan_stat(); // seems necessary;
    vector<SUM_STAT> fSumStat(nch);
    vector<JAPAN_STAT> fJStat(nch,fzero);
    vector<POSTPAN_STAT> fPStat(nch,fzero_pp);
    fSumStatMap[fTreeName] = fSumStat;
    fJStatMap[fTreeName] = fJStat;
    fPStatMap[fTreeName] = fPStat;
    for(int ich=0;ich<nch;ich++){
      aTree->Branch(channel_list[ich],&(fSumStatMap[fTreeName][ich]),leaflist);
    }// end of branch loop
  } // end of tree loop
} 

void TaSumStat::collect_branchlist_from_input(TTree* input_tree){
  TString fTreeName = input_tree->GetName();
  Bool_t isNewTree = kFALSE;
  if(find(fTreeNameList.begin(),fTreeNameList.end(),fTreeName)==fTreeNameList.end()){
    fTreeNameList.push_back(fTreeName);
    isNewTree = kTRUE;
  }
  auto fBranchList = input_tree->GetListOfBranches();
  Int_t nbr = fBranchList->GetEntries();
  for(int ibr =0; ibr<nbr;ibr++){
    TBranch* fbranch = dynamic_cast<TBranch*> (fBranchList->At(ibr));
    if(fbranch->GetLeaf("hw_sum_m2")==NULL)
      continue;
    TString branch_name = fbranch->GetName();
    vector<TString> current_list = fBranchNameListMap[fTreeName];
    if(find(current_list.begin(),current_list.end(),branch_name)==current_list.end()){
      fBranchNameListMap[fTreeName].push_back(branch_name);
      if(!isNewTree)
	cout << " ++ Append new TBranch "
	     << branch_name 
	     << " to TTree " << fTreeName << endl;
    }
  }
}

void TaSumStat::collect_branchlist_from_postpan(TTree* input_tree){
  TString fTreeName = input_tree->GetName();
  Bool_t isNewTree = kFALSE;
  if(find(fTreeNameList.begin(),fTreeNameList.end(),fTreeName)==fTreeNameList.end()){
    fTreeNameList.push_back(fTreeName);
    isNewTree = kTRUE;
  }
  auto fBranchList = input_tree->GetListOfBranches();
  Int_t nbr = fBranchList->GetEntries();
  for(int ibr =0; ibr<nbr;ibr++){
    TBranch* fbranch = dynamic_cast<TBranch*> (fBranchList->At(ibr));
    if(fbranch->GetLeaf("mean")==NULL)
      continue;
    TString branch_name = fbranch->GetName();
    vector<TString> current_list = fBranchNameListMap[fTreeName];
    if(find(current_list.begin(),current_list.end(),branch_name)==current_list.end()){
      fBranchNameListMap[fTreeName].push_back(branch_name);
      if(!isNewTree)
	cout << " ++ Append new TBranch "
	     << branch_name 
	     << " to TTree " << fTreeName << endl;
    }
  }
}

void TaSumStat::write_trees_to_output(TFile* output){
  output->cd();
  Int_t nTree = fTreeArray.size();
  for(int it=0;it<nTree;it++)
    fTreeArray[it]->Write();
}

void TaSumStat::load_japan_stat_ptr(TTree* input_tree){
  TString fTreeName = input_tree->GetName();
  vector<TString> fBranchNameList = fBranchNameListMap[fTreeName];
  Int_t nbr = fBranchNameList.size();
  for(int ibr=0;ibr<nbr;ibr++){
    TBranch* branch_ptr = input_tree->GetBranch(fBranchNameList[ibr]);
    if(branch_ptr==NULL){
      fJStatMap[fTreeName][ibr] = init_japan_stat();
    }else{
      branch_ptr->GetLeaf("hw_sum")->SetAddress(&(fJStatMap[fTreeName][ibr].hw_sum));
      branch_ptr->GetLeaf("hw_sum_m2")->SetAddress(&(fJStatMap[fTreeName][ibr].hw_sum_m2));
      branch_ptr->GetLeaf("hw_sum_err")->SetAddress(&(fJStatMap[fTreeName][ibr].hw_sum_err));
      branch_ptr->GetLeaf("num_samples")->SetAddress(&(fJStatMap[fTreeName][ibr].num_samples));
    }
  }
}

void TaSumStat::load_postpan_stat_ptr(TTree* input_tree){
  TString fTreeName = input_tree->GetName();
  vector<TString> fBranchNameList = fBranchNameListMap[fTreeName];
  Int_t nbr = fBranchNameList.size();
  for(int ibr=0;ibr<nbr;ibr++){
    TBranch* branch_ptr = input_tree->GetBranch(fBranchNameList[ibr]);
    if(branch_ptr==NULL){
      fJStatMap[fTreeName][ibr] = init_japan_stat();
    }else{
      branch_ptr->GetLeaf("mean")->SetAddress(&(fPStatMap[fTreeName][ibr].mean));
      branch_ptr->GetLeaf("err")->SetAddress(&(fPStatMap[fTreeName][ibr].err));
      branch_ptr->GetLeaf("rms")->SetAddress(&(fPStatMap[fTreeName][ibr].rms));
    }
  }
}

void TaSumStat::load_null_stat_by_name(TString treename){
  Int_t nch = fJStatMap[treename].size();
  for(int ich=0;ich<nch;ich++)
    fJStatMap[treename][ich] = init_japan_stat();
}

void TaSumStat::load_null_postpan_stat_by_name(TString treename){
  Int_t nch = fPStatMap[treename].size();
  for(int ich=0;ich<nch;ich++)
    fPStatMap[treename][ich] = init_postpan_stat();
}

void TaSumStat::load_invalid_stat_by_name(TString treename){
  Int_t nch = fJStatMap[treename].size();
  for(int ich=0;ich<nch;ich++)
    fJStatMap[treename][ich] = invalid_japan_stat();
}
void TaSumStat::cache_japan_stat(TString treename){
  fJStatBuffer = fJStatMap[treename];  
}

void TaSumStat::merge_japan_stat(TString treename){
  Int_t nch = fJStatBuffer.size();
  for(int ich=0;ich<nch;ich++){
    merge_japan_stat(fJStatMap[treename][ich],fJStatBuffer[ich]);
  }
}

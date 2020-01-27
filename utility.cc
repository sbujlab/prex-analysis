TString get_basename(TString input){
  input.ReplaceAll("asym_","");
  input.ReplaceAll("diff_","");
  return input;
}

vector<Int_t> ParseRunList(TString list_name){
  vector<Int_t> fRet;
  FILE *runlist = fopen(list_name.Data(),"r");
  if(runlist==NULL)
    return fRet;
  
  while(!feof(runlist)){
    Int_t run_number=0;
    fscanf(runlist,"%d\n",&run_number);
    if(run_number!=0){
      fRet.push_back(run_number);
      cout << run_number << endl;
    }
  }
  fclose(runlist);
  return fRet;
}

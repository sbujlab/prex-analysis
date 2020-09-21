map<Int_t,Int_t> LoadArmMapBySlug(Int_t slug_id){
  map<Int_t,Int_t> fret;
  TString info_filename =Form("prex-runlist/slug%d_info.list",slug_id);
  ifstream slug_info;
  slug_info.open(info_filename.Data());
  TString sline;
  while(sline.ReadLine(slug_info)){
    TObjArray *token = sline.Tokenize(',');
    Int_t run_number = (((TObjString*)(token->At(0)))->GetString()).Atoi();
    Int_t arm_flag = (((TObjString*)(token->At(6)))->GetString()).Atoi();
    fret[run_number]=arm_flag;
  }
  slug_info.close();
  return fret;
}

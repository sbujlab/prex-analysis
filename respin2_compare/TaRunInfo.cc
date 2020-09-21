#include "TaRunInfo.hh"
ClassImp(TaRunInfo)
TaRunInfo::TaRunInfo(){
  run_number = -1;
  slug_number = -1;
  arm_flag = -1;
  target_id = -1;
  ihwp = "Unknown";
  wien = "Unknown";
  sign=0.0;
}

Bool_t TaRunInfo::ParseLine(TString sline){

  TObjArray *token = sline.Tokenize(',');

  if(token->GetEntries()!=7){
    cout << " -- Error: Not able to identify run info: " << endl;
    cout << sline << endl;
    return kFALSE;
  }

  run_number = (((TObjString*)(token->At(0)))->GetString()).Atoi();
  slug_number = (((TObjString*)(token->At(1)))->GetString()).Atoi();
  TString target_name = ((TObjString*)(token->At(2)))->GetString();
  target_id = ParseTargetName(target_name);
  run_flag = ((TObjString*)(token->At(3)))->GetString();
  ihwp = ((TObjString*)(token->At(4)))->GetString();
  wien = ((TObjString*)(token->At(5)))->GetString();
  arm_flag = (((TObjString*)(token->At(6)))->GetString()).Atoi();
  
  if( ihwp == "IN" )
    sign = 1;
  else if (ihwp == "OUT" )
    sign =-1;

  if( wien == "FLIP-RIGHT" )
    sign *= 1;
  else if ( wien == "FLIP-LEFT" )
    sign *=-1;

  return kTRUE;
}

Int_t TaRunInfo::ParseTargetName(TString input){
  Ssiz_t begin = input.First('b');
  Ssiz_t end = input.Last('-');
  Ssiz_t length = end-begin;

  TString targetID_str = input(begin+1,length);
  return targetID_str.Atoi();
}


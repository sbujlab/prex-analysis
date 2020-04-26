#include "lib/TaRunInfo_v2.hh"

map<Int_t, TaRunInfo> LoadRunInfoMap(){
  map< Int_t, TaRunInfo > fRunInfoMap;  
  for(int islug=1;islug<=94;islug++){
    TString info_filename = Form("./prex-runlist/slug%d_info.list",islug);
    ifstream slug_info;
    slug_info.open(info_filename.Data());

    if(!slug_info.is_open()){
      cout << " slug info " << info_filename << " does not exist." << endl;
      continue;
    }
    
    TString sline;
    Int_t slug_id=-1;
    Int_t arm_flag=-1;
    TString ihwp_stat;
    TString wien_mode;

    while(sline.ReadLine(slug_info)){
      TaRunInfo fRunInfo;
      if(!fRunInfo.ParseLine(sline))
	continue;
      if(fRunInfo.GetRunFlag()=="Bad")
      	continue;
      
      Int_t run_number = fRunInfo.GetRunNumber();
      fRunInfoMap[run_number]=fRunInfo;
    }
    slug_info.close();
  }
  cout << " -- done " << endl;
  return fRunInfoMap;
}

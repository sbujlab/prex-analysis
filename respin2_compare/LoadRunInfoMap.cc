#include "TaRunInfo.hh"

map<Int_t, TaRunInfo> LoadRunInfoMapBySlug(Int_t slug,Int_t arm_switch){
  map< Int_t, TaRunInfo > fRunInfoMap;  

  TString info_filename = Form("./prex-runlist/slug%d_info.list",slug);
  ifstream slug_info;
  slug_info.open(info_filename.Data());

  if(!slug_info.is_open()){
    cout << " slug info " << info_filename << " does not exist." << endl;
    return fRunInfoMap;
  }
    
  TString sline;

  while(sline.ReadLine(slug_info)){
    TaRunInfo fRunInfo;
    if(!fRunInfo.ParseLine(sline))
      continue;
    if(fRunInfo.GetRunFlag()!="Good")
      continue;
    if(fRunInfo.GetArmFlag()!=arm_switch)
      continue;
    Int_t run_number = fRunInfo.GetRunNumber();
    fRunInfoMap[run_number]=fRunInfo;
  }
  slug_info.close();
  return fRunInfoMap;
}

TString LoadArmFlagCuts(map<Int_t,TaRunInfo> InfoMap, TString det_name){
  TString fRetCut="0";
  auto iter =  InfoMap.begin();
  while(iter!=InfoMap.end()){
    Int_t run_number = (*iter).first;
    TaRunInfo myRunInfo = (*iter).second;
    if(det_name.Contains("avg") || det_name.Contains("dd") ){
      if(myRunInfo.GetArmFlag()==0){
	if(fRetCut!="")
	  fRetCut+="||";
	fRetCut += Form("run==%d",run_number);
      }
    }else if (det_name.Contains("l")){
      if(myRunInfo.GetArmFlag()!=1){
	if(fRetCut!="")
	  fRetCut+="||";
	fRetCut += Form("run==%d",run_number);
      }

    }else if (det_name.Contains("r")){
      if(myRunInfo.GetArmFlag()!=2){
	if(fRetCut!="")
	  fRetCut+="||";
	fRetCut += Form("run==%d",run_number);
      }
    }
    iter++;
  }
  return fRetCut;
}

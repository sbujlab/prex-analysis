/*
   author: Tao Ye
   Last update: July 2019
 */
map<Int_t, Int_t> LoadArmFlag(Int_t slug); // return <run, arm_flag>
map<Int_t, Int_t> LoadRunFlag(Int_t slug); // return <run, run_flag> : (1)Good (2)Tgt (-1) others
map< pair<Int_t,Int_t> , vector<Double_t> >GetLagrangeSlope(Int_t slug); // return <run,burst> => slope
map< pair<Int_t,Int_t> , vector<Double_t> >GetRegressionSlope(Int_t slug); // return <run,burst> => slope
map< pair<Int_t,Int_t> , Int_t > LoadCyc2BurstMap(Int_t slug, Int_t run); // return <cyc,coil> =>BurstCounter
map<Int_t,vector<Int_t> > LoadBadCycleList(); // return <cycle> => coil array;
Bool_t InBadCycleList(map<Int_t,vector<Int_t> > fMap, Int_t cyc, Int_t coil);

map<Int_t, Int_t> LoadArmFlag(Int_t slug){ // return <run, arm_flag>
  map<Int_t,Int_t> fret;
  TString info_filename =Form("prex-runlist/slug%d_info.list",slug);
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

map<Int_t, Int_t> LoadRunFlag(Int_t slug){ // return <run, run_flag> : (1)Good (2)Tgt (-1) others
  map<Int_t,Int_t> fret;
  TString info_filename =Form("prex-runlist/slug%d_info.list",slug);
  ifstream slug_info;
  slug_info.open(info_filename.Data());
  TString sline;
  while(sline.ReadLine(slug_info)){
    TObjArray *token = sline.Tokenize(',');
    Int_t run_number = (((TObjString*)(token->At(0)))->GetString()).Atoi();
    TString run_flag = ((TObjString*)(token->At(3)))->GetString();
    if(run_flag=="Good")
      fret[run_number]=1;
    else if(run_flag=="TargetDegradation")
      fret[run_number]=2;
    else
      fret[run_number]=-1;
  }
  slug_info.close();
  return fret;

}

map<Int_t,vector<Int_t> > LoadBadCycleList(){
  map<Int_t,vector<Int_t > > fRet;
  FILE *black_list;
  black_list= fopen("badcycle.list","r");

  if(black_list==NULL){
    cerr << " -- Error: bad cycle list is not found " << endl;
    return fRet;
  }
  while(!feof(black_list)){
    Int_t cycle=0;
    Int_t coil=0;
    fscanf(black_list,"%d,%d\n",&cycle,&coil);
    if(cycle!=0){
      if(coil==-1)
	for(int i=1;i<=7;i++)
	  fRet[cycle].push_back(i);
      else
	fRet[cycle].push_back(coil);
    }
  }
  fclose(black_list);
  return fRet;

}

Bool_t InBadCycleList(map<Int_t,vector<Int_t> > fMap, Int_t cycID, Int_t coil_idx=-1){

  if(fMap.find(cycID)==fMap.end())
    return kFALSE;
  
  if(coil_idx==-1) // only checking cycle number
    return kTRUE;
  vector<Int_t> coil_list = fMap[cycID];
  if(find(coil_list.begin(),coil_list.end(),coil_idx) == coil_list.end())
    return kFALSE;
  else
    return kTRUE;
}

map< pair<Int_t,Int_t>, vector<Double_t> > GetRegressionSlope(Int_t slug){
  map< pair<Int_t,Int_t>, vector<Double_t> > fSlopeMap;
  vector<TString> fDet_name = { "usl","usr","us_avg","us_dd"};
  vector<TString> fBPM_name;
  vector<TString> fBPM_name1 = {"11X","11Y","12X","12Y","16X","16Y",
				"1X","1Y","4aX","4aY","4eX","4eY"};
  vector<TString> fBPM_name2 = {"8X","8Y","12X","12Y",
				"1X","1Y","4aX","4aY","4eX","4eY"};
  if(slug<=2)
    fBPM_name = fBPM_name2;
  else
    fBPM_name = fBPM_name1;

  TFile *mini_merge = TFile::Open(Form("../lagrange_residual/respin2_treeMergeOutput/MergedLagrange_slug%d.root",slug));
  TTree* mini_tree =(TTree*) mini_merge->Get("mini_regall");
  mini_tree->AddFriend("mini");
  Int_t nDet = fDet_name.size();
  Int_t nBPM = fBPM_name.size();
  Int_t kRun, kBurst;
  vector<Double_t> fSlopes(nDet*nBPM,0);
  for(int idet=0;idet<nDet;idet++){
    for(int ibpm=0;ibpm<nBPM;ibpm++){
      mini_tree->SetBranchAddress(Form("%s_bpm%s",fDet_name[idet].Data(),fBPM_name[ibpm].Data()),
				  &fSlopes[idet*nBPM+ibpm]);
    }
  }
  mini_tree->SetBranchAddress("run",&kRun);
  mini_tree->SetBranchAddress("mini",&kBurst);
  Int_t nEntries = mini_tree->GetEntries();
  for(int ievt=0;ievt<nEntries;ievt++){
    mini_tree->GetEntry(ievt);
    fSlopeMap[ make_pair(kRun,kBurst) ] = fSlopes;
  }
  mini_merge->Close();
  return fSlopeMap;
}

map< pair<Int_t,Int_t>, vector<Double_t> > GetLagrangeSlope(Int_t slug){
  map< pair<Int_t,Int_t>, vector<Double_t> > fSlopeMap;
  vector<TString> fDet_name = { "usl","usr","us_avg","us_dd"};
  vector<TString> fBPM_name;
  vector<TString> fBPM_name1 = {"11X","11Y","12X","12Y","16X","16Y",
				"1X","1Y","4aX","4aY","4eX","4eY"};
  vector<TString> fBPM_name2 = {"8X","8Y","12X","12Y",
				"1X","1Y","4aX","4aY","4eX","4eY"};
  if(slug<=2)
    fBPM_name = fBPM_name2;
  else
    fBPM_name = fBPM_name1;

  TFile *mini_merge = TFile::Open(Form("../lagrange_residual/respin2_treeMergeOutput/MergedLagrange_slug%d.root",slug));
  TTree* mini_tree =(TTree*) mini_merge->Get("mini_lagrall");
  mini_tree->AddFriend("mini");
  Int_t nDet = fDet_name.size();
  Int_t nBPM = fBPM_name.size();
  Int_t kRun, kBurst;
  vector<Double_t> fSlopes(nDet*nBPM,0);
  for(int idet=0;idet<nDet;idet++){
    for(int ibpm=0;ibpm<nBPM;ibpm++){
      mini_tree->SetBranchAddress(Form("%s_bpm%s",fDet_name[idet].Data(),fBPM_name[ibpm].Data()),
				  &fSlopes[idet*nBPM+ibpm]);
    }
  }
  mini_tree->SetBranchAddress("run",&kRun);
  mini_tree->SetBranchAddress("mini",&kBurst);
  Int_t nEntries = mini_tree->GetEntries();
  for(int ievt=0;ievt<nEntries;ievt++){
    mini_tree->GetEntry(ievt);
    fSlopeMap[ make_pair(kRun,kBurst) ] = fSlopes;
  }
  mini_merge->Close();
  return fSlopeMap;
}

map< pair<Int_t,Int_t> , Int_t > LoadCyc2BurstMap(Int_t slug, Int_t run){
  map< pair<Int_t,Int_t> , Int_t > fRet; // <cycle,coil> -> burst;
  TFile* input = TFile::Open(Form("../lagrange_residual/rootfiles/slug%d_cycle_burst.root",slug));
  TTree* cyc2b_tree = (TTree*)input->Get("cyc2b");
  
  Int_t npt = cyc2b_tree->Draw("bmwcycnum:bmwobj:BurstCounter",Form("run==%d",run),"goff");
  double *fCycle = cyc2b_tree->GetV1();
  double *fCoil = cyc2b_tree->GetV2();
  double *fBurst = cyc2b_tree->GetV3();
  
  for(int ipt=0;ipt<npt;ipt++){
    Int_t kCyc = fCycle[ipt];
    Int_t kCoil = fCoil[ipt];
    Int_t kBurst = fBurst[ipt];
    fRet[ make_pair(kCyc,kCoil)] = kBurst;
  }
  
  input->Close();
  return fRet;
}

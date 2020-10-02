map<Int_t, vector<Int_t> > LoadRunInfo();

void GenerateRunInfoTreeBySlug(int slug){
  map<Int_t , vector<Int_t> > fRunInfo  = LoadRunInfo();  
  TFile* input = TFile::Open(Form("rootfiles/slug%d_sorted_eigenvector_allbpm.root",slug));
  TTree* parent_tree = (TTree*)input->Get("eig");
  TFile* output = TFile::Open(Form("rootfiles/slug%d_runinfo.root",slug),"RECREATE");
  TTree *info_tree = new TTree("mini_info"," Run Information by minirun ");
  Int_t run, burst_counter, arm_flag;
  Int_t kHelicitySign;
  Int_t isGoodRun;
  info_tree->Branch("slug",&slug);
  info_tree->Branch("arm_flag",&arm_flag);
  info_tree->Branch("run",&run);
  info_tree->Branch("mini",&burst_counter);
  info_tree->Branch("spin",&kHelicitySign);
  info_tree->Branch("kGood",&isGoodRun);
  
  parent_tree->SetBranchAddress("run",&run);
  parent_tree->SetBranchAddress("mini",&burst_counter);
  
  Int_t nevt = parent_tree->GetEntries();
  for(int ievt=0;ievt<nevt;ievt++){
    parent_tree->GetEntry(ievt);
    vector<Int_t> myInfo = fRunInfo[run];
    if(myInfo.size()==0){
      cout << " **  Run Info not found for run " <<run << ", skip " << endl;
      continue;
    }
    slug = myInfo[0];
    isGoodRun = myInfo[1];
    arm_flag = myInfo[2];
    kHelicitySign = myInfo[3];
    info_tree->Fill();
  } // end of minirun loop 
  
  input->Close();
  output->cd();
  info_tree->Write();
  output->Close();
}


map<Int_t, vector<Int_t> > LoadRunInfo(){
  map<Int_t, vector<Int_t> > fRet;
  TString info_filename ="../runlist/all_production.list";
  ifstream slug_info;
  slug_info.open(info_filename.Data());
  TString sline;
  while(sline.ReadLine(slug_info)){
    vector<Int_t> myInfo;
    TObjArray *token = sline.Tokenize(',');
    Int_t run_number = (((TObjString*)(token->At(0)))->GetString()).Atoi();
    Int_t slug_number = (((TObjString*)(token->At(1)))->GetString()).Atoi();
    myInfo.push_back(slug_number);
    TString run_flag = ((TObjString*)(token->At(3)))->GetString();
    if(run_flag=="Good")
      myInfo.push_back(1);
    else if(run_flag=="TargetDegradation")
      myInfo.push_back(2);
    else
      myInfo.push_back(0);
    Int_t arm_flag = (((TObjString*)(token->At(6)))->GetString()).Atoi();
    myInfo.push_back(arm_flag);
    TString ihwp = ((TObjString*)(token->At(4)))->GetString();
    TString wien = ((TObjString*)(token->At(5)))->GetString();
    Int_t kHelicitySign;
    if( ihwp == "IN" )
      kHelicitySign = 1;
    else if (ihwp == "OUT" )
      kHelicitySign =-1;
    if( wien == "FLIP-RIGHT" )
      kHelicitySign *= 1;
    else if ( wien == "FLIP-LEFT" )
      kHelicitySign *=-1;
    
    myInfo.push_back(kHelicitySign);
    fRet[run_number] = myInfo;
  }
  slug_info.close();
  return fRet;
}

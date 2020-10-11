map<Int_t, vector<Int_t> > LoadRunInfo();
void SaveMulHistogram(Int_t run_start=3897, Int_t run_end = 4864){
  TString label = Form("run%d-%d",run_start,run_end);
  TString filename  = label+".list";
  TString path = "../runlist/good_list/";
  FILE *runlist ;
  filename = "all.list";
  runlist = fopen((path+filename).Data(),"r");
  map<Int_t,vector<Int_t> > fInfoMap = LoadRunInfo();
  TFile *output = TFile::Open(Form("mulhist_%s.root",label.Data()),"RECREATE");
  TH1D *h1dreg = new TH1D("h1dreg","",1000,-3000,3000);
  TH1D *h1dlagr = new TH1D("h1dlagr","",1000,-3000,3000);
  TH1D *h1dregdd = new TH1D("h1dregdd","",1000,-1e3,1e3);
  TH1D *h1dlagrdd = new TH1D("h1dlagrdd","",1000,-1e3,1e3);

  gDirectory->ls();
  while(!feof(runlist)){
    Int_t run_number=0;
    Int_t seg_number=0;
    fscanf(runlist,"%d/n",&run_number);
    if(run_number==0 || run_number<run_start )
      continue;
    if(run_number>run_end)
      break;
    vector<Int_t> myInfo = fInfoMap[run_number];
    if(myInfo.size()!=4){
      cout << " ** run"<<run_number<<" info is not found and  will skip" << endl;
      continue;
    }
    if(myInfo[1]!=1){
      cout << " ** run" << run_number << " is not a Good run and will skip" << endl;
      continue;
    }
    if(myInfo[2]!=0){
      cout << " ** run" << run_number << " is a Single Arm run and will skip" << endl;
      continue;
    }
    Int_t kSign = myInfo[3];
    // TString qw_path ="/lustre19/expphy/volatile/halla/parity/prex-respin2/japanOutput/";
    // TString rootfile_name = Form("prexPrompt_pass2_%d.%03d.root",
    // 				 run_number,seg_number);
    TString lagrange_path= "/lustre/expphy/volatile/halla/parity/prex-respin2/LagrangeOutput/rootfiles/";
    TString rootfile_name = Form("prexRespin2_lagrange_eigen_%d.%03d.root",
				 run_number,seg_number);
    TFile *this_file;
    if(gSystem->AccessPathName(lagrange_path+rootfile_name)==0){
      this_file = TFile::Open(lagrange_path+rootfile_name);
      // cout << " -- Reading ROOT file ";
      // cout << this_file->GetName() << endl;
      TTree *lagrall = (TTree*) this_file->Get("lagrall");
      lagrall->AddFriend("mul");
      TTree *regall = (TTree*) this_file->Get("regall");
      regall->AddFriend("mul");

      output->cd();
      Int_t npt_lagr = lagrall->Project("+h1dlagr",
					Form("%f*lagr_asym_us_avg/ppm",(double)kSign),
					"ErrorFlag==0");
      lagrall->Project("+h1dlagrdd",
		       Form("%f*lagr_asym_us_dd/ppm",(double)kSign),
		       "ErrorFlag==0");

      Int_t npt_reg= regall->Project("+h1dreg",
				     Form("%f*reg_asym_us_avg/ppm",(double)kSign),
				     "ErrorFlag==0");
      regall->Project("+h1dregdd",
		      Form("%f*reg_asym_us_dd/ppm",(double)kSign),
		      "ErrorFlag==0");
      cout << " -- Done with ROOT file ";
      cout << this_file->GetName() << endl;
      this_file->Close();
    } // end of file exists
  } // end of run number loop

  output->cd();
  h1dreg->Write();
  h1dlagr->Write();
  h1dregdd->Write();
  h1dlagrdd->Write();
  cout << " -- Closing " << output->GetName() << endl;
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

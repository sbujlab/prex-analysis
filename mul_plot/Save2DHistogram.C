map<Int_t, vector<Int_t> > LoadRunInfo();
void Save2DHistogram(Int_t run_start=3305, Int_t run_end = 4980){
  TString label = Form("run%d-%d",run_start,run_end);
  TString filename  = label+".list";
  TString path = "../runlist/good_list/";
  FILE *runlist ;
  filename = "all.list";
  runlist = fopen((path+filename).Data(),"r");
  map<Int_t,vector<Int_t> > fInfoMap = LoadRunInfo();


  TFile *output = TFile::Open(Form("mul2d_%s.root",label.Data()),"RECREATE");
  TH2D *h2dreg = new TH2D("h2dreg","",run_end-run_start+1,run_start-0.5,run_end+0.5,1000,-3000,3000);
  TH2D *h2dlagr = new TH2D("h2dlagr","",run_end-run_start+1,run_start-0.5,run_end+0.5,1000,-3000,3000);
  TH2D *h2dregdd = new TH2D("h2dregdd","",run_end-run_start+1,run_start-0.5,run_end+0.5,1000,-1e3,1e3);
  TH2D *h2dlagrdd = new TH2D("h2dlagrdd","",run_end-run_start+1,run_start-0.5,run_end+0.5,1000,-1e3,1e3);

  Int_t slug_start =0;
  slug_start= fInfoMap[run_start][0];
  Int_t slug_end=0;
  slug_end  = fInfoMap[run_end][0];
  if(slug_start==0)
    slug_start = 1;
  if(slug_end==0)
    slug_end = 94;
  Int_t nslug = slug_end-slug_start+1;
  TH2D *h2dreg_slug = new TH2D("h2dreg_slug","",nslug,slug_start-0.5,slug_end+0.5,1000,-3000,3000);
  TH2D *h2dlagr_slug = new TH2D("h2dlagr_slug","",nslug,slug_start-0.5,slug_end+0.5,1000,-3000,3000);
  TH2D *h2dregdd_slug = new TH2D("h2dregdd_slug","",nslug,slug_start-0.5,slug_end+0.5,1000,-1e3,1e3);
  TH2D *h2dlagrdd_slug = new TH2D("h2dlagrdd_slug","",nslug,slug_start-0.5,slug_end+0.5,1000,-1e3,1e3);

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
    Int_t slug_number = myInfo[0];
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

      output->cd(); // Change directory so that histogram name can be identified
      lagrall->Project("+h2dlagr",
		       Form("%f*lagr_asym_us_avg/ppm:%d",(double)kSign,run_number),
		       "ErrorFlag==0");
      lagrall->Project("+h2dlagrdd",
      		       Form("%f*lagr_asym_us_dd/ppm:%d",(double)kSign,run_number),
      		       "ErrorFlag==0");
      regall->Project("+h2dreg",
		      Form("%f*reg_asym_us_avg/ppm:%d",(double)kSign,run_number),
		      "ErrorFlag==0");
      regall->Project("+h2dregdd",
      		      Form("%f*reg_asym_us_dd/ppm:%d",(double)kSign,run_number),
      		      "ErrorFlag==0");
      // by Slugs
      lagrall->Project("+h2dlagr_slug",
		       Form("%f*lagr_asym_us_avg/ppm:%d",(double)kSign,slug_number),
		       "ErrorFlag==0");
      lagrall->Project("+h2dlagrdd_slug",
      		       Form("%f*lagr_asym_us_dd/ppm:%d",(double)kSign,slug_number),
      		       "ErrorFlag==0");
      regall->Project("+h2dreg_slug",
		      Form("%f*reg_asym_us_avg/ppm:%d",(double)kSign,slug_number),
		      "ErrorFlag==0");
      regall->Project("+h2dregdd_slug",
      		      Form("%f*reg_asym_us_dd/ppm:%d",(double)kSign,slug_number),
      		      "ErrorFlag==0");

      cout << " -- Done with ROOT file ";
      cout << this_file->GetName() << endl;
      this_file->Close();
    } // end of file exists
  } // end of run number loop

  output->cd();
  h2dreg->Write();
  h2dlagr->Write();
  h2dregdd->Write();
  h2dlagrdd->Write();
  h2dreg_slug->Write();
  h2dlagr_slug->Write();
  h2dregdd_slug->Write();
  h2dlagrdd_slug->Write();

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

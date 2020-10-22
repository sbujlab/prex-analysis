map<Int_t, vector<Int_t> > LoadRunInfo();
void Save2DHistogramByList(TString list_name);

void Save2DHistogramByList(){
  TString list_array[]={"50uA","50uA_1arm"};
			// "70uA_120Hz","70uA_120Hz_1arm",
			// "70uA_240Hz","70uA_240Hz_1arm",
			// "85uA"};
  Int_t npt=  sizeof(list_array)/sizeof(*list_array);
  for(int i=0;i<npt;i++){
    Save2DHistogramByList(list_array[i]);
  }
}
void Save2DHistogramByList(TString list_name){
  TString label = list_name;
  FILE *runlist = fopen(("./list/"+list_name).Data(),"r") ;
  
  map<Int_t,vector<Int_t> > fInfoMap = LoadRunInfo();
  TFile *output = TFile::Open(Form("mul2d_%s.root",label.Data()),"RECREATE");
  Int_t run_start = 3305;
  Int_t run_end = 4980;

  TH2D *h2dreg = new TH2D("h2dreg","",run_end-run_start+1,run_start-0.5,run_end+0.5,1000,-3000,3000);
  TH2D *h2dlagr = new TH2D("h2dlagr","",run_end-run_start+1,run_start-0.5,run_end+0.5,1000,-3000,3000);
  TH2D *h2dbcm = new TH2D("h2dbcm","",run_end-run_start+1,run_start-0.5,run_end+0.5,1000,0,100);

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
  TH2D *h2dbcm_slug = new TH2D("h2dbcm_slug","",nslug,slug_start-0.5,slug_end+0.5,1000,0,100);
  
  gDirectory->ls();
  while(!feof(runlist)){
    Int_t run_number=0;
    Int_t seg_number=0;
    fscanf(runlist,"%d/n",&run_number);
    if(run_number==0)
      continue;
    vector<Int_t> myInfo = fInfoMap[run_number];
    if(myInfo.size()!=4){
      cout << " ** run"<<run_number<<" info is not found and  will skip" << endl;
      continue;
    }
    if(myInfo[1]!=1){
      cout << " ** run" << run_number << " is not a Good run and will skip" << endl;
      continue;
    }
    Int_t arm_flag = myInfo[2];
    TString det_name = "asym_us_avg";
    if(arm_flag==1)
      det_name = "asym_usr";
    else if(arm_flag==2)
      det_name = "asym_usl";
    
    Int_t kSign = myInfo[3];
    Int_t slug_number = myInfo[0];
    TString qw_path ="/lustre19/expphy/volatile/halla/parity/prex-respin2/japanOutput/";
    TString japan_rootfile_name = Form("prexPrompt_pass2_%d.%03d.root",
    				 run_number,seg_number);
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

      TFile* japan_file = TFile::Open(qw_path+japan_rootfile_name);
      TTree *mul_tree = (TTree*)japan_file->Get("mul"); // for beam current
      
      output->cd(); // Change directory so that histogram name can be identified
      lagrall->Project("+h2dlagr",
		       Form("%f*lagr_%s/ppm:%d",(double)kSign,det_name.Data(),run_number),
		       "ErrorFlag==0");
      regall->Project("+h2dreg",
		      Form("%f*reg_%s/ppm:%d",(double)kSign,det_name.Data(),run_number),
		      "ErrorFlag==0");

      mul_tree->Project("+h2dbcm",
		      Form("yield_bcm_target:%d",run_number),
		      "ErrorFlag==0");

      // by Slugs
      lagrall->Project("+h2dlagr_slug",
		       Form("%f*lagr_%s/ppm:%d",(double)kSign,det_name.Data(),slug_number),
		       "ErrorFlag==0");
      regall->Project("+h2dreg_slug",
		      Form("%f*reg_%s/ppm:%d",(double)kSign,det_name.Data(),slug_number),
		      "ErrorFlag==0");
      mul_tree->Project("+h2dbcm_slug",
		      Form("yield_bcm_target:%d",slug_number),
		      "ErrorFlag==0");
      
      cout << " -- Done with ROOT file ";
      cout << this_file->GetName() << endl;
      this_file->Close();
      japan_file->Close();
    } // end of file exists
  } // end of run number loop

  output->cd();
  h2dreg->Write();
  h2dlagr->Write();
  h2dbcm->Write();
  h2dreg_slug->Write();
  h2dlagr_slug->Write();
  h2dbcm_slug->Write();

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

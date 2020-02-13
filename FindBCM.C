// Normalizing BCM Searching for JAPAN Output
// Author - Tao Ye
#include "TaAccumulator.cc"
TString device_list[]={"bcm_an_us","bcm_an_ds",
		       "bcm_an_ds3","bcm_an_ds10",
		       "bcm_dg_us","bcm_dg_ds",
		       "unser"};
void FindBCM(TString listname="all_production.list"){
  TStopwatch tsw;
  Int_t nthreads=10;
  ROOT::EnableImplicitMT(nthreads);
  TFile *output = TFile::Open("prex_bcm_status.root","RECREATE");
  TTree *run_tree = new TTree("T","BCM run status Tree");
  Int_t ndev = sizeof(device_list)/sizeof(*device_list);

  TString detector="sam1";
  std::vector<Double_t> bcm_value(ndev);
  std::vector<Int_t> bcm_flag(ndev);
  typedef struct{Double_t hw_sum,raw,nsamples;} LeafList;
  LeafList detector_leaf;
  Int_t run_number;
  Int_t slug_number;
  Int_t bcm_index;
  Double_t best_chi;
  run_tree->Branch("run",&run_number);
  run_tree->Branch("slug",&slug_number);
  run_tree->Branch("bcm_index",&bcm_index);
  run_tree->SetMarkerStyle(20);
  for(int idev=0;idev<ndev;idev++){
    run_tree->Branch(device_list[idev]+"_flag",&(bcm_flag[idev]));
    run_tree->Branch("best_chi",&best_chi);
  }

  Int_t counts=0;
  FILE *runlist = fopen(listname,"r");
  ofstream error_log;
  error_log.open("find_bcm.error");
  while(!feof(runlist)){
    TStopwatch tsw_run;
    fscanf(runlist,"%d\t%d\n",&run_number,&slug_number);
    TFile *input_rf = TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass1_%d.000.root",run_number));
    if(input_rf==NULL){
      error_log << " -- skip " << run_number  << endl;
      error_log << " -- rootfile does not exsits "  << endl;
      cout << " -- skip " << run_number  << endl;
      cout << " -- rootfile does not exsits "  << endl;
      continue;
    }
    TTree* evt_tree = (TTree*)input_rf->Get("evt");
    if(evt_tree==NULL){
      error_log << " -- skip " << run_number  << endl;
      error_log << " -- evt tree not found"  << endl;
      cout << " -- skip " << run_number  << endl;
      cout << " -- evt tree not found"  << endl;

      continue;
    }
    cout << " -- Processing run  " << run_number ;
    if(run_number>=3883 && run_number<=3893){
      detector="usl";
      cout << " -- Using usl for this run " << endl ;
    }
    else
      detector="sam1";
    evt_tree->SetBranchStatus("*",0);
    evt_tree->SetBranchStatus("ErrorFlag",1);
    evt_tree->SetBranchStatus(detector,1);
    Double_t kErrorFlag;
    evt_tree->SetBranchAddress("ErrorFlag",&kErrorFlag);
    TBranch *det_branch = evt_tree->GetBranch(detector);
    det_branch->GetLeaf("hw_sum")->SetAddress(&detector_leaf.hw_sum);
    det_branch->GetLeaf("hw_sum_raw")->SetAddress(&detector_leaf.raw);
    det_branch->GetLeaf("num_samples")->SetAddress(&detector_leaf.nsamples);
    for(int idev=0;idev<ndev;idev++){
      TString chname = device_list[idev];
      TBranch *this_branch = evt_tree->GetBranch(chname);
      if(this_branch!=NULL){
	evt_tree->SetBranchStatus(chname,1);
	this_branch->GetLeaf("hw_sum")->SetAddress(&bcm_value[idev]);
	bcm_flag[idev]= 0;
      }
      else{
	bcm_value[idev]=0;
	bcm_flag[idev]=-1;
      }
    }
    AccumulatorArray fChiAccumulator(ndev);  // Chi Covariance Array
    AccumulatorArray fDetAccumulator(ndev);  // Re-normalized Detector Variance Array
    AccumulatorArray fBCMAccumulator(ndev);  // BCM Variance
    TaAccumulator fRawVar;  // Raw Detector Volt Variance

    Bool_t kSkip=kFALSE;
    Int_t ncounts=0;
    Int_t ievt =0;
    while(ncounts<20){
      Int_t kStatus=evt_tree->GetEntry(ievt);
      if(kStatus<=0){
	error_log << " -- skip " << run_number << endl;
	error_log << " -- Not enough events in this run " << endl;
	cout << " -- skip " << run_number << endl;
	cout << " -- Not enough events in this run " << endl;
	kSkip = kTRUE;
	break;
      }
      ievt++;
      Int_t nGoodBCM =0;
      for(int idev=0;idev<ndev;idev++){
	if(bcm_value[idev]>10){
	  nGoodBCM++;
	}
      }
      if(nGoodBCM>=3)
	ncounts++;
      else
	continue;
      Double_t raw_volt = (detector_leaf.raw/detector_leaf.nsamples)*7.6e-5;
      Double_t yield = detector_leaf.hw_sum; // still normalized here 
      fRawVar.Update(raw_volt);
      for(int idev=0;idev<ndev;idev++){
	if(bcm_flag[idev]==0){
	  Double_t current = bcm_value[idev];
	  Double_t yield_in_volt = yield*current;
	  fChiAccumulator[idev].Update(yield_in_volt-raw_volt);
	  fDetAccumulator[idev].Update(yield_in_volt);
	  fBCMAccumulator[idev].Update(current);
	}
      }
    } // end of event loop
    if(!kSkip){
      best_chi = 1e6;
      bcm_index = -1;
      for(int idev=0;idev<ndev;idev++){
	if(bcm_flag[idev]==0
	   && fBCMAccumulator[idev].GetMean1()>2.5){
	  Double_t this_chi = fChiAccumulator[idev].GetM2();
	  if(this_chi < best_chi){
	    bcm_index = idev;
	    best_chi =this_chi;
	  }
	}
      }
      if(bcm_index>=0){
	bcm_flag[bcm_index]=1;
	cout << " -- Found normalizing BCM " << device_list[bcm_index] << endl;
      }
      run_tree->Fill();
    }
    input_rf->Close();
    counts++;
    cout << " in real time " << tsw_run.RealTime()*1000 << " msec. " << endl;
  }
  output->cd();
  run_tree->Write();
  output->Close();
  error_log.close();
  cout << "Got " << counts << " runs" << endl;
  tsw.Print();
}

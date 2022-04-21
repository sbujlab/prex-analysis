/*
  author: Tao Ye <tao.ye@stonybrook.edu>
  last update: July 2020

 */
#include "utilities_eigen.cc"

void SortEigenvectors_CREX_superstabilize(TString treename = "mini_eigen_reg_allbpms", TString pass = "", Int_t nCheck = 5, Int_t min = 7499, Int_t max = 8559){ //(int start, int end)
  //map<Int_t , pair<Int_t, Int_t> > fRunInfo  = LoadRunInfo();
  //TString filename_tag = "allbpms";
  vector< TString > IVlist = {
      "bpm12X",
      "bpm11X",
      "bpm4eX",
      "bpm1X" ,
      "bpm4aX",
      "bpm16X",
      "bpm4aY",
      "bpm4eY",
      "bpm1Y" ,
      "bpm16Y",
      "bpm12Y",
      "bpm11Y",
  };
  if (treename.Contains("5bpms")) {
    IVlist = {
      "bpm12X",
      "bpm1X" ,
      "bpm4eX",
      "bpm4aY",
      "bpm4eY",
    };
  }

  vector<TString> DVlist ={"us_avg","us_dd","usl","usr"};//,"us_dd"};
  Int_t nDV = DVlist.size();
  Int_t nIV = IVlist.size();
  //// TFile* output = TFile::Open(Form("rootfiles/slug%d-%d_sorted_eigenvector_%s.root",
  //// 				   start,end,filename_tag.Data()),"RECREATE");
  TFile* output = TFile::Open(Form("rootfiles/%s_sorted_super%s.root",
  				   treename.Data(),pass.Data()),"RECREATE");

  //TTree *lagr_slope_tree = new TTree("lagr","lagrange slopes in eigenvectors basis");
  //TTree *reg_slope_tree = new TTree("reg","regression slopes in eigenvectors basis");
  TTree *eig_sort_tree = new TTree(treename+"_sorted", "resorted eigenvectors tree");
  //TTree *eig_sort_tree = new TTree("mini_eigen_reg_allbpms_tr", "resorted eigenvectors tree");
  //TTree *eig_sort_tree = new TTree("mini_eigen_reg_5bpms", "resorted eigenvectors tree");
  //vector<Double_t> fEigenLagrangeSlopes(nDV*nIV);
  //vector<Double_t> fEigenRegressionSlopes(nDV*nIV);
  vector<Double_t> fSign(nIV);
  vector<Double_t> fLocker(nIV*nIV,0);
  vector<Double_t> fEigenVectorSorted(nIV*nIV);
  vector<Double_t> fEigRegSlopeSorted(nDV*nIV);
  typedef struct {Double_t mean,err,rms,m2,nsamp;} STAT;
  vector<STAT> fEvMon(nIV);
  
  for(int j=0;j<nIV;j++)
    eig_sort_tree->Branch(Form("diff_evMon%d",j),&fEvMon[j],"mean/D:err:rms:m2:nsamp");
  
  for(int i=0;i<nIV;i++)
    for(int j=0;j<nIV;j++)
      eig_sort_tree->Branch(Form("evMon%d_%s",j, IVlist[i].Data()),
			 &fEigenVectorSorted[i*nIV+j]); // be carefull with row x col

  for(int i=0;i<nDV;i++)
    for(int j=0;j<nIV;j++){
      eig_sort_tree->Branch(DVlist[i]+Form("_evMon%d",j),
          &fEigRegSlopeSorted[i*nIV+j]); // be carefull with row x col
      //lagr_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			      //&fEigenLagrangeSlopes[i*nIV+j]);
      //reg_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			      //&fEigenRegressionSlopes[i*nIV+j]);

    }

  //TChain *agg_tree = new TChain("agg");
  TChain *mini_tree       = new TChain("mini","original mini tree");
  TChain *eig_all_tree    = new TChain(treename,"original eigen regression tree");
  ////TChain *eig_all_tr_tree = new TChain("mini_eigen_reg_allbpms_tr");
  ////TChain *eig_5_tree      = new TChain("mini_eigen_reg_5bpms");
  ////TChain *eig_5_tr_tree   = new TChain("mini_eigen_reg_5bpms_tr");
  //TChain *lagr_tree = new TChain("mini_lagrall");
  //TChain *reg_tree = new TChain("mini_regall");
  //for(int i=start;i<=end;i++){
  //agg_tree->Add(Form("dataRootfiles/CREX-All-miniruns.root"));
  //if (pass != "") {
    //mini_tree       -> Add(Form("dataRootfiles/rcdb_eigenvectors_sorted.root"));
    //eig_all_tree    -> Add(Form("dataRootfiles/rcdb_eigenvectors_sorted.root"));
    mini_tree       -> Add(Form("rootfiles/rcdb_eigenvectors%s.root",pass.Data()));
    eig_all_tree    -> Add(Form("rootfiles/rcdb_eigenvectors%s.root",pass.Data()));
  //}
  //else {
    ////mini_tree       -> Add(Form("dataRootfiles/rcdb_eigenvectors.root"));
    ////eig_all_tree    -> Add(Form("dataRootfiles/rcdb_eigenvectors.root"));
    //mini_tree       -> Add(Form("rootfiles/rcdb_eigenvectors.root"));
    //eig_all_tree    -> Add(Form("rootfiles/rcdb_eigenvectors.root"));
  //}
  ////eig_all_tr_tree -> Add(Form("dataRootfiles/rcdb_eigenvectors.root"));
  ////eig_5_tree      -> Add(Form("dataRootfiles/rcdb_eigenvectors.root"));
  ////eig_5_tr_tree   -> Add(Form("dataRootfiles/rcdb_eigenvectors.root"));

  //eig_all_tree->AddFriend(lagr_tree);
  //eig_all_tree->AddFriend(reg_tree);
  eig_all_tree->AddFriend(mini_tree);
  //eig_all_tree->AddFriend(eig_all_tr_tree);
  //eig_all_tree->AddFriend(eig_5_tree);
  //eig_all_tree->AddFriend(eig_5_tr_tree);

  vector<STAT> fEvMonRaw(nIV);
  vector<Double_t> fEigenVectorRaw(nIV*nIV);
  vector<Double_t> fEigRegSlope(nDV*nIV);
  //vector<Double_t> fLagrangeSlope(nDV*nIV);
  //vector<Double_t> fRegressionSlope(nDV*nIV);  

  Int_t run, burst_counter, arm_flag, slug, run_type, run_flag;
  mini_tree->SetBranchAddress("run",&run);
  mini_tree->SetBranchAddress("mini",&burst_counter);
  mini_tree->SetBranchAddress("rcdb_slug",&slug);
  mini_tree->SetBranchAddress("rcdb_arm_flag",&arm_flag);
  mini_tree->SetBranchAddress("rcdb_run_type",&run_type);
  mini_tree->SetBranchAddress("rcdb_run_flag",&run_flag);

  // Adding new branches?
  //eig_sort_tree->Branch();
  //lagr_slope_tree->Branch();
  //reg_slope_tree->Branch();

  for(int i=0;i<nIV;i++) {
    eig_all_tree->SetBranchAddress(Form("diff_evMon%d",i),
        &fEvMonRaw[i]);
  }
  for(int i=0;i<nIV;i++) {
    for(int j=0;j<nIV;j++) {
      eig_all_tree->SetBranchAddress(Form("evMon%d_%s",j,IVlist[i].Data()),
          &fEigenVectorRaw[i*nIV+j]);
    }
  }
  // [row]:bpm ; [col]:eigv
  Int_t local_run  = 0;
  Int_t local_mini = 0;
  eig_all_tree->SetBranchAddress(Form("mini.run"), &local_run);
  eig_all_tree->SetBranchAddress(Form("mini.mini"),&local_mini);
  for(int i=0;i<nDV;i++) {
    for(int j=0;j<nIV;j++) {
      eig_all_tree->SetBranchAddress(DVlist[i]+Form("_evMon%d",j),
         &fEigRegSlope[i*nIV+j]);
      //lagr_tree->SetBranchAddress(DVlist[i]+"_"+IVlist[j],
      //   &fLagrangeSlope[i*nIV+j]);
      //reg_tree->SetBranchAddress(DVlist[i]+"_"+IVlist[j],
         //&fRegressionSlope[i*nIV+j]);
    }
  }

  Int_t nEntries = eig_all_tree->GetEntries();
  Int_t last_run = -1;

  vector< vector<Double_t> > fEVRing; // [RingSize] x [nIV*nIV]
  vector< Double_t > fRingAvg(nIV*nIV,0);

  //eig_all_tree->GetEntry(0);
  //eig_all_tree->GetEntry(810); // Run 6036 starting guess
  //fEVRing.push_back(fEigenVectorRaw);
  /*
  for (Int_t ll = 363 ; ll <= 1315 ; ll++) {
    eig_all_tree->GetEntry(ll);
    fEVRing.push_back(fEigenVectorRaw);
    fRingAvg = GetRingAverage(fEVRing);
  }
  */
  //fRingAvg = GetRingAverage(fEVRing);

  //for(int ievt=nEntries-1;ievt>=0;ievt--)
  for(int ievt=0;ievt<nEntries;ievt++){
    eig_all_tree->GetEntry(ievt);
    if (local_run == 6328 || local_run == 7626) {
      fEVRing.clear();
      if (local_run == 6328) {
        eig_all_tree->GetEntry(1514); // Run 6516 = entry 1514... stable point
      }
      if (local_run == 7626) {
        eig_all_tree->GetEntry(6670); // Run 8003
      }
      fEVRing.push_back(fEigenVectorRaw);
      fRingAvg = GetRingAverage(fEVRing);
      eig_all_tree->GetEntry(ievt);
    }
    //if (local_run > min && local_run < max) {
    //////// lagr_tree->GetEntry(ievt);
    if(run!=last_run){
      last_run  = run;
      cout << " -----  " << endl;
      cout << " -- Run " << run << endl;
      cout << " -----  " << endl;
    }
    // Load raw BPM slopes to matrix
    //TMatrixD raw_eigreg_slope_matrix(nIV,nDV);
    //TMatrixD raw_lagr_slope_matrix(nIV,nDV);
    //TMatrixD raw_reg_slope_matrix(nIV,nDV);
    /*for(int i=0;i<nDV;i++) {
      for(int j=0;j<nIV;j++) {
        //raw_eigreg_slope_matrix[j][i] = fEigRegSlope[i*nIV+j];
        //raw_lagr_slope_matrix[j][i] = fLagrangeSlope[i*nIV+j];
        //raw_reg_slope_matrix[j][i] = fRegressionSlope[i*nIV+j];
      }
    }*/

    // Checking eigenvector identities 
    // re-allocation and sign lock-in to  fEigenVectorSorted

    Int_t thirdevt = ievt;
    Int_t fourthevt = ievt;
    if (local_run < 6328) {
      thirdevt = 1315;
      Printf("Reset to do backwards loop");
      fourthevt = eig_all_tree->GetEntry(thirdevt); // Reset's local_run == the 1315 event = run 6321
      fEVRing.clear();
      fEVRing.push_back(fEigenVectorRaw);
      fRingAvg = GetRingAverage(fEVRing);
    }
    // Loop if thirdevt == ievt       -> we are in normal loop
    // Loop if local_run != ievt  -> we are in backwards first section loop
    while ( fourthevt != ievt || thirdevt == ievt) {
    if (local_run < 6328 && (fourthevt!=ievt && thirdevt!=ievt)) {
      // We are in backwards first section loop
      Printf("Event %d",thirdevt);
      thirdevt--;
      fourthevt = eig_all_tree->GetEntry(thirdevt);
    }
    else {
      thirdevt++;
      fourthevt = ievt;
    }
      Printf("Event %d",thirdevt);

    vector<Int_t> fMap = CheckStabilizeIdentity(fRingAvg, fEigenVectorRaw, nCheck);
    vector<Double_t> fReMapped = RemapVectors(fEigenVectorRaw, fMap);
    // NEW: Need to do fMap realignment for Eigen Vector contents and for Eigen Vector Slopes contents
    Printf("Test 1");
    vector<Double_t> fReMappedEigRegSlopes = RemapColumns(fEigRegSlope, fMap, nIV);
    Printf("Test 2");
    vector<Int_t> fSign = CheckSignWithRingAvg( fRingAvg, fReMapped);
    Printf("Test 3");
    vector<Double_t> fFlipped = FlipVectors(fReMapped, fSign); 
    Printf("Test 4");
    vector<Double_t> fSwappedEigRegSlopes = FlipColumns(fReMappedEigRegSlopes, fSign, nIV); 
    Printf("Test 5");


    for(int imon=0;imon<nIV;imon++){
      fEvMon[ fMap[imon] ].mean = fEvMonRaw[imon].mean * fSign[ fMap[imon] ];
      fEvMon[ fMap[imon] ].err = fEvMonRaw[imon].err;
      fEvMon[ fMap[imon] ].m2 = fEvMonRaw[imon].m2;
      fEvMon[ fMap[imon] ].rms = fEvMonRaw[imon].rms;
      fEvMon[ fMap[imon] ].nsamp = fEvMonRaw[imon].nsamp;
    }
    Printf("Test 6");


    for(int irow=0;irow<nIV;irow++)
      for(int icol=0;icol<nIV;icol++)
        fEigenVectorSorted[irow*nIV+icol]  = fFlipped[irow*nIV+icol];
    Printf("Test 7");

    for(int irow=0;irow<nDV;irow++)
      for(int icol=0;icol<nIV;icol++)
        fEigRegSlopeSorted[irow*nIV+icol] = fSwappedEigRegSlopes[irow*nIV+icol];
    Printf("Test 8");

    // Update Template and Ring averages after flipped and remapped. 
    if(ievt!=0){
      Printf("Adding onto Event Ring for further calculations");
      // For the first chunk of crex_part1 ignore the first 363 entries.... don't add to ring
      if (ievt > 363) {
        fEVRing.push_back(fEigenVectorSorted);
      }
      Printf("Averaging");
      fRingAvg = GetRingAverage(fEVRing);
      Printf("Returned");
    }
    Printf("Test 9");

    //  Remix Slopes in eigenvector 
    //TMatrixD eig_vec_matrix(nIV,nIV);   // [row]:bpm ; [col]:eigv
    //for(int i=0;i<nIV;i++)
      //for(int j=0;j<nIV;j++)
        //eig_vec_matrix[i][j] = fEigenVectorSorted[i*nIV+j];

    //TMatrixD inv_eig_vec(eig_vec_matrix);
    //inv_eig_vec.Invert();
    //TMatrixD eig_lagr_slope_matrix(nIV,nDV);
    //eig_lagr_slope_matrix = inv_eig_vec * raw_lagr_slope_matrix;

    //TMatrixD eig_reg_slope_matrix(nIV,nDV);
    //eig_reg_slope_matrix = inv_eig_vec * raw_reg_slope_matrix;

    //for(int i=0;i<nDV;i++)
      //for(int j=0;j<nIV;j++){
        //fEigenLagrangeSlopes[i*nIV+j] = eig_lagr_slope_matrix[j][i];
        //fEigenRegressionSlopes[i*nIV+j] = eig_reg_slope_matrix[j][i];
      //}

    //if( fRunInfo.find(run)!= fRunInfo.end()){
      //slug = (fRunInfo[run]).first;
      //arm_flag = (fRunInfo[run]).second;
      //lagr_slope_tree->Fill();
      //reg_slope_tree->Fill();
    //}

    }
    eig_sort_tree->Fill();
  }
  vector<Double_t> eigen_width(nIV);
  vector<Int_t> fGlobalRank(nIV);
  for(int iev=0;iev<nIV;iev++){
    eig_sort_tree->Draw(Form("diff_evMon%d.rms*1e3",iev),"","goff");
    TH1D* hptr = (TH1D*)gDirectory->FindObject("htemp");
    if (hptr) {
      eigen_width[iev] = hptr->GetMean();
    }
    else {
      eigen_width[iev] = 1.0e9;
      Printf("Eigen vector %d failed",iev);
    }
    fGlobalRank[iev] = iev;
  }
  for(int i=0;i<nIV;i++){
    for(int j=i+1;j<nIV;j++){
      if(eigen_width[j] > eigen_width[i] ){
        double swap_val = eigen_width[i];
        eigen_width[i] = eigen_width[j];
        eigen_width[j] = swap_val;
        double swap_idx = fGlobalRank[i];
        fGlobalRank[i] = fGlobalRank[j];
        fGlobalRank[j] = swap_idx;
      }
    }
  }
  cout << " -- EigenWidth : ";
  for(int i=0;i<nIV;i++){
    printf("%.1f,",eigen_width[i] );
  }
  cout << endl;
  cout << " -- fGlobal Rank : ";
  for(int i=0;i<nIV;i++){
    cout << fGlobalRank[i] << "," ;
  }
  cout << endl;
  vector<Int_t> fGlobalMap(nIV);
  for(int i=0;i<nIV;i++)
  fGlobalMap[ fGlobalRank[i] ] = i;

  fRingAvg  = RemapVectors(fRingAvg, fGlobalMap);
  TMatrixD skeleton_m(nIV,nIV);
  for(int irow=0;irow<nIV;irow++){
    for(int icol=0;icol<nIV;icol++){
      skeleton_m[irow][icol] = fRingAvg[irow*nIV+icol];
    }
  }

  cout << " -- Global Ring Averages " << endl;
  PrintVector(fRingAvg);
  output->cd();
  eig_sort_tree->Write();
  mini_tree->Write();
  eig_all_tree->Write();
  //lagr_slope_tree->Write();
  //reg_slope_tree->Write();
  output->Close();

  TFile *output_avg = TFile::Open(Form("%s_skeleton_matrix.root",treename.Data()),"RECREATE");
  output_avg->cd();
  output_avg->WriteObject( &fRingAvg, "fRingAvg");
  output_avg->WriteObject( &skeleton_m, "skeleton");
  output_avg->Close();
}

/*
  author: Tao Ye <tao.ye@stonybrook.edu>
  last update: Oct 2020

 */
#include "utilities_eigen.cc"

void SortEigenvectorsBySlug_updownDD(int slug, Bool_t kForcedAvg=kFALSE){
  // Note: use kForcedAvg for slug 10,21,39,49,50,60,63 ( ignore this ... )
  map<Int_t , pair<Int_t, Int_t> > fRunInfo  = LoadRunInfo();
  map<Int_t , Int_t > fSpinInfo  = LoadSpinInfo();
  TString filename_tag = "allbpm";
  vector< TString > IVlist =  {"bpm4aX","bpm4eX","bpm1X",
			       "bpm11X","bpm12X","bpm16X",
			       "bpm4aY","bpm4eY","bpm1Y",
			       "bpm11Y","bpm12Y","bpm16Y"};
  
  vector< TString > IVlist1 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm8X","bpm12X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm8Y","bpm12Y"};
  if(slug<=2){
    IVlist=IVlist1;
    // filename_tag = "10bpm";
  }
  vector<TString> DVlist ={"us_avg","usl","usr",
			   "ds_avg","dsl","dsr",
			   "lr_avg_dd","lr_dd_dd","l_dd","r_dd"};//FIXME later  

  Int_t nDV = DVlist.size();
  Int_t nIV = IVlist.size();
  TFile* output = TFile::Open(Form("rootfiles/slug%d_updownDD_sorted_eigenvector_%s.root",
  				   slug,filename_tag.Data()),"RECREATE");

  TTree *lagr_slope_tree = new TTree("lagr","lagrange slopes in eigenvectors basis");
  TTree *reg_slope_tree = new TTree("reg","regression slopes in eigenvectors basis");
  TTree *eig_sort_tree = new TTree("eig", "sorted eigenvectors tree");
  vector<Double_t> fEigenLagrangeSlopes(nDV*nIV);
  vector<Double_t> fEigenRegressionSlopes(nDV*nIV);
  vector<Double_t> fSign(nIV);
  vector<Double_t> fLocker(nIV*nIV,0);
  vector<Double_t> fEigenVectorSorted(nIV*nIV);
  typedef struct {Double_t mean,err,rms,nsamp;} STAT;
  vector<STAT> fEvMon(nIV);
  vector<STAT> fDetLagr(nDV);
  vector<STAT> fDetReg(nDV);
  Int_t run, burst_counter, arm_flag;
  Int_t kHelicitySign;
  eig_sort_tree->Branch("slug",&slug);
  eig_sort_tree->Branch("arm_flag",&arm_flag);
  eig_sort_tree->Branch("run",&run);
  eig_sort_tree->Branch("mini",&burst_counter);
  eig_sort_tree->Branch("spin",&kHelicitySign);
  for(int j=0;j<nIV;j++)
    eig_sort_tree->Branch(Form("diff_evMon%d",j),&fEvMon[j],"mean/D:err:rms:num_samples");

  for(int i=0;i<nIV;i++)
    for(int j=0;j<nIV;j++)
      eig_sort_tree->Branch(Form("evMon%d_%s",j, IVlist[i].Data()),
			 &fEigenVectorSorted[i*nIV+j]); // be carefull with row x col

  for(int i=0;i<nDV;i++){
    lagr_slope_tree->Branch("lagr_asym_"+DVlist[i], &fDetLagr[i],"mean/D:err:rms:num_samples");
    reg_slope_tree->Branch("reg_asym_"+DVlist[i], &fDetReg[i],"mean/D:err:rms:num_samples");
    for(int j=0;j<nIV;j++){
      lagr_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			      &fEigenLagrangeSlopes[i*nIV+j]);
      reg_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			      &fEigenRegressionSlopes[i*nIV+j]);

    }
  }
  TFile* input_matrix;
  if(slug>=3)
    input_matrix= TFile::Open("skeleton_matrix.root");
  else
    input_matrix= TFile::Open("skeleton_matrix_10BPM.root");
  
  TMatrixD* skeleton_matrix = (TMatrixD*) input_matrix->Get("skeleton");
  vector<Double_t> fRingAvg(nIV*nIV);
  for(int irow=0;irow<nIV;irow++)
    for(int icol=0;icol<nIV;icol++)
      fRingAvg[irow*nIV+icol] = (*skeleton_matrix)[irow][icol];
  
  TFile* input_merged = TFile::Open( Form("treeMergeOutput/MergedLagrange_updownDD_slug%d.root",slug));
  TTree *mini_tree = (TTree*)input_merged->Get("mini");
  TTree *reg_tree = (TTree*)input_merged->Get("mini_regall");
  TTree *lagr_tree = (TTree*)input_merged->Get("mini_lagrall");

  reg_tree->AddFriend(lagr_tree);
  reg_tree->AddFriend(reg_tree);
  reg_tree->AddFriend(mini_tree);

  vector<STAT> fEvMonRaw(nIV);
  vector<Double_t> fEigenVectorRaw(nIV*nIV);
  vector<Double_t> fLagrangeSlope(nDV*nIV);
  vector<Double_t> fRegressionSlope(nDV*nIV);  
  mini_tree->SetBranchAddress("run",&run);
  mini_tree->SetBranchAddress("mini",&burst_counter);
  
  for(int i=0;i<nIV;i++)
    reg_tree->SetBranchAddress(Form("diff_evMon%d",i),
  			       &fEvMonRaw[i]);
  
  for(int i=0;i<nIV;i++)
    for(int j=0;j<nIV;j++)
      reg_tree->SetBranchAddress(Form("evMon%d_%s",j,IVlist[i].Data()),
				 &fEigenVectorRaw[i*nIV+j]);
  // [row]:bpm ; [col]:eigv
  for(int i=0;i<nDV;i++)
    for(int j=0;j<nIV;j++){
      lagr_tree->SetBranchAddress(Form("%s_evMon%d",DVlist[i].Data(),j),
				  &fLagrangeSlope[i*nIV+j]);

      reg_tree->SetBranchAddress(Form("%s_evMon%d",DVlist[i].Data(),j),
      				  &fRegressionSlope[i*nIV+j]);
      lagr_tree->SetBranchAddress("lagr_asym_"+DVlist[i], &fDetLagr[i]);
      reg_tree->SetBranchAddress("reg_asym_"+DVlist[i], &fDetReg[i]);
    }
  
  Int_t nEntries = reg_tree->GetEntries();
  Int_t last_run = -1;
  const int fRingSize=5;
  vector< vector<Double_t> > fEVRing; // [RingSize] x [nIV*nIV]
  vector<Double_t> fSlugAvg = fRingAvg;
  // cout << " -- Skeleton "; PrintVector(fSlugAvg);

  for(int ievt=0;ievt<nEntries;ievt++){
    reg_tree->GetEntry(ievt);
    vector<Int_t> fMap = CheckIdentity(fRingAvg, fEigenVectorRaw); // Checking with fRingAvg(skeleton)
    vector<Double_t> fReMapped = RemapVectors(fEigenVectorRaw, fMap);
    vector<Int_t> fSign = CheckSignWithRingAvg( fRingAvg, fReMapped);
    vector<Double_t> fFlipped = FlipVectors(fReMapped, fSign); 
    fEVRing.push_back(fFlipped);
  }
  fSlugAvg = GetRingAverage(fEVRing);
  fEVRing.clear();

  for(int ievt=0;ievt<nEntries;ievt++){
    reg_tree->GetEntry(ievt);
    vector<Int_t> fMap = CheckIdentity(fSlugAvg, fEigenVectorRaw); // Checking with fSlugAvg
    vector<Double_t> fReMapped = RemapVectors(fEigenVectorRaw, fMap);
    vector<Int_t> fSign = CheckSignWithRingAvg( fSlugAvg, fReMapped);
    vector<Double_t> fFlipped = FlipVectors(fReMapped, fSign); 
    fEVRing.push_back(fFlipped);
  }
  fSlugAvg = GetRingAverage(fEVRing);
  // cout << " -- fSlugAvg "; PrintVector(fSlugAvg);
  if(slug==94){
    vector<Int_t> fslug94map(12);
    vector<Int_t> fslug94sign(12,1);
    for(int i=0;i<12;i++)
      fslug94map[i] =i;
    fslug94map[8] = 7;
    fslug94map[7] = 8;
    fSlugAvg = RemapVectors(fSlugAvg,fslug94map);
    fslug94sign[8] = -1;
    fSlugAvg = FlipVectors(fSlugAvg,fslug94sign);
  }
  fEVRing.clear();
  
  for(int ievt=0;ievt<nEntries;ievt++){
    reg_tree->GetEntry(ievt);
    
    if(fEvMonRaw[0].nsamp<4500) // short run cut
      continue;
    if(run!=last_run){
      last_run  = run;
      cout << " -----  " << endl;
      cout << " -- Run " << run << endl;
      cout << " -----  " << endl;
    }
    
    // Checking eigenvector identities 
    // re-allocation and sign lock-in to  fEigenVectorSorted
    // cout << "-- fEigenVectorRaw: " ;PrintVector(fEigenVectorRaw);

    vector<Int_t> fMap = CheckIdentityWithSlug(fSlugAvg, fEigenVectorRaw);
    vector<Double_t> fReMapped = RemapVectors(fEigenVectorRaw, fMap);
    vector<Int_t> fSign = CheckSignWithRingAvg( fSlugAvg, fReMapped);
    vector<Double_t> fFlipped = FlipVectors(fReMapped, fSign);
    
    for(int imon=0;imon<nIV;imon++){
      fEvMon[ fMap[imon] ].mean = fEvMonRaw[imon].mean * fSign[ fMap[imon] ];
      fEvMon[ fMap[imon] ].err = fEvMonRaw[imon].err;
      fEvMon[ fMap[imon] ].rms = fEvMonRaw[imon].rms;
      fEvMon[ fMap[imon] ].nsamp = fEvMonRaw[imon].nsamp;
    }
    
    for(int irow=0;irow<nIV;irow++)
      for(int icol=0;icol<nIV;icol++)
    	fEigenVectorSorted[irow*nIV+icol]  = fFlipped[irow*nIV+icol];
    
    for(int i=0;i<nDV;i++)
      for(int j=0;j<nIV;j++){
    	fEigenLagrangeSlopes[i*nIV+fMap[j]] = fSign[fMap[j]] * fLagrangeSlope[i*nIV + j ];
	fEigenRegressionSlopes[i*nIV+fMap[j]] = fSign[fMap[j]] *fRegressionSlope[i*nIV + j ];
      }

    if( fRunInfo.find(run)!= fRunInfo.end()){
      slug = (fRunInfo[run]).first;
      arm_flag = (fRunInfo[run]).second;
      kHelicitySign = fSpinInfo[run];
      eig_sort_tree->Fill();
      lagr_slope_tree->Fill();
      reg_slope_tree->Fill();
    }else{
      cout << " **  run info not found " << endl;
    }

  }
  output->cd();
  eig_sort_tree->Write();
  lagr_slope_tree->Write();
  reg_slope_tree->Write();
  output->Close();
}

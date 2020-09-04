/*
  author: Tao Ye <tao.ye@stonybrook.edu>
  last update: July 2020

 */
#include "utilities_eigen.cc"

void SortEigenvectorsBySlug(Int_t slug){
  map<Int_t , pair<Int_t, Int_t> > fRunInfo  = LoadRunInfo();
  TString filename_tag = "allbpm";
  vector< TString > IVlist =  {"bpm4aX","bpm4eX","bpm1X",
			       "bpm11X","bpm12X","bpm16X",
			       "bpm4aY","bpm4eY","bpm1Y",
			       "bpm11Y","bpm12Y","bpm16Y"};

  vector<TString> DVlist ={"us_avg","usl","usr","us_dd"};
  Int_t nDV = DVlist.size();
  Int_t nIV = IVlist.size();
  TFile* output = TFile::Open(Form("rootfiles/slug%d_sorted_eigenvector_%s.root",
				   slug,filename_tag.Data()),"RECREATE");
  TTree *lagr_slope_tree = new TTree("lagr_slope","lagrange slopes in eigenvectors basis");
  TTree *reg_slope_tree = new TTree("reg_slope","regression slopes in eigenvectors basis");
  TTree *eig_vec_tree = new TTree("eig","eigenvectors");
  
  vector<Double_t> fRegEigSlopes(nDV*nIV);
  vector<Double_t> fLagrEigSlopes(nDV*nIV);
  vector<Double_t> fSign(nIV);
  vector<Double_t> fLocker(nIV*nIV,0);
  vector<Double_t> fEigenVectorSorted(nIV*nIV);
  typedef struct {Double_t mean,err,rms,m2,nsamp;} STAT;
  vector<STAT> fEvMon(nIV);
  
  Int_t run, burst_counter, arm_flag;
  eig_vec_tree->Branch("run",&run);
  eig_vec_tree->Branch("mini",&burst_counter);
  eig_vec_tree->Branch("arm_flag",&arm_flag);
  for(int j=0;j<nIV;j++)
    eig_vec_tree->Branch(Form("diff_evMon%d",j),&fEvMon[j],"mean/D:err:rms:m2:nsamp");
  
  for(int i=0;i<nIV;i++)
    for(int j=0;j<nIV;j++)
      eig_vec_tree->Branch(Form("evMon%d_%s",j, IVlist[i].Data()),
			 &fEigenVectorSorted[i*nIV+j]); // be carefull with row x col

  for(int i=0;i<nDV;i++)
    for(int j=0;j<nIV;j++){
      lagr_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			&fLagrEigSlopes[i*nIV+j]);
      reg_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
		       &fRegEigSlopes[i*nIV+j]);
    }

  TFile *input = TFile::Open(Form("treeMergeOutput/MergedLagrange_slug%d.root",slug));
  TTree *mini_tree = (TTree*)input->Get("mini");
  TTree *eig_tree = (TTree*)input->Get("mini_eigall");
  TTree *lagr_tree = (TTree*)input->Get("mini_lagrall");
  TTree *reg_tree = (TTree*)input->Get("mini_regall");

  eig_tree->AddFriend(lagr_tree);
  eig_tree->AddFriend(mini_tree);
  eig_tree->AddFriend(reg_tree);
  
  vector<STAT> fEvMonRaw(nIV);
  vector<Double_t> fEigenVectorRaw(nIV*nIV);
  vector<Double_t> fLagrangeSlope(nDV*nIV);
  vector<Double_t> fRegressionSlope(nDV*nIV);  // FIXME
  mini_tree->SetBranchAddress("run",&run);
  mini_tree->SetBranchAddress("mini",&burst_counter);
  
  for(int i=0;i<nIV;i++)
    eig_tree->SetBranchAddress(Form("diff_evMon%d",i),
  			       &fEvMonRaw[i]);
  
  for(int i=0;i<nIV;i++)
    for(int j=0;j<nIV;j++)
      eig_tree->SetBranchAddress(Form("evMon%d_%s",j,IVlist[i].Data()),
				 &fEigenVectorRaw[i*nIV+j]);
  // [row]:bpm ; [col]:eigv
  for(int i=0;i<nDV;i++)
    for(int j=0;j<nIV;j++){
      lagr_tree->SetBranchAddress(DVlist[i]+"_"+IVlist[j],
				  &fLagrangeSlope[i*nIV+j]);
      reg_tree->SetBranchAddress(DVlist[i]+"_"+IVlist[j],
				 &fRegressionSlope[i*nIV+j]);
    }
  Int_t nEntries = eig_tree->GetEntries();
  Int_t last_run = -1;
  
  const Int_t nRingsize = 5;
  vector< vector<Double_t> > fTemplate; // [nIV] x [nIV]
  vector< vector<Double_t> > fEVRing; // [RingSize] x [nIV*nIV]
  vector< Double_t > fRingAvg(nIV*nIV,0);
  
  eig_tree->GetEntry(0);
  fEVRing.push_back(fEigenVectorRaw);
  fRingAvg = GetRingAverage(fEVRing);
  fTemplate = GetTemplate(fRingAvg);
  
  for(int ievt=0;ievt<nEntries;ievt++){
    eig_tree->GetEntry(ievt);
    if(run!=last_run){
      last_run  = run;
      cout << " -----  " << endl;
      cout << " -- Run " << run << endl;
      cout << " -----  " << endl;
    }
    // Load raw BPM slopes to matrix
    TMatrixD lagr_raw_slope_matrix(nIV,nDV);
    TMatrixD reg_raw_slope_matrix(nIV,nDV);
    for(int i=0;i<nDV;i++)
      for(int j=0;j<nIV;j++){
	lagr_raw_slope_matrix[j][i] = fLagrangeSlope[i*nIV+j];
	reg_raw_slope_matrix[j][i] = fRegressionSlope[i*nIV+j];
      }
    
    // Checking eigenvector identities 
    // re-allocation and sign lock-in to  fEigenVectorSorted

    vector<Int_t> fSign = CrossCheckWithTemplate( fTemplate, fEigenVectorRaw);
    vector<Int_t> fMap = CheckIdentity(fRingAvg, fEigenVectorRaw);
    vector<Double_t> fFlipped = FlipVectors(fEigenVectorRaw, fSign); // Flip before remap
    vector<Double_t> fReMapped = RemapVectors(fFlipped, fMap);

    vector<STAT> fFlippedMon(nIV);
    for(int imon=0;imon<nIV;imon++){
      fFlippedMon[imon].mean = fEvMonRaw[imon].mean * fSign[imon];
      fFlippedMon[imon].err = fEvMonRaw[imon].err;
      fFlippedMon[imon].m2 = fEvMonRaw[imon].m2;
      fFlippedMon[imon].rms = fEvMonRaw[imon].rms;
      fFlippedMon[imon].nsamp = fEvMonRaw[imon].nsamp;
    }

    for(int imon=0;imon<nIV;imon++){
      fEvMon[ fMap[imon] ].mean = fFlippedMon[imon].mean;
      fEvMon[ fMap[imon] ].err = fFlippedMon[imon].err;
      fEvMon[ fMap[imon] ].m2 = fFlippedMon[imon].m2;
      fEvMon[ fMap[imon] ].rms = fFlippedMon[imon].rms;
      fEvMon[ fMap[imon] ].nsamp = fFlippedMon[imon].nsamp;
    }
    
    for(int irow=0;irow<nIV;irow++)
      for(int icol=0;icol<nIV;icol++)
	fEigenVectorSorted[irow*nIV+icol]  = fReMapped[irow*nIV+icol];

    
    if(ievt!=0){
      fEVRing.push_back(fEigenVectorSorted);
      if(fEVRing.size()> nRingsize)
	fEVRing.erase( fEVRing.begin() );
      fRingAvg = GetRingAverage(fEVRing);
      fTemplate = GetTemplate(fRingAvg);
    }
    
    //  Remix Slopes in eigenvector 
    TMatrixD eig_vec_matrix(nIV,nIV);   // [row]:bpm ; [col]:eigv
    for(int i=0;i<nIV;i++)
      for(int j=0;j<nIV;j++)
    	eig_vec_matrix[i][j] = fEigenVectorSorted[i*nIV+j];
    
    TMatrixD inv_eig_vec(eig_vec_matrix);
    inv_eig_vec.Invert();
    TMatrixD reg_eig_slope_matrix(nIV,nDV);
    TMatrixD lagr_eig_slope_matrix(nIV,nDV);
    
    lagr_eig_slope_matrix = inv_eig_vec * lagr_raw_slope_matrix;
    reg_eig_slope_matrix = inv_eig_vec * reg_raw_slope_matrix;
    
    for(int i=0;i<nDV;i++)
      for(int j=0;j<nIV;j++){
    	fLagrEigSlopes[i*nIV+j] = lagr_eig_slope_matrix[j][i];
	fRegEigSlopes[i*nIV+j]  = reg_eig_slope_matrix[j][i];
      }
    
    if( fRunInfo.find(run)!= fRunInfo.end()){
      arm_flag = (fRunInfo[run]).second;
      
      lagr_slope_tree->Fill();
      reg_slope_tree->Fill();
      eig_vec_tree->Fill();
    }


  }
  input->Close();
  output->cd();
  lagr_slope_tree->Write();
  reg_slope_tree->Write();
  eig_vec_tree->Write();
  output->Close();
}

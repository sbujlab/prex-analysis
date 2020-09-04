/*
  author: Tao Ye <tao.ye@stonybrook.edu>
  last update: July 2020

 */
#include "utilities_eigen.cc"

void SortEigenvectors(int start, int end){
  map<Int_t , pair<Int_t, Int_t> > fRunInfo  = LoadRunInfo();
  TString filename_tag = "allbpm";
  vector< TString > IVlist =  {"bpm4aX","bpm4eX","bpm1X",
			       "bpm11X","bpm12X","bpm16X",
			       "bpm4aY","bpm4eY","bpm1Y",
			       "bpm11Y","bpm12Y","bpm16Y"};

  vector<TString> DVlist ={"us_avg","usl","usr","us_dd"};
  Int_t nDV = DVlist.size();
  Int_t nIV = IVlist.size();
  // TFile* output = TFile::Open(Form("rootfiles/slug%d-%d_sorted_eigenvector_%s.root",
  // 				   start,end,filename_tag.Data()),"RECREATE");
  TFile* output = TFile::Open(Form("rootfiles/sorted_eigenvector_%s.root",
  				   filename_tag.Data()),"RECREATE");

  TTree *lagr_slope_tree = new TTree("lagr","lagrange slopes in eigenvectors basis");
  TTree *reg_slope_tree = new TTree("reg","regression slopes in eigenvectors basis");
  TTree *eig_sort_tree = new TTree("eig", "sorted eigenvectors tree");
  vector<Double_t> fEigenLagrangeSlopes(nDV*nIV);
  vector<Double_t> fEigenRegressionSlopes(nDV*nIV);
  vector<Double_t> fSign(nIV);
  vector<Double_t> fLocker(nIV*nIV,0);
  vector<Double_t> fEigenVectorSorted(nIV*nIV);
  typedef struct {Double_t mean,err,rms,m2,nsamp;} STAT;
  vector<STAT> fEvMon(nIV);
  
  Int_t run, burst_counter, arm_flag, slug;
  eig_sort_tree->Branch("slug",&slug);
  eig_sort_tree->Branch("arm_flag",&arm_flag);
  eig_sort_tree->Branch("run",&run);
  eig_sort_tree->Branch("mini",&burst_counter);
  for(int j=0;j<nIV;j++)
    eig_sort_tree->Branch(Form("diff_evMon%d",j),&fEvMon[j],"mean/D:err:rms:m2:nsamp");
  
  for(int i=0;i<nIV;i++)
    for(int j=0;j<nIV;j++)
      eig_sort_tree->Branch(Form("evMon%d_%s",j, IVlist[i].Data()),
			 &fEigenVectorSorted[i*nIV+j]); // be carefull with row x col

  for(int i=0;i<nDV;i++)
    for(int j=0;j<nIV;j++){
      lagr_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			      &fEigenLagrangeSlopes[i*nIV+j]);
      reg_slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			      &fEigenRegressionSlopes[i*nIV+j]);

    }

  TChain *mini_tree = new TChain("mini");
  TChain *eig_tree = new TChain("mini_eigall");
  TChain *reg_tree = new TChain("mini_regall");
  TChain *lagr_tree = new TChain("mini_lagrall");
  for(int i=start;i<=end;i++){
    mini_tree->Add(Form("treeMergeOutput/MergedLagrange_slug%d.root",i));
    eig_tree->Add(Form("treeMergeOutput/MergedLagrange_slug%d.root",i));
    lagr_tree->Add(Form("treeMergeOutput/MergedLagrange_slug%d.root",i));
    reg_tree->Add(Form("treeMergeOutput/MergedLagrange_slug%d.root",i));
  }
  eig_tree->AddFriend(lagr_tree);
  eig_tree->AddFriend(reg_tree);
  eig_tree->AddFriend(mini_tree);

  vector<STAT> fEvMonRaw(nIV);
  vector<Double_t> fEigenVectorRaw(nIV*nIV);
  vector<Double_t> fLagrangeSlope(nDV*nIV);
  vector<Double_t> fRegressionSlope(nDV*nIV);  
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
  
  vector< vector<Double_t> > fEVRing; // [RingSize] x [nIV*nIV]
  vector< Double_t > fRingAvg(nIV*nIV,0);
  
  eig_tree->GetEntry(0);
  fEVRing.push_back(fEigenVectorRaw);
  fRingAvg = GetRingAverage(fEVRing);

  for(int ievt=0;ievt<nEntries;ievt++){
    eig_tree->GetEntry(ievt);
    // lagr_tree->GetEntry(ievt);
    if(run!=last_run){
      last_run  = run;
      cout << " -----  " << endl;
      cout << " -- Run " << run << endl;
      cout << " -----  " << endl;
    }
    // Load raw BPM slopes to matrix
    TMatrixD raw_lagr_slope_matrix(nIV,nDV);
    TMatrixD raw_reg_slope_matrix(nIV,nDV);
    for(int i=0;i<nDV;i++)
      for(int j=0;j<nIV;j++){
	raw_lagr_slope_matrix[j][i] = fLagrangeSlope[i*nIV+j];
	raw_reg_slope_matrix[j][i] = fRegressionSlope[i*nIV+j];
      }
    
    // Checking eigenvector identities 
    // re-allocation and sign lock-in to  fEigenVectorSorted
    
    vector<Int_t> fMap = CheckIdentity(fRingAvg, fEigenVectorRaw);
    vector<Double_t> fReMapped = RemapVectors(fEigenVectorRaw, fMap);
    vector<Int_t> fSign = CheckSignWithRingAvg( fRingAvg, fReMapped);
    vector<Double_t> fFlipped = FlipVectors(fReMapped, fSign); 

    for(int imon=0;imon<nIV;imon++){
      fEvMon[ fMap[imon] ].mean = fEvMonRaw[imon].mean * fSign[ fMap[imon] ];
      fEvMon[ fMap[imon] ].err = fEvMonRaw[imon].err;
      fEvMon[ fMap[imon] ].m2 = fEvMonRaw[imon].m2;
      fEvMon[ fMap[imon] ].rms = fEvMonRaw[imon].rms;
      fEvMon[ fMap[imon] ].nsamp = fEvMonRaw[imon].nsamp;
    }
    
    for(int irow=0;irow<nIV;irow++)
      for(int icol=0;icol<nIV;icol++)
	fEigenVectorSorted[irow*nIV+icol]  = fFlipped[irow*nIV+icol];

    // Update Template and Ring averages after flipped and remapped. 
    if(ievt!=0){
      fEVRing.push_back(fEigenVectorSorted);
      fRingAvg = GetRingAverage(fEVRing);
    }
    
    //  Remix Slopes in eigenvector 
    TMatrixD eig_vec_matrix(nIV,nIV);   // [row]:bpm ; [col]:eigv
    for(int i=0;i<nIV;i++)
      for(int j=0;j<nIV;j++)
    	eig_vec_matrix[i][j] = fEigenVectorSorted[i*nIV+j];
    
    TMatrixD inv_eig_vec(eig_vec_matrix);
    inv_eig_vec.Invert();
    TMatrixD eig_lagr_slope_matrix(nIV,nDV);
    eig_lagr_slope_matrix = inv_eig_vec * raw_lagr_slope_matrix;

    TMatrixD eig_reg_slope_matrix(nIV,nDV);
    eig_reg_slope_matrix = inv_eig_vec * raw_reg_slope_matrix;

    for(int i=0;i<nDV;i++)
      for(int j=0;j<nIV;j++){
    	fEigenLagrangeSlopes[i*nIV+j] = eig_lagr_slope_matrix[j][i];
	fEigenRegressionSlopes[i*nIV+j] = eig_reg_slope_matrix[j][i];
      }
    
    if( fRunInfo.find(run)!= fRunInfo.end()){
      slug = (fRunInfo[run]).first;
      arm_flag = (fRunInfo[run]).second;
      eig_sort_tree->Fill();
      lagr_slope_tree->Fill();
      reg_slope_tree->Fill();
    }

  }
  vector<Double_t> eigen_width(nIV);
  vector<Int_t> fGlobalRank(nIV);
  for(int iev=0;iev<nIV;iev++){
    eig_sort_tree->Draw(Form("diff_evMon%d.rms*1e3",iev),"","goff");
    TH1D* hptr = (TH1D*)gDirectory->FindObject("htemp");
    eigen_width[iev] = hptr->GetMean();
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
  lagr_slope_tree->Write();
  reg_slope_tree->Write();
  output->Close();

  TFile *output_avg = TFile::Open("skeleton_matrix.root","RECREATE");
  output_avg->cd();
  output_avg->WriteObject( &fRingAvg, "fRingAvg");
  output_avg->WriteObject( &skeleton_m, "skeleton");
  output_avg->Close();
}

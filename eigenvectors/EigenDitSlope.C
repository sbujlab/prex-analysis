void EigenDitSlope(Int_t slug=94){
  
  vector< TString > IVlist1 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm8X","bpm12X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm8Y","bpm12Y"};
  vector< TString > IVlist2 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm11X","bpm12X","bpm16X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm11Y","bpm12Y","bpm16Y"};
  vector<TString> IVlist;
  if(slug>=3)
    IVlist = IVlist2;
  else
    IVlist = IVlist1;

  vector<TString> DVlist ={"us_avg","usl","usr"};
  Int_t nDV = DVlist.size();
  Int_t nIV = IVlist.size();
  TFile* output = TFile::Open(Form("rootfiles/dit_eigslopes_slug%d.root",slug),"RECREATE");
  TTree *slope_tree = new TTree("slope","lagrange slopes in eigenvectors basis");
  vector<Double_t> fEigenSlopes(nDV*nIV);
  vector<Double_t> fSign(nIV);
  vector<Double_t> fSignLocker(nIV*nIV);
  for(int j=0;j<nIV;j++)
    slope_tree->Branch(Form("sign%d",j),&fSign[j]);
    
  for(int i=0;i<nDV;i++)
    for(int j=0;j<nIV;j++)
      slope_tree->Branch(DVlist[i]+Form("_evMon%d",j),
			 &fEigenSlopes[i*nIV+j]);
  
  TFile* input = TFile::Open(Form("rootfiles/MergedLagrange_slug%d.root",slug));
  TTree *eig_tree = (TTree*)input->Get("mini_eig");
  TTree *lagr_tree = (TTree*)input->Get("mini_lagr");
  eig_tree->AddFriend(lagr_tree);
  vector<Double_t> fEigenVector(nIV*nIV);
  vector<Double_t> fLagrangeSlope(nDV*nIV);
    
  for(int i=0;i<nIV;i++)
    for(int j=0;j<nIV;j++)
      eig_tree->SetBranchAddress(Form("evMon%d_%s",j,IVlist[i].Data()),
				 &fEigenVector[i*nIV+j]);

  for(int i=0;i<nDV;i++)
    for(int j=0;j<nIV;j++)
      lagr_tree->SetBranchAddress(DVlist[i]+"_"+IVlist[j],
				  &fLagrangeSlope[i*nIV+j]);
  
  Int_t nEntries = eig_tree->GetEntries();
  for(int ievt=0;ievt<nEntries;ievt++){
    eig_tree->GetEntry(ievt);

    if(ievt==0){
      for(int i=0;i<nIV;i++)
	for(int j=0;j<nIV;j++)
	 fSignLocker[i*nIV+j] = fEigenVector[i*nIV+j];
    }
    
    TMatrixD raw_slope_matrix(nIV,nDV);
    for(int i=0;i<nDV;i++)
      for(int j=0;j<nIV;j++)
	raw_slope_matrix[j][i] = fLagrangeSlope[i*nIV+j];
    
    for(int icol=0;icol<nIV;icol++){
      Double_t prod_sum=0;
      for(int jrow=0;jrow<nIV;jrow++){
	prod_sum += (fSignLocker[jrow*nIV+icol] * fEigenVector[jrow*nIV+icol]);
      }
      if(prod_sum>0)
	fSign[icol] =1.0;
      else
	fSign[icol] =-1.0;
    }
    
    TMatrixD eig_vec_matrix(nIV,nIV);
    for(int i=0;i<nIV;i++)
      for(int j=0;j<nIV;j++)
	eig_vec_matrix[i][j] = fEigenVector[i*nIV+j];

    
    TMatrixD inv_eig_vec(eig_vec_matrix);
    inv_eig_vec.Invert();
    TMatrixD eig_slope_matrix(nIV,nDV);
    eig_slope_matrix = inv_eig_vec * raw_slope_matrix;

    for(int i=0;i<nDV;i++)
      for(int j=0;j<nIV;j++)
	fEigenSlopes[i*nIV+j] = eig_slope_matrix[j][i];
    slope_tree->Fill();
  }


  input->Close();
  output->cd();
  slope_tree->Write();
  output->Close();
}

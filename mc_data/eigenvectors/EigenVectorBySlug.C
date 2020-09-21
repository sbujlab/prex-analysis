void EigenVectorBySlug(Int_t slug=94){

  vector< TString > IVlist = {"diff_bpm4aX","diff_bpm4eX","diff_bpm1X",
			      "diff_bpm11X","diff_bpm12X","diff_bpm16X",
			      "diff_bpm4aY","diff_bpm4eY","diff_bpm1Y",
			      "diff_bpm11Y","diff_bpm12Y","diff_bpm16Y"};
  TFile* output = TFile::Open(Form("./rootfiles/slug%d.root",slug),"RECREATE");
  Int_t nBPM = IVlist.size();;
  vector<Int_t> fIndex(nBPM,-1);
  vector<Int_t> fSignLockID(nBPM,0);
  TTree *eig_tree = new TTree("eig","Eigenvector Monitors");
  Double_t *fProto = new Double_t[nBPM];
  Int_t fVectorID;
  Int_t fMyStat;
  Int_t fMini=0;
  Int_t fRun;
  Double_t fEigenValue;
  Int_t fNdim = nBPM;
  Int_t fSign =1;
  eig_tree->Branch("MyStat",&fMyStat);
  eig_tree->Branch("eigID",&fVectorID);
  eig_tree->Branch("run",&fRun);
  eig_tree->Branch("minirun",&fMini);
  eig_tree->Branch("ndim",&fNdim);
  eig_tree->Branch("eigval",&fEigenValue);
  eig_tree->Branch("sign",&fSign);
  TString leaflist="";
  for(int i=0;i<nBPM;i++){
    TString tag = IVlist[i];
    tag.ReplaceAll("diff_bpm","");
    leaflist += tag+"/D";
    if(i!=nBPM-1)
      leaflist +=":";
  }
  TBranch *fBranch = eig_tree->Branch("eigvec",0,leaflist);
  for(int i=0;i<nBPM;i++){
    TString tag = IVlist[i];
    tag.ReplaceAll("diff_bpm","");
    fBranch->GetLeaf(tag)->SetAddress(&fProto[i]);
  }

  TString listname = Form("prex-runlist/simple_list/slug%d.list",94);
  FILE *runlist = fopen(listname.Data(),"r");
  Bool_t kFirstRun = kTRUE; // used for vector sign lock
  while(!feof(runlist)){
    fRun=0;
    fscanf(runlist,"%d\n",&fRun);
    if(fRun==0)
      break;
    
    TString lrb_filename = Form("./LRBoutput/burst_blueR%d.000all.slope.root",fRun);
    TFile *lrb_file = TFile::Open(lrb_filename);
    if(lrb_file==NULL)
      continue;
    cout << fRun << endl;

    if((TMatrixD*)lrb_file->Get("MyStat")==NULL)
      continue;
    TH1D* hist_iv = (TH1D*)lrb_file->Get("IVname");
    TH1D* hist_dv = (TH1D*)lrb_file->Get("DVname");
    TAxis *ivAxis = hist_iv->GetXaxis();
    TAxis *dvAxis = hist_dv->GetXaxis();
    Int_t nIV = ivAxis->GetLast();
    Int_t nDV = dvAxis->GetLast();


    for(int i=0;i<nIV;i++){
      const char* char_buff = ivAxis->GetBinLabel(i+1);
      auto itf= find( IVlist.begin(),IVlist.end(),TString(char_buff));
      if(itf!=IVlist.end()){
	fIndex[ itf - IVlist.begin() ] = i;
      }
    }

    for(int i=0;i<nBPM;i++){
      cout << i << ":" << IVlist[i] << endl;
    }
    
    fMini=0;
    Int_t cycle=1;
    TString cycle_tag = Form(";%d",cycle);
    while( ((TMatrixD*)lrb_file->Get("MyStat"+cycle_tag)!=NULL)){

      TMatrixT<double> MyStat = *((TMatrixT<double>*)lrb_file->Get("MyStat"+cycle_tag));
      TMatrixT<double> IV_IV_normVar = *((TMatrixT<double>*)lrb_file->Get("IV_IV_normVariance"+cycle_tag));
      fMyStat = MyStat[0][0];
    
      TMatrixDSym S_IV(nBPM);
      for(int i=0;i<nBPM;i++)
	for(int j=0;j<nBPM;j++)
	  S_IV[i][j] = IV_IV_normVar[fIndex[i]][fIndex[j]];

      TMatrixDSymEigen S_IV_eig(S_IV);
      TMatrixD eigen_vector = S_IV_eig.GetEigenVectors();
      TVectorD eigen_values = S_IV_eig.GetEigenValues();
      // eigen_vector.Print();
      // eigen_values.Print();
      // TMatrixD eigen_vector_trans(eigen_vector);
      // eigen_vector_trans.T();
      // (eigen_vector*eigen_vector_trans).Print();

      if(kFirstRun){
	for(int i=0;i<nBPM;i++){
	  Double_t prim_comp = fabs(eigen_vector[0][i]);
	  Int_t prim_index = 0;
	  for(int j=1;j<nBPM;j++){
	    if(fabs(eigen_vector[j][i])>prim_comp){
	      prim_comp = fabs(eigen_vector[j][i]);
	      prim_index = j;
	    }
	  } // component loop
	  fSignLockID[i] = prim_index;
	} // eigenvector loop
	kFirstRun = kFALSE;
      }
      
      for(int i=0;i<nBPM;i++){
	fVectorID = i;
	fEigenValue = eigen_values[i];

	if( eigen_vector[fSignLockID[i]][i]<0)
	  fSign = -1;
	else
	  fSign = 1;
	
	for(int j=0;j<nBPM;j++)
	  fProto[j] = fSign*eigen_vector[j][i];
	
	eig_tree->Fill();      
      }

      fMini++;    
      cycle++;
      cycle_tag = Form(";%d",cycle);
    
    }
    lrb_file->Close();
    
  } // end of run loop

  output->cd();
  eig_tree->Write();
  output->Close();
}

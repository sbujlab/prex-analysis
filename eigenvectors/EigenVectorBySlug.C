void EigenVectorBySlug(Int_t slug=94){

  vector< TString > IVlist = {"diff_bpm4aX","diff_bpm4eX","diff_bpm1X",
  			      "diff_bpm11X","diff_bpm12X","diff_bpm16X",
  			      "diff_bpm4aY","diff_bpm4eY","diff_bpm1Y",
  			      "diff_bpm11Y","diff_bpm12Y","diff_bpm16Y"};
  // vector< TString > IVlist = {"diff_bpm4aX","diff_bpm4eX",
  // 			      "diff_bpm11X",
  // 			      "diff_bpm4aY","diff_bpm4eY"};
  
  vector<TString> DVlist={"asym_usl","asym_usr","asym_us_avg","asym_us_dd"};
  TFile* output = TFile::Open(Form("./rootfiles/slug%d.root",slug),"RECREATE");
  Int_t nBPM = IVlist.size();
  vector<Int_t> fIVIndex(nBPM,-1);
  Int_t nDet = DVlist.size();
  vector<Int_t> fDVIndex(nDet,-1);

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
  // Contruct Eigen-Tree
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
  // Contruct Slope Tree
  TTree *slope_tree = new TTree("slope","eigenvector regression slope");
  slope_tree->Branch("MyStat",&fMyStat);
  slope_tree->Branch("run",&fRun);
  slope_tree->Branch("minirun",&fMini);
  vector<Double_t> fdummy(nBPM);
  vector< vector <Double_t> > fSlope(nDet,fdummy);
  for(int idet=0;idet<nDet;idet++){
    TString det_tag = DVlist[idet];
    det_tag.ReplaceAll("asym_","");
    for(int ibpm=0;ibpm<nBPM;ibpm++){
      slope_tree->Branch(Form("%s_evMon%d",det_tag.Data(),ibpm), &fSlope[idet][ibpm]);
    }
  }
    
  TString listname = Form("prex-runlist/simple_list/slug%d.list",slug);
  FILE *runlist = fopen(listname.Data(),"r");
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
	fIVIndex[ itf - IVlist.begin() ] = i;
      }
    }

    for(int i=0;i<nIV;i++){
      const char* char_buff = dvAxis->GetBinLabel(i+1);
      auto itf= find( DVlist.begin(),DVlist.end(),TString(char_buff));
      if(itf!=DVlist.end()){
	fDVIndex[ itf - DVlist.begin() ] = i;
      }
    }

    for(int i=0;i<nBPM;i++){
      cout << i << ":" << IVlist[i] << "@" << fIVIndex[i] << endl;
    }
    for(int i=0;i<nDet;i++){
      cout << i << ":" << DVlist[i] << "@" << fDVIndex[i] << endl;
    }
    
    fMini=0;
    Int_t cycle=1;
    TString cycle_tag = Form(";%d",cycle);
    while( ((TMatrixD*)lrb_file->Get("MyStat"+cycle_tag)!=NULL)){

      TMatrixT<double> MyStat = *((TMatrixT<double>*)lrb_file->Get("MyStat"+cycle_tag));
      TMatrixT<double> IV_IV_normVar = *((TMatrixT<double>*)lrb_file->Get("IV_IV_normVariance"+cycle_tag));
      TMatrixT<double> IV_DV_normVar = *((TMatrixT<double>*)lrb_file->Get("IV_DV_normVariance"+cycle_tag));
      fMyStat = MyStat[0][0];
    
      TMatrixDSym S_IV(nBPM);
      TMatrixD sub_IV_DV_normVar(nBPM,nDet);
      for(int i=0;i<nBPM;i++)
	for(int j=0;j<nBPM;j++)
	  S_IV[i][j] = IV_IV_normVar[fIVIndex[i]][fIVIndex[j]];

      for(int idet=0;idet<nDet;idet++)
	for(int ibpm=0;ibpm<nBPM;ibpm++)
	  sub_IV_DV_normVar[ibpm][idet] = IV_DV_normVar[fIVIndex[ibpm]][fDVIndex[idet]];

      TMatrixDSymEigen S_IV_eig(S_IV);
      TMatrixD eigen_vector = S_IV_eig.GetEigenVectors();
      TVectorD eigen_values = S_IV_eig.GetEigenValues();
      TMatrixD lambda_matrix(nBPM,nBPM);
      for(int i=0;i<nBPM;i++)
	for(int j=0;j<nBPM;j++)
	  if(i==j)
	    lambda_matrix[i][j]=eigen_values[i];
	  else
	    lambda_matrix[i][j]=0.0;
      
      for(int i=0;i<nBPM;i++){
	fVectorID = i;
	fEigenValue = eigen_values[i];
	Double_t test_val =0;
	
	if( fabs(eigen_vector[4][i]) > fabs(eigen_vector[10][i]) )
	  test_val = eigen_vector[4][i];
	else
	  test_val = eigen_vector[10][i];
	if(test_val>0)
	  fSign=1.0;
	else
	  fSign=-1.0;
	for(int j=0;j<nBPM;j++){
	  eigen_vector[j][i] = fSign*eigen_vector[j][i];
	  fProto[j] = eigen_vector[j][i];
	}
	
	eig_tree->Fill();      
      }

      TMatrixD eigen_vector_trans (eigen_vector);
      eigen_vector_trans.T();
      TMatrixD eigenIV_DV = eigen_vector_trans*sub_IV_DV_normVar;
      TMatrixD inv_lambda(lambda_matrix);
      inv_lambda.Invert();
      TMatrixD slope_matrix(nBPM,nDet);
      slope_matrix =  inv_lambda * eigenIV_DV;
      for(int idet=0;idet<nDet;idet++)
	for(int ibpm=0;ibpm<nBPM;ibpm++)
	  fSlope[idet][ibpm] = slope_matrix[ibpm][idet];
      
      slope_matrix.Print();

      slope_tree->Fill();
      fMini++;    
      cycle++;
      cycle_tag = Form(";%d",cycle);
    
    }
    lrb_file->Close();
    
  } // end of run loop

  output->cd();
  slope_tree->Write();
  eig_tree->Write();
  output->Close();
}

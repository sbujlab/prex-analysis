void EigenVectorByRun(Int_t fRun=4980){
  TString lrb_filename = Form("./LRBoutput/blueR%d.000all.slope.root",fRun);
  TFile *lrb_file = TFile::Open(lrb_filename);
  if(lrb_file==NULL)
    return;
  cout << fRun << endl;

  if((TMatrixD*)lrb_file->Get("MyStat")==NULL)
    return;
  TH1D* hist_iv = (TH1D*)lrb_file->Get("IVname");
  TH1D* hist_dv = (TH1D*)lrb_file->Get("DVname");
  TAxis *ivAxis = hist_iv->GetXaxis();
  TAxis *dvAxis = hist_dv->GetXaxis();
  Int_t nIV = ivAxis->GetLast();
  Int_t nDV = dvAxis->GetLast();

  // vector< TString > IVlist = {"diff_bpm4aX","diff_bpm4eX",
  // 			      "diff_bpm4aY","diff_bpm4eY",
  // 			      "diff_bpm11X","diff_bpm12X"};
  vector< TString > IVlist = {"diff_bpm4aX","diff_bpm4eX","diff_bpm1X",
  			      "diff_bpm11X","diff_bpm12X","diff_bpm16X",
  			      "diff_bpm4aY","diff_bpm4eY","diff_bpm1Y",
  			      "diff_bpm11Y","diff_bpm12Y","diff_bpm16Y"};

  if(IVlist.size()==0){
    for(int i=0;i<nIV;i++){
      const char* char_buff = ivAxis->GetBinLabel(i+1);
      IVlist.push_back( TString(char_buff));
    }
  }

  vector< TString > DVlist;
  Int_t nBPM = IVlist.size();
  vector<Int_t> fIndex(nBPM,-1);
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
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  gStyle->SetPalette(kBeach);
  TColor::InvertPalette();
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->SetFillStyle(4000);
  c1->Divide(2,1);
  // c1->Print(Form("run%d_diagonalized_matrices.pdf[",fRun));
  TFile* output = TFile::Open(Form("./rootfiles/prex_eigenvec_%d.root",fRun),"RECREATE");
  TTree *eig_tree = new TTree("eig","Eigenvector Monitors");
  Double_t *fProto = new Double_t[nBPM];
  Int_t fVectorID;
  Int_t fMyStat;
  Int_t fMini=0;
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
  // for(int i=0;i<nDV;i++){
  //   const char* char_buff = dvAxis->GetBinLabel(i+1);
  //   DVlist.push_back(TString(char_buff));
  // }

  ////////////////////////////////////                                                                   

  Int_t cycle=1;
  TString cycle_tag = Form(";%d",cycle);
  while( ((TMatrixD*)lrb_file->Get("MyStat"+cycle_tag)!=NULL)){

    TMatrixT<double> MyStat = *((TMatrixT<double>*)lrb_file->Get("MyStat"+cycle_tag));
    TMatrixT<double> IV_IV_normVar = *((TMatrixT<double>*)lrb_file->Get("IV_IV_normVariance"+cycle_tag))\
;
    fMyStat = MyStat[0][0];
    
    TMatrixDSym S_IV(nBPM);
    for(int i=0;i<nBPM;i++)
      for(int j=0;j<nBPM;j++)
	S_IV[i][j] = IV_IV_normVar[fIndex[i]][fIndex[j]];

    TMatrixDSymEigen S_IV_eig(S_IV);
    TMatrixD eigen_vector = S_IV_eig.GetEigenVectors();
    TVectorD eigen_values = S_IV_eig.GetEigenValues();
    eigen_vector.Print();
    eigen_values.Print();
    TMatrixD eigen_vector_trans(eigen_vector);
    eigen_vector_trans.T();
    (eigen_vector*eigen_vector_trans).Print();

    for(int i=0;i<nBPM;i++){
      fVectorID = i;
      fEigenValue = eigen_values[i];

	
      for(int j=0;j<nBPM;j++)
	fProto[j] = eigen_vector[j][i];
	
      eig_tree->Fill();      
    }

    fMini++;    
    cycle++;
    cycle_tag = Form(";%d",cycle);

    c1->cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetFillStyle(4000);
    TH2D* raw_cov = new TH2D("raw_cov","Covariance (um^{2})",
			     nBPM,-0.5,nBPM-0.5,
			     nBPM,-0.5,nBPM-0.5); 
    for(int i=0;i<nBPM;i++)
      for(int j=0;j<nBPM;j++)
	raw_cov->Fill(j,nBPM-1-i, S_IV[i][j]*1e6); // this is really a tricky
    // raw_cov->GetXaxis()->SetNdivisions(nBPM);
    for(int i=0;i<nBPM;i++){
      TString label = IVlist[i];
      label.ReplaceAll("diff_bpm","");
      raw_cov->GetXaxis()->SetBinLabel(i+1,label);
      raw_cov->GetYaxis()->SetBinLabel(nBPM-i,label);
    }
    raw_cov->GetXaxis()->SetTickLength(0.01);
    raw_cov->GetYaxis()->SetTickLength(0.01);
    raw_cov->Draw("COLZ TEXT");

    c1->cd(2);
    gPad->SetRightMargin(0.15);
    gPad->SetFillStyle(4000);
    TH2D* diag_cov = new TH2D("diag_cov","Diagonalized (um^{2} )",
			     nBPM,-0.5,nBPM-0.5,
			     nBPM,-0.5,nBPM-0.5); 

    TMatrixD diagS_IV(S_IV);
    diagS_IV  = eigen_vector_trans*S_IV*eigen_vector;
    for(int i=0;i<nBPM;i++)
      for(int j=0;j<nBPM;j++){
	if(fabs(diagS_IV[i][j]*1e6)>0.001)
	  diag_cov->Fill(j,nBPM-1-i, diagS_IV[i][j]*1e6);
	else
	  diag_cov->Fill(j,nBPM-1-i, 0);
      }

    for(int i=0;i<nBPM;i++){
      diag_cov->GetXaxis()->SetBinLabel(i+1,Form("%d",i));
      diag_cov->GetYaxis()->SetBinLabel(nBPM-i,Form("%d",i));
    }
    diag_cov->GetXaxis()->SetTickLength(0.01);
    diag_cov->GetYaxis()->SetTickLength(0.01);
    gPad->SetLogz();
    diag_cov->Draw("COLZ TEXT");
    
    
    c1->Print(Form("run%d_diagonalized_matrices.pdf",fRun));
  }
  // c1->Print(Form("run%d_diagonalized_matrices.pdf]",fRun));
  lrb_file->Close();
  output->cd();
  eig_tree->Write();
  output->Close();
}

void calculator(){
  TFile *input = TFile::Open("prex_sens_matrix.4066-4080.root");
  TMatrixD *sens_matrix = (TMatrixD *) input->Get("sens_matrix");
  Int_t coil_idx[7] = {1,3,4,6,7,2,5};
  TString bpm_name[12] = { "4aX","4aY","4eX","4eY",
			   "1X","1Y",
			   "11X","11Y","12X","12Y",
			   "16X","16Y"};

  Int_t bpm_idx[12] = {0,1,2,3,
		       11,12,
		       9,10,7,8,
		       5,6  };

  TString det_name[2]={"usl","usr"};
  Int_t det_idx[2]={13,14};

  TMatrixD bpm_coil(12,7);
  for(int irow=0;irow<12;irow++){
    for(int icol=0;icol<7;icol++){
      bpm_coil[irow][coil_idx[icol]-1] = (*sens_matrix)[bpm_idx[irow]][icol] * 1e3; // [um/counts]
    }
  }

  TMatrixD det_coil(2,7);
  for(int irow=0;irow<2;irow++){
    for(int icol=0;icol<7;icol++){
      det_coil[irow][coil_idx[icol]-1] = (*sens_matrix)[det_idx[irow]][icol]*1e6;  // [ppm/counts]
    }
  }

  // beam parameters : posX, angleX, posY, angleY, energy
  // unit: [um/100counts] 
  double par_coil[5][7]={{-2.6, 0,-80.2,0,-24.3, 0 , 1.9},
			 {-162.25, 0 ,-51.18, 0 ,188.05, 0, 20.63},
			 { 0 ,-117.77, 0,-65.35, 0 ,132.44,-25.04},
			 {0 ,13.43,0,-47.14,0,-11.40,0},
			 {0, 0, 0, 0, 0, 0,-124.87}  };

  TMatrixD dAlpha_dcoil(5,7);
  for(int irow=0;irow<5;irow++){
    printf("{");
    for(int icol=0;icol<7;icol++){
      dAlpha_dcoil[irow][icol]  = par_coil[irow][icol]/ 100; //convert to [um/counts]
      printf("%.2g",dAlpha_dcoil[irow][icol]);
      if(icol!=6)
	printf(",");
    }
    printf("}");
    if(irow!=4)
      printf(",");
    printf("\n");
  }
  TMatrixD dAlpha_dcoil_trans(dAlpha_dcoil);
  dAlpha_dcoil_trans.T();
  TMatrixD CCT(5,5);
  CCT = dAlpha_dcoil * dAlpha_dcoil_trans;
  TMatrixD inv_CCT(CCT);
  inv_CCT.Invert();
  cout << "dBPM_dAlpha " << endl;
  TMatrixD dBPM_dAlpha(12,5); // bpm_coil = dBPM_dAlpha * dAlpha_dcoil
  dBPM_dAlpha = bpm_coil*dAlpha_dcoil_trans*inv_CCT;
  for(int i=0;i<12;i++){
    cout << "{";
    for(int j=0;j<5;j++){
      printf("%.2g",dBPM_dAlpha[i][j]);
      if(j!=4)
	printf(",");
    }
    cout << "}" ;
    if(i!=11)
      cout << ",";
    cout << endl;
  }
      
  cout << "dDet_dAlpha " << endl;
  TMatrixD dDet_dAlpha(2,5);
  dDet_dAlpha = det_coil*dAlpha_dcoil_trans*inv_CCT;
  dDet_dAlpha.Print();
  
  cout << "{";
  for(int j=0;j<5;j++){
    printf("%.2g",dDet_dAlpha[0][j] );
    if(j!=4)
      printf(",");
  }
  cout << "}" << endl;;

  cout << "{";
  for(int j=0;j<5;j++){
    printf("%.2g",dDet_dAlpha[1][j] );
    if(j!=4)
      printf(",");
  }
  cout << "}" << endl;
  
  cout << "{";
  for(int j=0;j<5;j++){
    printf("%.2g",0.5*dDet_dAlpha[0][j] +0.5*dDet_dAlpha[1][j] );
    if(j!=4)
      printf(",");
  }
  cout << "}" << endl;;

  cout << "beta" << endl;
  TRandom3  *frng = new TRandom3();
  for(int i=0;i<12;i++){
    cout << "{" ;
    for(int j=0;j<12;j++){
      if( i==j)
	printf("0.0");
      else
	printf("%.2f",frng->Uniform(-1,1));
      if(j!=11)
	printf(",");
    }
    printf("}");
    if(i!=11)
      printf(",");
    printf("\n");
    
  }
}

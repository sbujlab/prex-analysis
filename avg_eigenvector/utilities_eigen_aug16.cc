/*
  author: Tao Ye <tao.ye@stonybrook.edu>
  last update: July 2020

 */
#include "TaAccumulator.cc"

void PrintMatrixD(TMatrixD input);
void PrintVector(vector<Double_t> input);
void PrintVector(vector<Int_t> input);

vector<Int_t> CheckIdentity(vector<Double_t> fLocker,vector<Double_t> fRaw);
vector<Int_t> CheckSignWithRingAvg( vector<Double_t>  fRingAvg,
				     vector<Double_t> fRaw );
vector<Double_t> RemapVectors(vector<Double_t> fRaw, vector<Int_t> fmap);
vector<Double_t> FlipVectors(vector<Double_t> fRaw, vector<Int_t> fSign);
map<Int_t, pair<Int_t,Int_t> > LoadRunInfo();
vector<Double_t> GetRingAverage( vector< vector<Double_t> > fRing) ;

vector<Int_t> CheckIdentity(vector<Double_t> fLocker,
			    vector<Double_t> fRaw){

  Int_t nLength = fRaw.size();
  Int_t nIV  = (Int_t)sqrt((Double_t)nLength);
  vector<Int_t> fMap(nIV);
  vector<Double_t> fCorrelation_max(nIV);
  vector<Double_t> fCorrelation(nIV*nIV,0.0);
  vector<Double_t> fDistance(nIV*nIV,0.0);

  vector<Int_t> fLockerRanks(nIV*nIV);
  vector<Double_t> fAbsValue(nIV*nIV);
  for(int irow=0;irow<nIV;irow++)
    for(int icol=0;icol<nIV;icol++){
      fAbsValue[irow*nIV+icol] = fabs(fLocker[irow*nIV+icol]);
      fLockerRanks[irow*nIV +icol] = irow;
    }

  for(int icol=0; icol<nIV;icol++){
    for(int irow=0;irow<nIV;irow++){
      double max_val = fAbsValue[irow*nIV+icol];
      for(int jrow=irow+1;jrow<nIV;jrow++){
	if( fAbsValue[jrow*nIV+icol] > max_val) {
	  max_val = fAbsValue[jrow*nIV+icol];
	  fAbsValue[jrow*nIV+icol] = fAbsValue[irow*nIV+icol];
	  fAbsValue[irow*nIV+icol] = max_val;
	  int swap_row = fLockerRanks[irow*nIV + icol];
	  fLockerRanks[irow*nIV + icol] = fLockerRanks[jrow*nIV+icol];
	  fLockerRanks[jrow*nIV+icol] = swap_row;
	}
      }
    }
  }
  cout << " fLockerRanks --" ;PrintVector(fLockerRanks);
  cout << " fAbsValue -- " ; PrintVector(fAbsValue);
  vector<Int_t> kUp(nIV);

  // Searching for the largest gap between fAbsValue
  for(int icol=0;icol<nIV;icol++){
    double max_gap=-1;
    Int_t kMax = 1;
    for(int irow=1;irow<nIV;irow++){
      if( fAbsValue[(irow-1)*nIV+icol] - fAbsValue[irow*nIV+icol] > max_gap){
	max_gap = fAbsValue[(irow-1)*nIV+icol] - fAbsValue[irow*nIV+icol];
	kMax = irow;
      }
    }
    kUp[icol] = kMax;
  }
  cout << " -- kUp " ;
  for(int icol=0;icol<nIV;icol++){
    cout << kUp[icol] << ",";
  }
  cout << endl;
  for(int irow=0;irow<nIV;irow++){
    for(int icol=0;icol<nIV;icol++){
      double row_length =0.0;
      double col_length =0.0;
      double sign_flipper=0.0;
      for(int iele=0;iele<kUp[icol];iele++){
	int rank_id = fLockerRanks[iele*nIV+icol] ;
	sign_flipper+=fRaw[rank_id*nIV+irow]*fLocker[rank_id*nIV+icol];
      }
      double fSign = 1.0;
      if(sign_flipper>0.0)
	fSign = 1.0;
      else
	fSign = -1.0;

      double fScaling =0.0;
      for(int iele=0;iele<kUp[icol];iele++){
	int rank_id = fLockerRanks[iele*nIV+icol] ;
	fDistance[irow*nIV+icol]+= pow(fSign*fRaw[rank_id*nIV+irow] - fLocker[rank_id*nIV+icol],2); //
	fCorrelation[irow*nIV+icol]+=fRaw[rank_id*nIV+irow]*fLocker[rank_id*nIV+icol]; //
	fScaling += pow( fLocker[rank_id*nIV+icol],2);
      }
      fCorrelation[irow*nIV+icol]/=fScaling;
    }
  }
  cout << "fRaw -- " ; PrintVector(fRaw);
  cout << "fLocker -- " ; PrintVector(fLocker);
  cout << "fCorrelation -- "; PrintVector(fCorrelation);
  cout << "fDistance -- "; PrintVector(fDistance);
  // = = === V3
  // vector<Int_t> fRow(nIV);
  // vector<Int_t> fCol(nIV);
  // for(int i=0; i<nIV;i++){
  //   fRow[i] = i; 
  //   fCol[i] = i; 
  // }
  // while(fRow.size()!=0){
  //   auto iter_row = fRow.begin();
  //   auto iter_col = fCol.begin();
  //   auto max_row = iter_row;
  //   auto max_col = iter_col;
  //   double max_correlation = 1e5;
  //   double my_correlation;
  //   while( iter_row!=fRow.end() ) {
  //     while( iter_col!=fCol.end()){
  // 	my_correlation = fDistance[ (*iter_row )*nIV +(*iter_col) ];
  // 	if(fabs(my_correlation)<fabs(max_correlation)){
  // 	  max_correlation = my_correlation;
  // 	  max_row = iter_row;
  // 	  max_col = iter_col;
  // 	}
  // 	iter_col++;
  //     }
  //     iter_row++;
  //   }
  //   fMap[ *max_row ] = *max_col;
  //   fRow.erase(max_row);
  //   fCol.erase(max_col);
  // }
  // = = === V1  
  vector<Int_t> fCandidate(nIV);
  for(int i=0; i<nIV;i++)
    fCandidate[i] = i; // row index

  for(int icol=0;icol<nIV;icol++){
    auto iter=fCandidate.begin();
    double max_correlation = 0.0;
    double min_distance = 1e5;
    auto max_iter = iter;
    while(iter!=fCandidate.end()){
      int my_idx = *iter;
      // if( fabs(fCorrelation[my_idx*nIV + icol] ) >fabs(max_correlation)){
      if( fDistance[my_idx*nIV + icol] < min_distance){
  	// max_correlation = fCorrelation[my_idx*nIV + icol];
  	min_distance = fDistance[my_idx*nIV+icol];
  	max_iter  = iter;
      }
      iter++;
    }
    fMap[ *max_iter ]  = icol;
    fCandidate.erase( max_iter);
  }
  
  // = = === v0  
  // vector<Int_t> fCandidate(nIV);
  // for(int i=0; i<nIV;i++)
  //   fCandidate[i] = i;
  // for(int irow=0;irow<nIV;irow++){
  //   double correlation_max = 0.0;
  //   int icol_candidate = -1;
  //   int iter_counter = 0;
  //   auto iter = fCandidate.begin();
  //   while(iter!=fCandidate.end()){
  //     int cur_icol = (*iter);
  //     if( fabs(fCorrelation[irow*nIV+cur_icol]) > fabs(correlation_max) ){
  // 	icol_candidate  = cur_icol;
  // 	iter_counter = iter - fCandidate.begin();
  // 	correlation_max = fCorrelation[irow*nIV+cur_icol];
  //     }
  //     iter++;
  //   }
  //   if(correlation_max==0.0){
  //     iter_counter = 0;
  //     icol_candidate = fCandidate[iter_counter];
  //   }
    
  //   fCandidate.erase( fCandidate.begin() + iter_counter );
  //   fMap[irow] = icol_candidate;
  //   fCorrelation_max[irow ] = correlation_max;
  // }
  cout << " -- Map:[" ;
  for(int irow=0;irow<nIV;irow++)
    cout << fMap[irow] << "," ;
  cout << "]" << endl;
  // cout << " -- Corr:[" ;
  // for(int irow=0;irow<nIV;irow++)
  //   printf("%.1f,", fCorrelation_max[irow] );
  // cout << "]" << endl;
  return fMap;
}


vector<Double_t> RemapVectors(vector<Double_t> fRaw, vector<Int_t> fmap){
  Int_t nLength = fRaw.size();
  vector< Double_t> fRet(nLength,0);
  Int_t nIV  = (Int_t)sqrt((Double_t)nLength);
  // remap by columns 
  for(int icol=0;icol<nIV;icol++)
    for(int irow=0;irow<nIV;irow++)
      fRet[ irow*nIV + fmap[icol] ] = fRaw[ irow*nIV+ icol ];

  return fRet;
}

vector<Double_t> FlipVectors(vector<Double_t> fRaw, vector<Int_t> fSign){
  Int_t nLength = fRaw.size();
  vector< Double_t> fRet(nLength,0);
  Int_t nIV  = (Int_t)sqrt((Double_t)nLength);
  
  for(int icol=0;icol<nIV;icol++)
    for(int irow=0;irow<nIV;irow++)
      fRet[ irow*nIV + icol ] = fSign[icol]* fRaw[ irow*nIV+ icol ];

  return fRet;
}

void PrintMatrixD(TMatrixD input){
  Int_t nRow = input.GetNrows();
  Int_t nCol = input.GetNcols();
  printf("------Printing MatrixD \n");
  for(int irow=0;irow<nRow;irow++){
    for(int icol=0;icol<nCol;icol++){
      if(fabs(input[irow][icol])<0.1)
	printf("    |");
      else{
	if(input[irow][icol]>=0)
	  printf(" ");
	printf( "%.1f|" , input[irow][icol] );
      }
    }
    printf("\n");
  }
  printf("------ \n");

}

void PrintVector(vector<Double_t> input){
  Int_t nLength = input.size();
  Int_t nRow = (Int_t) sqrt( nLength);
  Int_t nCol = nRow;
  printf("------Printing vector2d \n");
  for(int irow=0;irow<nRow;irow++){
    for(int icol=0;icol<nCol;icol++){
      if(input[irow*nCol+icol]>0)
	printf(" ");
      printf( "%.2f | " , input[irow*nCol + icol] );
    }
    printf("\n");
  }
  printf("------ \n");

}

void PrintVector(vector<Int_t> input){
  Int_t nLength = input.size();
  Int_t nRow = (Int_t) sqrt( nLength);
  Int_t nCol = nRow;
  printf("------Printing vector2d Integer \n");
  for(int irow=0;irow<nRow;irow++){
    for(int icol=0;icol<nCol;icol++){
      if(input[irow*nCol+icol]>0)
	printf(" ");
      printf( "%d | " , input[irow*nCol + icol] );
    }
    printf("\n");
  }
  printf("------ \n");

}


vector<Int_t> CheckSignWithRingAvg( vector<Double_t>  fRingAvg, 
				     vector<Double_t> fRaw ){
  double fThreshold =0.1;
  Int_t nLength = fRaw.size();
  Int_t nDim = (Int_t) sqrt( nLength);
  vector<Int_t> fSign(nDim,1);
  for(int icol=0;icol<nDim;icol++){
    double my_prod = 0;
    for(int irow=0;irow<nDim;irow++){
      my_prod  += fRaw[irow*nDim + icol] * fRingAvg[irow*nDim+icol];
    }
    if(fabs(my_prod)>fThreshold){
      if(my_prod>0)
	fSign[icol]=1;
      else
	fSign[icol]=-1;
    }else{
      fSign[icol]=1;
    }
  }
  return fSign;
}

map<Int_t, pair<Int_t,Int_t> > LoadRunInfo(){
  map<Int_t,pair<Int_t,Int_t> > fret;
  TString info_filename ="prex-runlist/all_production.list";
  ifstream slug_info;
  slug_info.open(info_filename.Data());
  TString sline;
  while(sline.ReadLine(slug_info)){
    TObjArray *token = sline.Tokenize(',');
    Int_t run_number = (((TObjString*)(token->At(0)))->GetString()).Atoi();
    Int_t slug_number = (((TObjString*)(token->At(1)))->GetString()).Atoi();
    Int_t arm_flag = (((TObjString*)(token->At(6)))->GetString()).Atoi();
    TString run_flag = ((TObjString*)(token->At(3)))->GetString();
    if(run_flag=="Good")
      fret[run_number]=make_pair(slug_number,arm_flag);
  }
  slug_info.close();
  return fret;

}

vector<Double_t> GetRingAverage( vector< vector<Double_t> > fRing) {

  Int_t nLength = fRing.size();
  Int_t nElement = fRing[0].size();
  vector<Double_t> fRingAvg(nElement,0.0);
  
  for(int ievt=0;ievt<nLength;ievt++){
    for(int iele=0;iele<nElement;iele++){
      fRingAvg[iele] += fRing[ievt][iele] ;
    }
  }
  
  for(int iele=0;iele<nElement;iele++){
    fRingAvg[iele] = fRingAvg[iele]/nLength ;
  }

  Int_t nDim = (int)sqrt(nElement);
  for(int icol=0;icol<nDim;icol++){
    double norm_factor=0.0;
    for(int irow=0;irow<nDim;irow++){
      norm_factor += pow(fRingAvg[irow*nDim +icol],2);
    }
    norm_factor = sqrt(norm_factor);
    for(int irow=0;irow<nDim;irow++){
      fRingAvg[irow*nDim +icol] = fRingAvg[irow*nDim+icol]/norm_factor;
    }
  }
  // PrintVector(fRingAvg);
  return fRingAvg;
}

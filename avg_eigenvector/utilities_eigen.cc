/*
  author: Tao Ye <tao.ye@stonybrook.edu>
  last update: July 2020

 */
#include "TaAccumulator.cc"

void PrintMatrixD(TMatrixD input);
void PrintVector(vector<Double_t> input);
void PrintVector(vector<Int_t> input);

vector<Int_t> CheckIdentity(vector<Double_t> fLocker,vector<Double_t> fRaw);
vector<Int_t> CheckIdentityWithSlug(vector<Double_t> fLocker,vector<Double_t> fRaw);
vector<Int_t> CheckSignWithRingAvg( vector<Double_t>  fRingAvg,
				     vector<Double_t> fRaw );
vector<Double_t> RemapVectors(vector<Double_t> fRaw, vector<Int_t> fmap);
vector<Double_t> FlipVectors(vector<Double_t> fRaw, vector<Int_t> fSign);
map<Int_t, pair<Int_t,Int_t> > LoadRunInfo();
vector<Double_t> GetRingAverage( vector< vector<Double_t> > fRing) ;
vector<vector<Int_t> > GeneratePermutations(Int_t nLength);

vector<Int_t> CheckIdentityWithSlug(vector<Double_t> fLocker,
				    vector<Double_t> fRaw){
  Int_t nLength = fRaw.size();
  Int_t nIV  = (Int_t)sqrt((Double_t)nLength);
  vector<Int_t> fMap(nIV);
  vector<Double_t> fDistance(nIV*nIV,0.0);
  vector<Int_t> fLockerMask(nIV*nIV,0);
  
  for(int icol=0;icol<nIV;icol++){
    vector<Int_t> fLockerRanks(nIV,0);
    vector<Double_t> fAbsValue(nIV);
    for(int irow=0;irow<nIV;irow++){
      fAbsValue[irow] = fabs(fLocker[irow*nIV+icol]);
      fLockerRanks[irow] = irow;
    }
    for(int irow=0;irow<nIV;irow++){
      double max_val = fAbsValue[irow];
      for(int jrow=irow+1;jrow<nIV;jrow++){
	if( fAbsValue[jrow] > max_val) {
	  max_val = fAbsValue[jrow];
	  fAbsValue[jrow] = fAbsValue[irow];
	  fAbsValue[irow] = max_val;
	  int swap_row = fLockerRanks[irow];
	  fLockerRanks[irow] = fLockerRanks[jrow];
	  fLockerRanks[jrow] = swap_row;
	}
      }
    }

    for(int irow=0;irow<nIV;irow++){
      fLockerMask[fLockerRanks[irow]*nIV + icol] =  1;
    }
  }
  for(int irow=0;irow<nIV;irow++){
    for(int icol=0;icol<nIV;icol++){
      double sign_flipper=0.0;
      for(int iele=0;iele<nIV;iele++){
      	sign_flipper+=fRaw[iele*nIV+irow]*fLocker[iele*nIV+icol]*fLockerMask[iele*nIV+icol];
      }
      double fSign = 1.0;
      if(sign_flipper>0.0)
      	fSign = 1.0;
      else
      	fSign = -1.0;
     
      for(int iele=0;iele<nIV;iele++){
	fDistance[irow*nIV+icol]+= pow(fSign*fRaw[iele*nIV+irow] - fLocker[iele*nIV+icol],2) *fLockerMask[iele*nIV+icol];
      }

    }
  }
  // cout << "fLocker -- "; PrintVector(fLocker);  
  // cout << "fRaw -- "; PrintVector(fRaw);  
  // cout << "fDistance -- "; PrintVector(fDistance);  
  vector<Int_t> fRow(nIV);
  vector<Int_t> fCol(nIV);
  for(int i=0; i<nIV;i++){
    fRow[i] = i;
    fCol[i] = i;
  }
  vector<Int_t> idx_min_by_col(nIV,-1);
  for(int icol=0;icol<nIV;icol++){
    double min_val=1e5;
    for(int irow=0;irow<nIV;irow++){
      if(fDistance[irow*nIV+icol] < min_val){
	min_val = fDistance[irow*nIV+icol];
	idx_min_by_col[icol] = irow;
      }
    }
  }

  for(int idx=0;idx<nIV;idx++){
    Bool_t kFound = kFALSE;
    for(int jdx=0;jdx<nIV;jdx++){
      if( jdx== idx )
	continue;
      if(idx_min_by_col[jdx] == idx_min_by_col[idx]){
	kFound = kTRUE;
	break;
      }
    }
    
    if(!kFound){
      fMap[ idx_min_by_col[idx] ] = idx;
      fRow.erase( find( fRow.begin(), fRow.end(), idx_min_by_col[idx] ) ) ;
      fCol.erase( find( fCol.begin(), fCol.end(), idx) ) ;
    }
  }
  
  cout << " -- fRow:[" ;
  auto iter_row = fRow.begin();
  while(iter_row!=fRow.end()){
    cout << *iter_row << ",";
    iter_row ++;
  }
  cout << "]" << endl;

  cout << " -- fCol:[" ;
  auto iter_col = fCol.begin();
  while(iter_col!=fCol.end()){
    cout << *iter_col << ",";
      iter_col ++;
  }
  cout << "]" << endl;

  double min_testing = 1e5;
  vector<Int_t> fSolution;
  Int_t nlength = fRow.size();
  // Loop over all possible combinations and find the optimized distance
  vector< vector<Int_t> > fCombinations = GeneratePermutations( fRow.size());
  auto iter_combo = fCombinations.begin();
  while(iter_combo!= fCombinations.end()){
    vector< Int_t> fArrayIdx = *iter_combo;
    double my_distance_sum  = 0.0;
    for(int iele=0; iele<nlength;iele++){
      int my_row = fRow[ iele ] ;
      int my_col = fCol[ fArrayIdx[iele] ];
      my_distance_sum += fDistance[ my_row*nIV + my_col ];
    }
    if( my_distance_sum <min_testing ) {
      min_testing = my_distance_sum;
      fSolution = fArrayIdx;
    }
    iter_combo++;
  }
  
  for(int iele=0;iele<nlength;iele++)
    fMap[ fRow[iele] ] = fCol[ fSolution[iele] ];
  cout << " -- Map:[" ;
  for(int irow=0;irow<nIV;irow++){
    cout << fMap[irow] << "," ;
  }
  cout << "]" << endl;
  return fMap;
}


vector<Int_t> CheckIdentity(vector<Double_t> fLocker,
			    vector<Double_t> fRaw){

  Int_t nLength = fRaw.size();
  Int_t nIV  = (Int_t)sqrt((Double_t)nLength);
  vector<Int_t> fMap(nIV);
  vector<Double_t> fDistance(nIV*nIV,0.0);

  for(int irow=0;irow<nIV;irow++){
    for(int icol=0;icol<nIV;icol++){
      double row_length =0.0;
      double col_length =0.0;
      double sign_flipper=0.0;
      for(int iele=0;iele<nIV;iele++){
	sign_flipper+=fRaw[iele*nIV+irow]*fLocker[iele*nIV+icol];
      }
      double fSign = 1.0;
      if(sign_flipper>0.0)
	fSign = 1.0;
      else
	fSign = -1.0;

      for(int iele=0;iele<nIV;iele++){
	fDistance[irow*nIV+icol]+= pow(fSign*fRaw[iele*nIV+irow] - fLocker[iele*nIV+icol],2);
      }

    }
  }
  // cout << "fLocker -- "; PrintVector(fLocker);  
  // cout << "fRaw -- "; PrintVector(fRaw);  
  // cout << "fDistance -- "; PrintVector(fDistance);  
  vector<Int_t> fRow(nIV);
  vector<Int_t> fCol(nIV);
  for(int i=0; i<nIV;i++){
    fRow[i] = i;
    fCol[i] = i;
  }
  vector<Int_t> idx_min_by_col(nIV,-1);
  for(int icol=0;icol<nIV;icol++){
    double min_val=1e5;
    for(int irow=0;irow<nIV;irow++){
      if(fDistance[irow*nIV+icol] < min_val){
	min_val = fDistance[irow*nIV+icol];
	idx_min_by_col[icol] = irow;
      }
    }
  }

  for(int idx=0;idx<nIV;idx++){
    Bool_t kFound = kFALSE;
    for(int jdx=0;jdx<nIV;jdx++){
      if( jdx== idx )
	continue;
      if(idx_min_by_col[jdx] == idx_min_by_col[idx]){
	kFound = kTRUE;
	break;
      }
    }
    
    if(!kFound){
      fMap[ idx_min_by_col[idx] ] = idx;
      fRow.erase( find( fRow.begin(), fRow.end(), idx_min_by_col[idx] ) ) ;
      fCol.erase( find( fCol.begin(), fCol.end(), idx) ) ;
    }
  }
  
  cout << " -- fRow:[" ;
  auto iter_row = fRow.begin();
  while(iter_row!=fRow.end()){
    cout << *iter_row << ",";
    iter_row ++;
  }
  cout << "]" << endl;

  cout << " -- fCol:[" ;
  auto iter_col = fCol.begin();
  while(iter_col!=fCol.end()){
    cout << *iter_col << ",";
      iter_col ++;
  }
  cout << "]" << endl;

  double min_testing = 1e5;
  vector<Int_t> fSolution;
  Int_t nlength = fRow.size();
  // Loop over all possible combinations and find the optimized distance
  vector< vector<Int_t> > fCombinations = GeneratePermutations( fRow.size());
  auto iter_combo = fCombinations.begin();
  while(iter_combo!= fCombinations.end()){
    vector< Int_t> fArrayIdx = *iter_combo;
    double my_distance_sum  = 0.0;
    for(int iele=0; iele<nlength;iele++){
      int my_row = fRow[ iele ] ;
      int my_col = fCol[ fArrayIdx[iele] ];
      my_distance_sum += fDistance[ my_row*nIV + my_col ];
    }
    if( my_distance_sum <min_testing ) {
      min_testing = my_distance_sum;
      fSolution = fArrayIdx;
    }
    iter_combo++;
  }
  
  for(int iele=0;iele<nlength;iele++)
    fMap[ fRow[iele] ] = fCol[ fSolution[iele] ];
  cout << " -- Map:[" ;
  for(int irow=0;irow<nIV;irow++){
    cout << fMap[irow] << "," ;
  }
  cout << "]" << endl;
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
map<Int_t, Int_t > LoadSpinInfo(){
  map<Int_t,Int_t > fret;
  TString info_filename ="prex-runlist/all_production.list";
  ifstream slug_info;
  slug_info.open(info_filename.Data());
  TString sline;
  while(sline.ReadLine(slug_info)){
    TObjArray *token = sline.Tokenize(',');
    Int_t run_number = (((TObjString*)(token->At(0)))->GetString()).Atoi();
    TString ihwp = ((TObjString*)(token->At(4)))->GetString();
    TString wien = ((TObjString*)(token->At(5)))->GetString();
    TString run_flag = ((TObjString*)(token->At(3)))->GetString();
    Int_t kHelicitySign;
    if( ihwp == "IN" )
      kHelicitySign = 1;
    else if (ihwp == "OUT" )
      kHelicitySign =-1;
    if( wien == "FLIP-RIGHT" )
      kHelicitySign *= 1;
    else if ( wien == "FLIP-LEFT" )
      kHelicitySign *=-1;
    
    if(run_flag=="Good")
      fret[run_number]= kHelicitySign;
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
  cout << "--fRingAvg";
  PrintVector(fRingAvg);
  return fRingAvg;
}


vector<vector<Int_t> > GeneratePermutations(Int_t nLength){
  
  vector< vector<Int_t> > fRet;
  if(nLength<=1)
    return fRet;
  
  Int_t index  =1 ;  
  vector<Int_t> fOne(1,0);
  vector< vector<Int_t> > fInput(1, fOne);
  while(index<nLength){
    vector< vector<Int_t> > fOutput;
    Int_t kRow = fInput.size();
    Int_t kLength = fInput[0].size();
    for(int irow=0;irow<kRow;irow++){
      for(int ipos=0;ipos<=kLength;ipos++){
	vector<Int_t> fCopy = fInput[irow];
	fCopy.insert( fCopy.begin()+ipos , index ) ;
	fOutput.push_back(fCopy);
      }
    }
    index++;
    fInput = fOutput;
  }

  fRet = fInput;
  return fRet;

}

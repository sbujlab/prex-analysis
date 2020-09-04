/*
  author: Tao Ye <tao.ye@stonybrook.edu>
  last update: July 2020

 */
#include "TaAccumulator.cc"

void PrintMatrixD(TMatrixD input);
void PrintVector(vector<Double_t> input);

vector<Int_t> CheckIdentity(vector<Double_t> fLocker,vector<Double_t> fRaw);
vector<Int_t> CrossCheckWithTemplate(vector< vector<Double_t> > fTemplate,
				     vector<Double_t> fRaw );
vector<Int_t> CrossCheckWithTemplate(vector< vector<TaAccumulator> > fTemplate,
				     vector<Double_t> fRaw );

vector<Double_t> RemapVectors(vector<Double_t> fRaw, vector<Int_t> fmap);
vector<Double_t> FlipVectors(vector<Double_t> fRaw, vector<Int_t> fSign);
map<Int_t, pair<Int_t,Int_t> > LoadRunInfo();
vector<Double_t> GetRingAverage( vector< vector<Double_t> > fRing) ;
void UpdateTemplate( vector<Double_t> &fTemplate, vector<Double_t> fVector);
void UpdateTemplate( vector< vector<TaAccumulator> > &fTemplate, vector<Double_t> fEigenVectorRaw);



vector<Int_t> CheckIdentity(vector<Double_t> fLocker,
			    vector<Double_t> fRaw){

  Int_t nLength = fRaw.size();
  Int_t nIV  = (Int_t)sqrt((Double_t)nLength);
  vector<Int_t> fMap(nIV);  
  TMatrixD mRaw(nIV,nIV);
  TMatrixD mLocker(nIV,nIV);
  
  for(int icol=0;icol<nIV;icol++){
    for(int irow=0;irow<nIV;irow++){
      mRaw[irow][icol] = fRaw[ irow*nIV+ icol ];
      mLocker[irow][icol] = fLocker[ irow*nIV+ icol ];
    }
  }
  
  TMatrixD mRaw_trans(mRaw);
  mRaw_trans.T();
  TMatrixD transform(nIV,nIV);
  transform = mRaw_trans*mLocker;
  // PrintMatrixD(transform);
  vector<Int_t> fCandidate(nIV);
  for(int i=0; i<nIV;i++)
    fCandidate[i] = i;

  for(int irow=0;irow<nIV;irow++){
    double correlation = 0.0;
    int icol_candidate = -1;
    int iter_counter = 0;
    auto iter = fCandidate.begin();
    while(iter!=fCandidate.end()){
      int cur_icol = (*iter);
      if( transform[irow][cur_icol] > correlation && transform[irow][cur_icol]>0.5  ){
	icol_candidate  = cur_icol;
	iter_counter = iter - fCandidate.begin();
	correlation = transform[irow][cur_icol];
      }
      iter++;
    }
    if(correlation==0){
      iter_counter = 0;
      icol_candidate = fCandidate[iter_counter];
    }
    
    fCandidate.erase( fCandidate.begin() + iter_counter );
    fMap[irow] = icol_candidate;
    if(irow==0){
      cout << " -- Largest correlation " << correlation << " at : " << icol_candidate  << endl;
    }
  }
  cout << " -- Re-Map:[" ;
  for(int irow=0;irow<nIV;irow++)
    cout << fMap[irow] << "," ;
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
      printf( "%.1f | " , input[irow*nCol + icol] );
    }
    printf("\n");
  }
  printf("------ \n");

}

vector<Int_t> CrossCheckWithTemplate(vector< vector<TaAccumulator> > fTemplate,
				     vector<Double_t> fRaw){

  double fThreshold =0.5;
  Int_t nLength = fRaw.size();
  Int_t nDim = (Int_t) sqrt( nLength);
  vector<Int_t> fSign(nDim,1);
  for(int icol=0;icol<nDim;icol++){
    Int_t nTemplate = fTemplate.size();
    Double_t prod_sum_max = 0;
    Int_t kTemplateIndex=-1;
    for(int itp=0;itp<nTemplate;itp++){
      double my_prod = 0;
      double my_scaling = 0.0;
      for(int iele=0;iele<nDim;iele++){
	my_prod  += fRaw[iele*nDim + icol] *  fTemplate[itp][iele].GetMean();
	my_scaling += pow( fTemplate[itp][iele].GetMean(), 2);
      }
      if(fabs(my_prod)/sqrt(my_scaling)> fabs(prod_sum_max)){
	prod_sum_max = my_prod/sqrt(my_scaling);
	kTemplateIndex = itp;
      }
    }
    if(icol==1){
      cout << " ^^ " << prod_sum_max << ":";
      for(int iele=0;iele<nDim;iele++)
	printf("%.1f,",fTemplate[kTemplateIndex][iele].GetMean());
      cout  << endl;
    }
    if(fabs(prod_sum_max)>fThreshold){
      if(prod_sum_max>0)
	fSign[icol]=1;
      else
	fSign[icol]=-1;
    }else{
      // even it is not strongly correlated with templates,
      // I let them align in direction.
      // if(fabs(prod_sum_max)<1e-10)
      // 	fSign[icol]=1;
      // else if(prod_sum_max>0) 
      // 	fSign[icol]=1;
      // else if(prod_sum_max<0)
      fSign[icol]=1;
      
      cout << " -- maximum_prod : " << prod_sum_max << endl;
      cout << " -- New Eigenvector: (";
      vector< Double_t > fBuff(nDim);
      for(int i=0;i<nDim;i++){
      	fBuff[i]  = fRaw[i*nDim+icol] * fSign[icol];
      	if(fBuff[i]>0)
      	  printf(" ");
      	printf("%.1f,",fBuff[i]);
      }
      cout << ")" << endl;
    }
  }
  return fSign;
}

vector<Int_t> CrossCheckWithTemplate(vector< vector<Double_t> > fTemplate, 
				     vector<Double_t> fRaw ){
  double fThreshold =0.8;
  Int_t nLength = fRaw.size();
  Int_t nDim = (Int_t) sqrt( nLength);
  vector<Int_t> fSign(nDim,1);
  for(int icol=0;icol<nDim;icol++){
    Int_t nTemplate = fTemplate.size();
    Double_t prod_sum_max = 0;
    Int_t kTemplateIndex=-1;
    for(int itp=0;itp<nTemplate;itp++){
      double my_prod = 0;
      for(int iele=0;iele<nDim;iele++)
	my_prod  += fRaw[iele*nDim + icol] * fTemplate[itp][iele];
      if(fabs(my_prod)> fabs(prod_sum_max)){
	prod_sum_max = my_prod;
	kTemplateIndex = itp;
      }
    }
    if(fabs(prod_sum_max)>fThreshold){
      if(prod_sum_max>0)
	fSign[icol]=1;
      else
	fSign[icol]=-1;
    }else{
      // even it is not strongly correlated with templates,
      // I let them align in direction.
      if(fabs(prod_sum_max)<1e-10)
	fSign[icol]=1;
      else if(prod_sum_max>0) 
	fSign[icol]=1;
      else if(prod_sum_max<0)
	fSign[icol]=-1;
      
      cout << " -- maximum_prod : " << prod_sum_max << endl;
      cout << " -- New Eigenvector: (";
      vector< Double_t > fBuff(nDim);
      for(int i=0;i<nDim;i++){
      	fBuff[i]  = fRaw[i*nDim+icol] * fSign[icol];
      	if(fBuff[i]>0)
      	  printf(" ");
      	printf("%.1f,",fBuff[i]);
      }
      cout << ")" << endl;
    }
  }
  return fSign;
}

void UpdateTemplate(vector<Double_t> &fTemplate, 
		    vector<Double_t> fVector ){
  // FIXME later: think about using running average update
  Int_t nDim  = fVector.size();

  for(int i=0;i<nDim;i++)
    fTemplate[i] = fVector[i];

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

  return fRingAvg;
}

vector< vector<Double_t> > GetTemplate( vector<Double_t> input){
  vector< vector<Double_t> > fRet;
  Int_t nLength = input.size();
  Int_t nRow = (Int_t) sqrt( nLength);
  Int_t nCol = nRow;

  for(int icol=0;icol<nCol;icol++){
    vector<Double_t> aVector;
    for(int irow=0;irow<nRow;irow++){
      aVector.push_back( input[irow*nCol+icol] );
    }
    fRet.push_back(aVector);
  }

  return fRet;
}


void UpdateTemplate( vector< vector<TaAccumulator> > &fTemplate,
		     vector<Double_t> fEigenVectorRaw){
  Int_t length =  fEigenVectorRaw.size();
  Int_t nDim = (int) sqrt(length);
  for(int icol=0;icol<nDim;icol++){
    Bool_t kFound = kFALSE;
    double kMaxProd = 0.0;
    auto itert = fTemplate.begin();
    while( itert != fTemplate.end()){
      double prod_sum =0.0;
      double scaling_factor= 0.0;
      for(int i=0;i<nDim;i++){
	prod_sum  += fEigenVectorRaw[i*nDim + icol] * ((*itert)[i].GetMean() );
	scaling_factor += pow((*itert)[i].GetMean(), 2);
      }
      if( fabs(prod_sum)/sqrt(scaling_factor) > kMaxProd )
	kMaxProd = prod_sum/sqrt(scaling_factor);
      
      if(prod_sum/sqrt(scaling_factor)>0.8){  // Assuming it is sign corrected
	kFound = kTRUE;
	for(int i=0;i<nDim;i++){
	  (*itert)[i].Update( fEigenVectorRaw[i*nDim + icol] );
	}
	break;
      }
    
      itert++;
    }
    if(!kFound){
      vector<TaAccumulator> fAccumulatorArray(length);
      for(int i=0;i<nDim;i++)
	fAccumulatorArray[i].Update( fEigenVectorRaw[i*nDim +icol] );
      fTemplate.push_back(fAccumulatorArray);
      cout << " -- Update : largest correlation "<< kMaxProd << endl;
      cout << " -- New Eigenvector: (";
      for(int i=0;i<nDim;i++){
      	if(fEigenVectorRaw[i*nDim + icol]>0)
      	  printf(" ");
      	printf("%.1f,",fEigenVectorRaw[i*nDim + icol]);
      }
      cout << ")" << endl;
      cout << " --  Template size:  " << fTemplate.size() << endl;
    }
  } // end of loop over eigenvectors
    
}
  

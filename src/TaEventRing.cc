#include "TaEventRing.hh"
TaEventRing::TaEventRing(){
  fRING_SIZE = 200;
  
  fHoldOff = 2000; // for SAM cooling down
  
  fThreshold = 0.05;
  fBeamOffLimit = 2.5; // default;
  
  fNextToBeFilled=0;
  fNextToRead=0;
  fNumberOfEvent=0;

  vector<Double_t> fDataElement; // [idet]
  for(int i =0; i<fRING_SIZE;i++){
    fDataArray.push_back(fDataElement);
    fBCMdata.push_back(0);
    fFlag.push_back(kTRUE);
  }
}

Bool_t TaEventRing::isReady(){
  if(fNumberOfEvent<fRING_SIZE)
    return kFALSE;
  return kTRUE;
}

void TaEventRing::PushBeamCurrent(Double_t input){
  fNumberOfEvent++;
  fBeamCurrent.Update(input);
  fBCMdata[fNextToBeFilled] = input;
  fFlag[fNextToBeFilled] = kTRUE;
  
  if(fCountDowns>0){
    fCountDowns--;
    for(int i=0;i<fRING_SIZE;i++)
      fFlag[i] = kFALSE;
  }

  if(input-fBeamCurrent.GetMean1()>fThreshold ||
     fBeamCurrent.GetMean1()>fBeamOffLimit ||
     fBeamCurrent.GetRMS()>fThreshold){
    
    for(int i=0;i<fRING_SIZE;i++)
      fFlag[i] = kFALSE;

    // for(int i=0;i<fPreCut;i++){
    //   Int_t iprev = fNextToBeFilled-i;
    //   iprev = iprev*fRING_SIZE;
    //   fFlag[iprev] = kFALSE;
    // }
    
    fCountDowns=fHoldOff;
  }
  
}

void TaEventRing::PushDetector(vector<Double_t> input){
  fDataArray[fNextToBeFilled] = input;
  // cout << "NTF: " << fNextToBeFilled ;  
  fNextToBeFilled++;
  fNextToBeFilled = fNextToBeFilled%fRING_SIZE;

}

Bool_t TaEventRing::Pop(vector<Double_t> &fOutput){
  // cout << " and NTR: " << fNextToRead << endl;  
  fOutput = fDataArray[fNextToRead];
  Bool_t fret = fFlag[fNextToRead];
  fBeamCurrent.DeAccumulate(fBCMdata[fNextToRead]);
  
  fNextToRead++;
  fNextToRead = (fNextToRead%fRING_SIZE);

  return fret;
}

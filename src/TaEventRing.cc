#include "TaEventRing.hh"
TaEventRing::TaEventRing(){
  fRING_SIZE = 50;
  fHoldOff = 1000;
  fThreshold = 0.01;
  fNextToBeFilled=0;
  fNextToRead=0;
  fNumberOfEvent=0;
  fBeamOffLimit = 2.5;

  vector<Double_t> fDataElement; // [idet]
  for(int i =0; i<fRING_SIZE;i++){
    fDataArray.push_back(fDataElement);
    fBCMdata.push_back(0);
  }
}

Bool_t TaEventRing::isReady(){
  if(fNumberOfEvent<fRING_SIZE)
    return kFALSE;
  
  if(fCountDowns>0){
    Pop();
    // cout << " count down:" << fCountDowns << endl;
    fCountDowns--;
    return kFALSE;
  }
  
  if(fBeamCurrent.GetMean1()>fBeamOffLimit ||
     fBCMdata[fNextToBeFilled-1]>fBeamOffLimit ||
     fabs(fBCMdata[fNextToBeFilled-1]-fBeamCurrent.GetMean1())>fThreshold){
    Pop();
    fCountDowns=fHoldOff;
    return kFALSE;
  }
  
  return kTRUE;
}

void TaEventRing::PushBeamCurrent(Double_t input){
  fNumberOfEvent++;
  fBeamCurrent.Update(input);
  fBCMdata[fNextToBeFilled] = input;
}

void TaEventRing::PushDetector(vector<Double_t> input){
  fDataArray[fNextToBeFilled] = input;
  // cout << "NTF: " << fNextToBeFilled ;  
  fNextToBeFilled++;
  fNextToBeFilled = fNextToBeFilled%fRING_SIZE;

}

vector<Double_t> TaEventRing::Pop(){
  // cout << " and NTR: " << fNextToRead << endl;  
  vector<Double_t> fRet = fDataArray[fNextToRead];
  fBeamCurrent.DeAccumulate(fBCMdata[fNextToRead]);
  fNextToRead++;
  fNextToRead = (fNextToRead%fRING_SIZE);

  return fRet;
}

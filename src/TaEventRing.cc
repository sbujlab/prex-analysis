#include "TaEventRing.hh"
TaEventRing::TaEventRing(){
  fRING_SIZE = 30;
  // fHoldOff = 1800; // for SAM cooling down
  fHoldOff = 300;
  
  fThreshold = 0.05;
  fBeamOffLimit = 2.5; // default;
  
  fNextToBeFilled=0;
  fNextToRead=0;
  fNumberOfEvent=0;

  vector<Double_t> fDataElement; // [idet]
  for(int i =0; i<fRING_SIZE;i++){
    fEventCounterArray.push_back(0);
    fBCMdata.push_back(0);
    fFlag.push_back(kGood);
  }
}

Bool_t TaEventRing::isReady(){
  if(fNumberOfEvent<fRING_SIZE)
    return kFALSE;
  return kTRUE;
}

void TaEventRing::PushBeamCurrent(Double_t input){
  fNumberOfEvent++;
  Int_t kErrorFlag =0;
  
  Int_t fRingErrorFlag=0;
  for(int i=0;i<fRING_SIZE;i++)
    fRingErrorFlag |=fFlag[i];

  if(fabs(input-fBeamCurrent.GetMean1())>fThreshold){
    kErrorFlag |= kBurpCut;
  }

  fBeamCurrent.Update(input);
  fBCMdata[fNextToBeFilled] = input;
  if(input>fBeamOffLimit || 
     fBeamCurrent.GetMean1()>fBeamOffLimit){
    fCountDowns=fHoldOff;
    kErrorFlag |=kAboveThreshold;
  }
  
  if(fBeamCurrent.GetRMS()>fThreshold){
    kErrorFlag |=kStabitlityError;
  }
  

  fFlag[fNextToBeFilled] = kErrorFlag;
  
  if(fRingErrorFlag==0 && (kErrorFlag&kAboveThreshold)==kAboveThreshold)
    fFlag[fNextToBeFilled] |= kBeamRecovery;

  if((fRingErrorFlag&kAboveThreshold)==kAboveThreshold &&
     (fRingErrorFlag&kBeamTripPoint)!=kBeamTripPoint &&
     kErrorFlag==0 )
    fFlag[fNextToBeFilled] |= kBeamTripPoint;

  
  for(int i=0;i<fRING_SIZE;i++)
    fFlag[i] |= kErrorFlag;

  if(fCountDowns>0){
    fCountDowns--;
    for(int i=0;i<fRING_SIZE;i++)
      fFlag[i] |= kHoldOff;
  }

  
}

// void TaEventRing::PushDetector(vector<Double_t> input){
//   fDataArray[fNextToBeFilled] = input;
//   // cout << "NTF: " << fNextToBeFilled ;  
//   fNextToBeFilled++;
//   fNextToBeFilled = fNextToBeFilled%fRING_SIZE;

// }

void TaEventRing::PushEventCounter(Double_t input){
  fEventCounterArray[fNextToBeFilled] = input;
  // cout << "NTF: " << fNextToBeFilled ;  
  fNextToBeFilled++;
  fNextToBeFilled = fNextToBeFilled%fRING_SIZE;

}

Int_t TaEventRing::PopEventCounter(Double_t &fOutput){
  // cout << " and NTR: " << fNextToRead << endl;  
  fOutput = fEventCounterArray[fNextToRead];
  Int_t fret = fFlag[fNextToRead];
  fBeamCurrent.DeAccumulate(fBCMdata[fNextToRead]);
  
  fNextToRead++;
  fNextToRead = (fNextToRead%fRING_SIZE);

  return fret;
}

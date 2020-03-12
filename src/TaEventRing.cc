#include "TaEventRing.hh"
TaEventRing::TaEventRing(){
  fRING_SIZE = 50;
  fHoldOff = fRING_SIZE;
}
TaEventRing::isReady(){
  if(fNumberOfEvent<fRING_SIZE)
    return kFALSE;

  if(fCountDowns>0){
    fCountDowns--;
    return kFALSE;
  }
    
  if(fBeamCurrent.GetRMS()>fThreshold){
    fCountDowns=fHoldOff;
    return kFALSE;
  }
    
  return kTRUE;
}


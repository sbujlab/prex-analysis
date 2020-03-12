#ifndef __TaEventRing_hh__
#define __TaEventRing_hh__

#include "TaAccumulator.hh"
#include "Rtypes.h"
class TaAccumulator;
class TaEventRing{
public:
  TaEventRing(){};
  virtual ~TaEventRing(){};
  oid PushBeamCurrent(Double_t input);
  Double_t Pop(Int_t i);
  Double_t Pop();
  Double_t Push(vector< Double_t> >);
  Bool_t isReady();

private:
  TaAccumulator fBeamCurrent;
  vector< vector<Double_t>  > 
  fRING_SIZE;
  fThreshold;
  fCountDowns;
  fHoldOff;
  fNumberOfEvent;
  
  ClassDef(TaEventRing,0);
};

#endif

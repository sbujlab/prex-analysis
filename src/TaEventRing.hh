#ifndef __TaEventRing_hh__
#define __TaEventRing_hh__

#include "TaAccumulator.hh"
#include "Rtypes.h"
class TaAccumulator;
class TaEventRing{
public:
  TaEventRing();
  virtual ~TaEventRing(){};
  void PushBeamCurrent(Double_t input);
  vector<Double_t> Pop();
  void PushDetector(vector< Double_t>);
  Bool_t isReady();

private:
  TaAccumulator fBeamCurrent;
  vector< vector<Double_t>  > fDataArray; // [ievt][idet]
  vector<Double_t> fBCMdata;
  Double_t fThreshold;
  
  Int_t fRING_SIZE;
  Int_t fCountDowns;
  Int_t fHoldOff;
  Int_t fNextToBeFilled;
  Int_t fNextToRead;
  Int_t fNumberOfEvent;
  ClassDef(TaEventRing,0);
};

#endif

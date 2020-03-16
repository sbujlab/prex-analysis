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
  Bool_t Pop(vector<Double_t> &);
  void PushDetector(vector< Double_t>);
  Bool_t isReady();
  inline void SetBeamOffLimit(double limit){ fBeamOffLimit = limit;};
private:
  TaAccumulator fBeamCurrent;
  vector< vector<Double_t>  > fDataArray; // [ievt][idet]
  vector<Double_t> fBCMdata;
  vector<Bool_t> fFlag; // [ievt]
  
  Double_t fThreshold;
  Double_t fBeamOffLimit;

  Int_t fRING_SIZE;
  Int_t fCountDowns;
  Int_t fHoldOff;
  Int_t fNextToBeFilled;
  Int_t fNextToRead;
  Int_t fNumberOfEvent;
  ClassDef(TaEventRing,0);
};

#endif

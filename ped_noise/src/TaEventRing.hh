#ifndef __TaEventRing_hh__
#define __TaEventRing_hh__

#include "TaAccumulator.hh"
#include "Rtypes.h"
#include "TaEventTypes.hh"
class TaAccumulator;
class TaEventRing{
public:
  TaEventRing();
  virtual ~TaEventRing(){};
  void PushBeamCurrent(Double_t input);
  // Bool_t Pop(vector<Double_t> &);
  // void PushDetector(vector< Double_t>);
  void PushEventCounter(Double_t);
  Int_t PopEventCounter(Double_t &);
  Bool_t isReady();
  
  inline void SetBeamOffLimit(double limit){ fBeamOffLimit = limit;};
  inline void SetRingSize(Int_t ring_size){ fRING_SIZE = ring_size;};
  inline void SetHoldOff(Int_t hold_off){ fHoldOff = hold_off;};
  
private:
  TaAccumulator fBeamCurrent;
  // vector< vector<Double_t>  > fDataArray; // [ievt][idet]
  vector<Double_t> fEventCounterArray;
  vector<Double_t> fBCMdata;
  vector<Int_t> fFlag; // [ievt]
  
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

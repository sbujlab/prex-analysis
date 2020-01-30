#ifndef __TaRunInfo_v2_hh_
#define __TaRunInfo_v2_hh_
class TaRunInfo: public TObject{
public:
  TaRunInfo();
  virtual ~TaRunInfo(){};
  Bool_t ParseLine(TString);
  Int_t ParseTargetName(TString);
  inline Int_t GetTargetID() const {return target_id;};
  inline TString GetRunFlag() const {return run_flag;};
  inline Int_t GetArmFlag() const {return arm_flag;};
  inline Int_t GetRunNumber() const {return run_number;};
  inline Int_t GetSlugNumber() const {return slug_number;};
  inline TString GetIHWPStatus() const {return ihwp;};
  inline TString GetWienMode() const {return wien;};
  inline Int_t GetSign() const {return sign;};
private:
  Int_t run_number;
  Int_t slug_number;
  Int_t target_id;
  Int_t arm_flag;
  TString run_flag;
  TString ihwp;
  TString wien;
  Int_t sign;
  ClassDef(TaRunInfo,0);
};

#endif

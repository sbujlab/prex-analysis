#ifndef __TaSumStat_hh__
#define __TaSumStat_hh__
typedef struct {Double_t hw_sum, hw_sum_m2,hw_sum_err,num_samples;} JAPAN_STAT;
typedef struct {Double_t mean, err,rms, num_samples;} POSTPAN_STAT;
typedef struct {Double_t mean,err,rms,num_samples;} SUM_STAT;
typedef struct {Double_t ppm,ppb,um,mm,nm;} UNIT;

using namespace std;

class TaSumStat{
public:
  TaSumStat();
  ~TaSumStat(){};
  JAPAN_STAT init_japan_stat();
  POSTPAN_STAT init_postpan_stat();
  JAPAN_STAT invalid_japan_stat();
  void merge_japan_stat(JAPAN_STAT&, JAPAN_STAT);
  void write_sum_stat(SUM_STAT&,JAPAN_STAT);
  void write_sum_stat(SUM_STAT&,POSTPAN_STAT);
  void write_sum_stat_by_name(TString);
  void write_sum_postpan_stat_by_name(TString);
  void construct_branches(TFile*);
  void collect_branchlist_from_input(TTree*);
  void collect_branchlist_from_postpan(TTree*);
  inline void set_run_number(Int_t input){ run_number = input;};
  inline void set_burst_counter(Int_t input){ burst_counter = input;};
  inline void set_minirun_counter(Int_t input){ minirun_counter = input;};
  void write_trees_to_output(TFile*);

  void load_null_stat_by_name(TString);
  void load_null_postpan_stat_by_name(TString);
  void load_invalid_stat_by_name(TString);
  void load_japan_stat_ptr(TTree* );
  void load_postpan_stat_ptr(TTree* );
  void fill_tree_by_name(TString);
  void cache_japan_stat(TString);
  void merge_japan_stat(TString);
private:
  map<TString,vector<TString> > fBranchNameListMap;
  vector<TString> fTreeNameList;
  map<TString, TTree*> fTreeMap;
  vector<TTree*> fTreeArray;
  map<TString, vector<SUM_STAT> > fSumStatMap;
  map<TString, vector<JAPAN_STAT> > fJStatMap;
  map<TString, vector<POSTPAN_STAT> > fPStatMap;
  vector<JAPAN_STAT> fJStatBuffer;
  UNIT aUnit;
  Int_t run_number;
  Int_t burst_counter;
  Int_t minirun_counter; // they are the same actually;
};

#endif

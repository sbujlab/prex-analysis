typedef struct {Double_t hw_sum, hw_sum_m2,hw_sum_err,num_samples;} JAPAN_STAT;
typedef struct {Double_t mean,err,rms,num_samples;} SUM_STAT;
typedef struct {Double_t ppm,ppb,um,mm,nm;} UNIT;

void init_japan_stat(JAPAN_STAT &this_stat){
  this_stat.hw_sum=0.0;
  this_stat.hw_sum_m2=0.0;
  this_stat.hw_sum_err=0.0;
  this_stat.num_samples=0.0;
}
void update_japan_stat(JAPAN_STAT &dest_stat,JAPAN_STAT in_stat){
  double mean_1 = dest_stat.hw_sum;
  double m2_1  = dest_stat.hw_sum_m2;
  double nsamp_1 = dest_stat.num_samples;

  double mean_2 = in_stat.hw_sum;
  double m2_2  = in_stat.hw_sum_m2;
  double nsamp_2 = in_stat.num_samples;
  double delta_mean = mean_2 - mean_1;

  m2_1 += m2_2;
  m2_1 += nsamp_1*nsamp_2*pow(delta_mean,2)/(nsamp_1+nsamp_2);

  dest_stat.hw_sum += nsamp_2*delta_mean/(nsamp_1+nsamp_2);
  dest_stat.num_samples+=nsamp_2;
  dest_stat.hw_sum_m2 = m2_1;
  dest_stat.hw_sum_err = TMath::Sqrt(dest_stat.hw_sum_m2)/dest_stat.num_samples;
}

void write_sum_stat(SUM_STAT &dest_stat,JAPAN_STAT in_stat){
  dest_stat.mean = in_stat.hw_sum;
  dest_stat.err = in_stat.hw_sum_err;
  dest_stat.rms = TMath::Sqrt(in_stat.hw_sum_m2/in_stat.num_samples);
  dest_stat.num_samples  = in_stat.num_samples;
}

void construct_branches(TTree* atree, vector<TString> fBranchNameList){
  TString leaflist = "mean/D:err:rms:num_samples";
  UNIT aUnit;
  aUnit.ppm = 1e-6;
  aUnit.ppb = 1e-9;
  aUnit.mm = 1.0;
  aUnit.um = 1e-3;
  aUnit.nm = 1e-6;
  atree->Branch("unit",&aUnit,"ppm/D:ppb:um:mm:nm");

    for(int i=0;i<ndet;i++)
      muls_tree->Branch(det_list[i],&sum_val[i],leaflist);

}


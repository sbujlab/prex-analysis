typedef struct {Double_t hw_sum, hw_sum_m2,hw_sum_err,num_samples;} SUM_STAT;
typedef struct {Double_t mean,err,rms;} MINI_STAT;

void null_sum_stat(SUM_STAT &this_stat){
  this_stat.hw_sum=0.0;
  this_stat.hw_sum_m2=0.0;
  this_stat.hw_sum_err=0.0;
  this_stat.num_samples=0.0;
}
void update_sum_stat(SUM_STAT &dest_stat,SUM_STAT in_stat){
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
void update_mini_stat(MINI_STAT &dest_stat,SUM_STAT in_stat){
  dest_stat.mean = in_stat.hw_sum;
  dest_stat.err = in_stat.hw_sum_err;
  dest_stat.rms = TMath::Sqrt(in_stat.hw_sum_m2/in_stat.num_samples);
}

void MergeMulTree(Int_t slug);
void MergeMulTree(TString label);
void MergeMulTree();

void MergeMulTree(){
  for(int i=0;i<=94;i++){
    MergeMulTree(i);
  }
}
void MergeMulTree(Int_t slug){
  TString label = Form("slug%d",slug);
  MergeMulTree(label);
}

void MergeMulTree(TString label){
  TString filename  = label+".list";
  TString path = "./prex-runlist/simple_list/";
  FILE *runlist = fopen((path+filename).Data(),"r");
  
  if(runlist==NULL){
    cerr << " -- Error: runlist is not found! " << endl;
    return;
  }

  TFile *output = TFile::Open("./rootfiles/MulMerged_"+label+".root","RECREATE");
  TTree *muls_tree = new TTree("muls","Mult summary");
  vector<TString> det_list={"asym_bcm_an_us","asym_bcm_dg_us",
			    "asym_bcm_an_ds","asym_bcm_dg_ds","asym_bcm_an_ds3"};
  const Int_t ndet = det_list.size();
  MINI_STAT mini_val[ndet];
  SUM_STAT sum_val[ndet];
  for(int i=0;i<ndet;i++)
    muls_tree->Branch(det_list[i],&mini_val[i],"mean/D:err:rms");
  
  typedef struct {Double_t ppm,ppb,um,mm,nm;} UNIT;
  UNIT aUnit;
  aUnit.ppm = 1e-6;
  aUnit.ppb = 1e-9;
  aUnit.mm = 1.0;
  aUnit.um = 1e-3;
  aUnit.nm = 1e-6;
  muls_tree->Branch("unit",&aUnit,"ppm/D:ppb:um:mm:nm");
  Int_t run_id;
  muls_tree->Branch("run",&run_id);
  while(!feof(runlist)){
    Int_t run_number=0;
    fscanf(runlist,"%d/n",&run_number);
    if(run_number==0)
      continue;

    TFile *this_file = TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass1_%d.000.root",run_number));
    if(this_file==NULL)
      continue;
    cout << this_file->GetName() << endl;

    TTree* burst = (TTree*)this_file->Get("burst");
    if(burst==NULL)
      continue;

    SUM_STAT temp_val[ndet];
    for(int i=0;i<ndet;i++){
      null_sum_stat(sum_val[i]);
      TBranch* branch_ptr = burst->GetBranch(det_list[i]);
      branch_ptr->GetLeaf("hw_sum")->SetAddress(&temp_val[i].hw_sum);
      branch_ptr->GetLeaf("hw_sum_m2")->SetAddress(&temp_val[i].hw_sum_m2);
      branch_ptr->GetLeaf("hw_sum_err")->SetAddress(&temp_val[i].hw_sum_err);
      branch_ptr->GetLeaf("num_samples")->SetAddress(&temp_val[i].num_samples);
    }
    Int_t nBurst = burst->GetEntries();    
    for(int ievt=0;ievt<nBurst;ievt++){
      burst->GetEntry(ievt);
      for(int i=0;i<ndet;i++)
	update_sum_stat(sum_val[i],temp_val[i]);
    }
    for(int i=0;i<ndet;i++)
      update_mini_stat(mini_val[i],sum_val[i]);
    run_id = run_number;
    muls_tree->Fill();

    this_file->Close();
  }
  output->cd();
  muls_tree->Write();
  cout << " -- Closing " << output->GetName() << endl;
  output->Close();
}



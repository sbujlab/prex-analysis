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

void MergeDitTree(Int_t slug);
void MergeDitTree(TString label);
void MergeDitTree();

void MergeDitTree(){
  for(int i=0;i<=94;i++){
    MergeDitTree(i);
  }
}
void MergeDitTree(Int_t slug){
  TString label = Form("slug%d",slug);
  MergeDitTree(label);
}

void MergeDitTree(TString label){
  TString filename  = label+".list";
  TString path = "./prex-runlist/simple_list/";
  FILE *runlist = fopen((path+filename).Data(),"r");
  
  if(runlist==NULL){
    cerr << " -- Error: runlist is not found! " << endl;
    return;
  }

  TFile *output = TFile::Open("./rootfiles/DitMerged_"+label+".root","RECREATE");
  TTree *dit_sum_tree = new TTree("dits","Dithering summary");
  vector<TString> det_list={"dit_asym_usl","dit_asym_usr",
			    "dit_asym_us_avg","dit_asym_us_dd"};
  vector<TString> at_list ={"dit_asym_atr1","dit_asym_atl1","dit_asym_atr2","dit_asym_atl2",
			    "dit_asym_atl_avg","dit_asym_atl_dd","dit_asym_atr_avg","dit_asym_atr_dd",
			    "dit_asym_at1_avg","dit_asym_at1_dd","dit_asym_at2_avg","dit_asym_at2_dd",
			    "dit_asym_atl1r2_avg","dit_asym_atl1r2_dd","dit_asym_atr1l2_avg","dit_asym_atr1l2_dd",
			    "dit_asym_atl_dd_atr_dd_avg","dit_asym_atl_dd_atr_dd_dd"};
  det_list.insert(det_list.end(),at_list.begin(),at_list.end());
  const Int_t ndet = det_list.size();
  MINI_STAT mini_val[ndet];
  SUM_STAT sum_val[ndet];
  for(int i=0;i<ndet;i++)
    dit_sum_tree->Branch(det_list[i],&mini_val[i],"mean/D:err:rms");
  
  typedef struct {Double_t ppm,ppb,um,mm,nm;} UNIT;
  UNIT aUnit;
  aUnit.ppm = 1e-6;
  aUnit.ppb = 1e-9;
  aUnit.mm = 1.0;
  aUnit.um = 1e-3;
  aUnit.nm = 1e-6;
  dit_sum_tree->Branch("unit",&aUnit,"ppm/D:ppb:um:mm:nm");
  Int_t run_id;
  dit_sum_tree->Branch("run",&run_id);
  while(!feof(runlist)){
    Int_t run_number=0;
    fscanf(runlist,"%d/n",&run_number);
    if(run_number==0)
      continue;

    TFile *this_file = TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass1_%d.000.root",run_number));
    if(this_file==NULL)
      continue;
    cout << this_file->GetName() << endl;

    TTree* burst_mulc_dit = (TTree*)this_file->Get("burst_mulc_dit");
    TTree* burst_mulc_dit_combo = (TTree*)this_file->Get("burst_mulc_dit_combo");
    if(burst_mulc_dit_combo==NULL || burst_mulc_dit==NULL)
      continue;

    burst_mulc_dit->AddFriend(burst_mulc_dit_combo);

    SUM_STAT temp_val[ndet];
    for(int i=0;i<ndet;i++){
      null_sum_stat(sum_val[i]);
      TBranch* branch_ptr = burst_mulc_dit->GetBranch(det_list[i]);
      if(branch_ptr==NULL)
	continue;
      branch_ptr->GetLeaf("hw_sum")->SetAddress(&temp_val[i].hw_sum);
      branch_ptr->GetLeaf("hw_sum_m2")->SetAddress(&temp_val[i].hw_sum_m2);
      branch_ptr->GetLeaf("hw_sum_err")->SetAddress(&temp_val[i].hw_sum_err);
      branch_ptr->GetLeaf("num_samples")->SetAddress(&temp_val[i].num_samples);
    }
    Int_t nBurst = burst_mulc_dit->GetEntries();    
    for(int ievt=0;ievt<nBurst;ievt++){
      burst_mulc_dit->GetEntry(ievt);
      for(int i=0;i<ndet;i++)
	update_sum_stat(sum_val[i],temp_val[i]);
    }
    for(int i=0;i<ndet;i++)
      update_mini_stat(mini_val[i],sum_val[i]);
    run_id = run_number;
    dit_sum_tree->Fill();

    this_file->Close();
  }
  output->cd();
  dit_sum_tree->Write();
  cout << " -- Closing " << output->GetName() << endl;
  output->Close();
}



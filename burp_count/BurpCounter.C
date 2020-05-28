void BurpCounter(Int_t slug);
void BurpCounter(){
  for(int i=1;i<=94;i++)
    BurpCounter(i);
}
void BurpCounter(Int_t slug){
  TFile *output = TFile::Open(Form("slug%d.root",slug),"RECREATE");
  TTree *cnt_tree = new TTree("T","event counter tree");
  Int_t nGood,nGoodIgnoringBPM,nGoodButFailedByBPM;
  Int_t run_number;
  cnt_tree->Branch("run",&run_number);
  cnt_tree->Branch("ngood",&nGood);
  cnt_tree->Branch("ngood_wo_bpm",&nGoodIgnoringBPM);
  cnt_tree->Branch("nBPM_burp",&nGoodButFailedByBPM);

  TString runlist_filename = Form("./prex-runlist/simple_list/slug%d.list",slug);
  FILE *runlist = fopen(runlist_filename.Data(),"r");
  
  while(!feof(runlist)){
    fscanf(runlist,"%d\n",&run_number);
  
    TFile* rootfile = TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass1_%d.000.root",run_number));
    TTree *evt_tree = (TTree*)rootfile->Get("evt");

    nGood = evt_tree->GetEntries("ErrorFlag==0");
    nGoodIgnoringBPM = evt_tree->GetEntries("(ErrorFlag&0xdbfefbff)==0") ;
    nGoodButFailedByBPM = evt_tree->GetEntries("(ErrorFlag&0xdbfefbff)==0 && (ErrorFlag&0x24010400)==0x24010400");
    
    cout << " Total Good Events : "<<  nGood << endl;
    cout << " Good Events ignoring BPM: "<< nGoodIgnoringBPM << endl;
    cout << " Good Events Failed By BPM: " << nGoodButFailedByBPM << endl;

    cnt_tree->Fill();
    rootfile->Close();
  }
  output->cd();
  cnt_tree->Write();
  output->Close();
}

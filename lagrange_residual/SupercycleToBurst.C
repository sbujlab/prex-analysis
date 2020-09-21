/*
   author: Tao Ye
   Last update: July 2019
 */

void SupercycleToBurst(Int_t slug){

  TFile *output = TFile::Open(Form("rootfiles/slug%d_cycle_burst.root",slug),"RECREATE");
  TTree *cyc2b_tree = new TTree("cyc2b","dithering supercycle to burst counter");
  Int_t kRun, kCycle, kBmwObj;
  Int_t kBurstCounter;
  
  cyc2b_tree->Branch("run",&kRun);
  cyc2b_tree->Branch("bmwcycnum",&kCycle);
  cyc2b_tree->Branch("bmwobj",&kBmwObj);
  cyc2b_tree->Branch("BurstCounter",&kBurstCounter);
  TString list_name = Form("prex-runlist/good_list/slug%d.list",slug);
  FILE *runlist = fopen(list_name.Data(),"r");
  while(!feof(runlist)){
    kRun = 0;
    fscanf(runlist,"%d\n",&kRun);
    if(kRun==0)
      continue;
    TFile *japan_file= TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass2_%d.000.root",kRun));
    if(japan_file==NULL){
      cerr<< " ** Error: " << japan_file->GetName() 
	  << " is not found. " << endl;
      continue;
    }else{
      cout<< " -- Found: " << japan_file->GetName()  << endl;
    }
    TTree *mul_tree = (TTree*)japan_file->Get("mul");
    Double_t fCycle, fBmwObj;
    Short_t fBurstCounter;
    Double_t fCodaEventNumber;
    mul_tree->SetBranchAddress("yield_bmwcycnum",&fCycle);
    mul_tree->SetBranchAddress("yield_bmwobj",&fBmwObj);
    mul_tree->SetBranchAddress("BurstCounter",&fBurstCounter);
    TEventList *evlist = new TEventList("evlist");
    Int_t nEntries = mul_tree->Draw(">>evlist","(ErrorFlag&0x19000)==0x19000");
    cout << " -- Found : " << nEntries << " entries." << endl;
    Int_t last_cycle = -1 ;
    Int_t last_bmwobj = -1;
    Int_t last_burst = -1;
    for(int ievt=0;ievt<nEntries;ievt++){
      Int_t kEntry = evlist->GetEntry(ievt);
      mul_tree->GetEntry(kEntry);
      if(fCycle>last_cycle && fCycle>0){
	last_cycle = fCycle;
	kCycle = fCycle;
	cout << endl;
	cout << " -- Found Supercycle " << last_cycle << endl;
	last_bmwobj = -1; 
	cout << " -- Map Bmod Obj: " ;
      }
      if(fBmwObj>last_bmwobj && fBmwObj>0){
	cout << fBmwObj << "(" << fBurstCounter << ")  " ;
	last_bmwobj = fBmwObj;
	kBmwObj = fBmwObj;
	cyc2b_tree->Fill();
      }
    }
    japan_file->Close();
    cout << endl;
  }
  output->cd();
  cyc2b_tree->Write();
  output->Close();

}

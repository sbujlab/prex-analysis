void CheckCoilAmplitude(Int_t run_number);
void CheckCoilAmplitude(){
  FILE *list = fopen("all.list","r");
  while(!feof(list)){
    Int_t run_number=0;
    fscanf(list,"%d\n",&run_number);
    if(run_number!=0)
      CheckCoilAmplitude(run_number);
    cout << run_number << endl;
  }
  fclose(list);
}
void CheckCoilAmplitude(Int_t run_number){
  TFile *input = TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass1_%d.000.root",run_number));
  TTree *evt_tree = (TTree*)input->Get("evt");
  if(evt_tree==NULL)
    return;
  
  Int_t nCoils=7;
  Double_t bmwcycnum;
  evt_tree->SetBranchAddress("bmwcycnum",&bmwcycnum);
  TEventList *elist = new TEventList("elist");
  Int_t npt = evt_tree->Draw(">>+elist","bmwcycnum>0 && bmod_ramp>0 && bmwobj>0");
  if(npt==0)
    return;
  evt_tree->SetEventList(elist);
  TH1D *hcyc = new TH1D("hcyc","",5000,-0.5,4999.5);
  evt_tree->Draw("bmwcycnum>>hcyc","","goff");
  Int_t bin_first = hcyc->FindFirstBinAbove(0);
  Int_t bin_last = hcyc->FindLastBinAbove(0);
  Int_t cycID_first = hcyc->GetBinCenter(bin_first);
  Int_t cycID_last = hcyc->GetBinCenter(bin_last);
  Int_t nCycle = cycID_last-cycID_first+1;
  

  TFile *output = TFile::Open(Form("rootfiles/bmod_amplitude_%d.root",run_number),"RECREATE");
  TTree *bmod = new TTree("bmod","amplitude tree");
  Int_t cycID ;
  Int_t run = run_number ;
  bmod->Branch("cycID",&cycID);
  bmod->Branch("run",&run_number);
  vector<Double_t> amplitude(nCoils,0);
  for(int icoil=0;icoil<nCoils;icoil++)
    bmod->Branch(Form("coil%d",icoil+1),&amplitude[icoil]);

  for(int icyc=0;icyc<nCycle;icyc++){
    cycID = cycID_first+icyc; 
    for(int icoil=0;icoil<nCoils;icoil++){
      evt_tree->Draw(Form("bmod_trim%d",icoil+1),
		     Form("bmwobj==%d && bmwcycnum==%d",icoil+1,cycID_first+icyc),
		     "goff");
      TH1D* htemp=(TH1D*)gDirectory->FindObject("htemp");
      double rms=htemp->GetRMS();
      amplitude[icoil] = rms;
    }
    bmod->Fill();
  }
  bmod->Write();
  output->Close();
  input->Close();
}

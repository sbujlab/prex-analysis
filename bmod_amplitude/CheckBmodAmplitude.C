void CheckBmodAmplitude(Int_t run_number);
void CheckBmodAmplitude(){
  // FILE *list = fopen("all.list","r");
  // while(!feof(list)){
  //   Int_t run_number=0;
  //   fscanf(list,"%d\n",&run_number);
  //   if(run_number!=0)
  //     CheckBmodAmplitude(run_number);
  //   cout << run_number << endl;
  // }
  // fclose(list);
}
void CheckBmodAmplitude(Int_t run_number){
  TFile *input = TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass1_%d.000.root",run_number));
  TTree *evt_tree = (TTree*)input->Get("evt");
  if(evt_tree==NULL)
    return;
  
  Int_t nCoils=7;
  Double_t bmwcycnum;
  evt_tree->SetBranchAddress("bmwcycnum",&bmwcycnum);
  TEventList *elist = new TEventList("elist");
  Int_t npt = evt_tree->Draw(">>+elist","bmwcycnum>0 && bmod_ramp>0 && bmwobj>0 && (ErrorFlag&0xda7e6bff)==0");
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
  Int_t coilID ;
  Int_t run = run_number ;
  Double_t deg_to_rad = TMath::DegToRad();
  Double_t unit_um = 1e-3;
  bmod->Branch("cycID",&cycID);
  bmod->Branch("coil",&coilID);
  bmod->Branch("run",&run_number);
  bmod->Branch("deg",&deg_to_rad);
  bmod->Branch("um",&unit_um);
  vector<TString> fChannelArray={"usl","usr","dsl","dsr",
				 "atl1","atl2","atr1","atr2",
				 "bpm4aX","bpm4eX","bpm4aY","bpm4eY",
				 "bpm11X","bpm12X"};
  Int_t nch = fChannelArray.size();
  typedef struct {Double_t amplitude, phase, offset, frequency, chisq,ndf,ndata;} BMODSTAT;
  TString leaflist = "amplitdue/D:phase:offset:frequency:chisq:ndf:ndata";
  vector<BMODSTAT> fBmodStatArray(nch);
  for(int ich = 0; ich<nch;ich++)
    bmod->Branch(fChannelArray[ich],&(fBmodStatArray[ich]),leaflist);

  TF1 *fsin = new TF1("fsin","[0]*TMath::Sin([1]*x +[2])+[3]",-1e5,1e5);
  Double_t par_init[] =  {0,1.0/1000,0,0};
  fsin->SetParameters(par_init);
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  for(int icyc=0;icyc<nCycle;icyc++){
    cycID = cycID_first+icyc; 
    for(int icoil=0;icoil<nCoils;icoil++){
      int npt = evt_tree->Draw("bmod_ramp",
			       Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0",
				    icoil+1,cycID_first+icyc),"goff");
      
      if(npt==0)
	continue;

      for(int ich=0;ich<nch;ich++){
	fsin->SetParameters(par_init);

	evt_tree->Draw(Form("%s",fChannelArray[ich].Data()),
		       Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 ",
			    icoil+1,cycID_first+icyc),"goff");

	TH1D *h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	fsin->SetParameter(3,h1dptr->GetMean());
	fsin->SetParameter(1,1.414*h1dptr->GetRMS());
	evt_tree->Draw(Form("%s:bmod_ramp",fChannelArray[ich].Data()),
		       Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 ",
			    icoil+1,cycID_first+icyc),"*");
	TGraph *gbuff = (TGraph*)gPad->FindObject("Graph");
	gbuff->Fit("fsin","Q");
	fBmodStatArray[ich].amplitude = fsin->GetParameter(0);
	fBmodStatArray[ich].frequency = fsin->GetParameter(1);
	fBmodStatArray[ich].phase = fsin->GetParameter(2);
	fBmodStatArray[ich].offset =fsin->GetParameter(3);
	fBmodStatArray[ich].chisq =fsin->GetChisquare();
	fBmodStatArray[ich].ndf =(Double_t)fsin->GetNDF();
	fBmodStatArray[ich].ndata =npt;
      } // ich loop
      
      coilID = icoil+1;
      bmod->Fill();
    } // icoil loop
  } // icyc loop
  bmod->Write();
  output->Close();
  input->Close();
}

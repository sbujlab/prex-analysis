#include "utility.cc"
void DrawResidualSensitivity(Int_t run_number=4621){
  int slug = 71;
  map< pair<Int_t,Int_t>, Int_t > fC2BMap;
  fC2BMap =  LoadCyc2BurstMap(slug,run_number);
  TFile *input = TFile::Open(Form("$QW_ROOTFILES/prexPrompt_pass2_%d.000.root",run_number));
  if(input==NULL)
    return;
  TTree *evt_tree = (TTree*)input->Get("evt_bmw");
  if(evt_tree==NULL)
    return;
  // gStyle->SetOptFit(1);
  Int_t fCoilArray[]={1,3,5};
  Int_t nCoils= sizeof(fCoilArray)/sizeof(*fCoilArray);
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

  vector<TString> fChannelArray={"usl","usr"};
  
  // "bmod_trim1","bmod_trim2","bmod_trim3","bmod_trim4",
  // "bmod_trim5","bmod_trim6","bmod_trim7"};
  Int_t nch = fChannelArray.size();
  typedef struct {Double_t amplitude, phase, offset, period, chisq,ndf,ndata;} BMODSTAT;
  TString leaflist = "amplitdue/D:phase:offset:period:chisq:ndf:ndata";
  vector<BMODSTAT> fBmodStatArray(nch);
  for(int ich = 0; ich<nch;ich++)
    bmod->Branch(fChannelArray[ich],&(fBmodStatArray[ich]),leaflist);

  TF1 *fsin = new TF1("fsin","[0]*TMath::Sin(2*TMath::Pi()*x/[1] +[2])+[3]",-1e5,1e5);
  fsin->SetLineColor(kBlue);
  fsin->SetLineStyle(2);

  Double_t par_init[] =  {0.01,1200,0.1,0};
  fsin->SetParameters(par_init);
  fsin->SetParName(0,"Amplitude");
  fsin->SetParName(1,"Period");
  fsin->SetParName(2,"Phase");
  fsin->SetParName(3,"Offset");
  // fsin->SetParLimits(0,0.0,10);
  // fsin->SetParLimits(2,0.0,2*TMath::Pi());

  TCanvas *c1 = new TCanvas("c1","c1",800,400);
  Int_t plot_counts=0;
  c1->cd();
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.15);
  c1->Print("residual.pdf[");
  map< pair<Int_t,Int_t>, vector<Double_t> > fSlopeMap = GetLagrangeSlope(slug);
  for(int icyc=0;icyc<nCycle;icyc++){
    cycID = cycID_first+icyc; 
    Int_t kBurstCounter = -1;

    for(int icoil=0;icoil<nCoils;icoil++){
      int npt = evt_tree->Draw("bmod_ramp",
			       Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0",
				    fCoilArray[icoil],cycID_first+icyc),"goff");

      if( fC2BMap.find( make_pair(cycID,fCoilArray[icoil])) != fC2BMap.end())
	kBurstCounter = fC2BMap[ make_pair( cycID, fCoilArray[icoil])  ];
      vector<Double_t> fSlopes = fSlopeMap[ make_pair(run_number,kBurstCounter) ] ;

      if(npt==0)
	continue;
      
      TH1D *h1dptr = (TH1D*)gDirectory->FindObject("htemp");
      double period = h1dptr->GetBinCenter(h1dptr->FindLastBinAbove(0));
      Int_t phase_max =(Int_t)(period+200);
      Int_t bin_size = 10;
      Int_t max_bin = phase_max/bin_size;
      TH1D hramp("hramp","",max_bin,-0.5,phase_max-0.5);
      evt_tree->Draw("bmod_ramp>>hramp",
		     Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0",
			  fCoilArray[icoil],cycID_first+icyc),"");
      vector<Double_t> fRampSteps;
      for(int ibin=1;ibin<=max_bin;ibin++){
	if(hramp.GetBinContent(ibin)<=5)
	  continue;
	if(hramp.GetBinContent(ibin)>hramp.GetBinContent(ibin-1) &&
	   hramp.GetBinContent(ibin)>hramp.GetBinContent(ibin+1) )
	  fRampSteps.push_back( hramp.GetBinCenter(ibin));
      }
      Int_t nSteps = fRampSteps.size();
      if(nSteps==0)
	continue;
      
      period = *(fRampSteps.end()-1);
	
      for(int ich=0;ich<nch;ich++){
	if(evt_tree->GetBranch(fChannelArray[ich])==NULL){
	  fBmodStatArray[ich].amplitude = 0;
	  fBmodStatArray[ich].period = 0;
	  fBmodStatArray[ich].phase = 0;
	  fBmodStatArray[ich].offset = 0;
	  fBmodStatArray[ich].chisq =0;
	  fBmodStatArray[ich].ndf = 0;
	  fBmodStatArray[ich].ndata = 0;
	  continue;
	}
	  
	fsin->SetParameters(par_init);

	vector<Double_t> fMean;
	vector<Double_t> fError;
	Double_t* fRamp_ptr = new Double_t[nSteps];

	cout << "---" << endl;
	double fNormalization = 1;
	TString correction = "";
	vector<TString> bpm_name = {"11X","11Y","12X","12Y","16X","16Y",
				    "1X","1Y","4aX","4aY","4eX","4eY"};
	if(!fChannelArray[ich].Contains("bpm")){
	  evt_tree->Draw(Form("%s",fChannelArray[ich].Data()),
			 Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 &&bmod_ramp>0 "
			      ,fCoilArray[icoil],cycID_first+icyc),"goff"); 
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  fNormalization = h1dptr->GetMean();
	}
	for(int ibpm=0;ibpm<12;ibpm++){
	  evt_tree->Draw(Form("bpm%s",bpm_name[ibpm].Data()),
			 Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 &&bmod_ramp>0 "
			      ,fCoilArray[icoil],cycID_first+icyc),"goff"); 
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");

	  correction += Form("-(%f*(bpm%s-%f))",
			     fSlopes[ich*12+ibpm],bpm_name[ibpm].Data(),h1dptr->GetMean());
	}

	for(int ir=0;ir<nSteps;ir++){
	  cout << fRampSteps[ir] << endl;
	  evt_tree->Draw(Form("%s/%f %s",fChannelArray[ich].Data(),fNormalization,correction.Data()),
			 Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 && (bmod_ramp>=%f && bmod_ramp<=%f) ",fCoilArray[icoil],cycID_first+icyc,fRampSteps[ir]-30, fRampSteps[ir]+30),"goff"); 
	  
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  fMean.push_back(h1dptr->GetMean());
	  fError.push_back(h1dptr->GetMeanError()); // FIXME later
	  evt_tree->Draw("bmod_ramp",
			 Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 && (bmod_ramp>=%f && bmod_ramp<=%f) ",fCoilArray[icoil],cycID_first+icyc,fRampSteps[ir]-30, fRampSteps[ir]+30),"goff"); 
	  
	  h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	  fRamp_ptr[ir] = h1dptr->GetMean();
	}
	double ymax = fMean[0];
	double ymin = fMean[0];
	double peak_phase ;
	for(int ir=1;ir<nSteps;ir++){
	  if(fMean[ir]>ymax){
	    ymax = fMean[ir];
	    peak_phase = fRampSteps[ir];
	  }
	  if(fMean[ir]<ymin)
	    ymin = fMean[ir];

	}
	// evt_tree->Draw(Form("%s",fChannelArray[ich].Data()),
	// 	       Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 ",
	// 		    fCoilArray[icoil],cycID_first+icyc),"goff");

	// h1dptr = (TH1D*)gDirectory->FindObject("htemp");
	fsin->SetParameter(3,(ymax+ymin)*0.5);
	double phase_init = TMath::Pi()/2.0 - peak_phase/period*2*TMath::Pi();
	if(phase_init<0)
	  phase_init += 2*TMath::Pi();
	fsin->SetParameter(2,phase_init);
	fsin->SetParameter(1,period/4.0);
	fsin->SetParameter(0,(ymax-ymin)*0.5);
	Double_t* fMean_ptr = new Double_t[nSteps];
	Double_t* fError_ptr = new Double_t[nSteps];
	for(int is=0;is<nSteps;is++){
	  fMean_ptr[is] = fMean[is];
	  fError_ptr[is] = fError[is];
	}
	TMultiGraph *mgall = new TMultiGraph();
	TGraphErrors mger(nSteps,fRamp_ptr,fMean_ptr,0,fError_ptr);
	mger.SetMarkerStyle(20);
	// mger.SetMarkerSize(1.2);
	mger.SetMarkerColor(kBlue);
	mger.SetLineColor(kBlue);
	mger.SetLineWidth(2);

	evt_tree->Draw(Form("%s/%f %s:bmod_ramp",
			    fChannelArray[ich].Data(),fNormalization,correction.Data()),
		       Form("bmwobj==%d && bmwcycnum==%d && (ErrorFlag&0xda7e6bff)==0 ",
			    fCoilArray[icoil],cycID_first+icyc),"goff");
	TGraph *gbuff = new TGraph(evt_tree->GetSelectedRows(),evt_tree->GetV2(),evt_tree->GetV1());
	gbuff->SetMarkerStyle(7);
	gbuff->SetMarkerSize(1);
	gbuff->SetMarkerColor(11);
	mgall->Add(gbuff,"p");
	mgall->Add(&mger,"P");
	TString unit = " Fractional Yield   ";
	if(fChannelArray[ich].Contains("bpm"))
	  unit = " Position(um)";
	mgall->SetTitle(Form("Residual Sensitivity %s vs coil %d at cycle %d;Modulation Phase; %s", 
			     fChannelArray[ich].Data(),fCoilArray[icoil],cycID,unit.Data()));
	mgall->GetYaxis()->SetLabelSize(0.05);
	mgall->GetXaxis()->SetLabelSize(0.05);
	mgall->GetYaxis()->SetTitleSize(0.05);
	mgall->GetXaxis()->SetTitleSize(0.05);
	mgall->Draw("A");
	mger.Fit("fsin","Q");
	fBmodStatArray[ich].amplitude = fsin->GetParameter(0);
	fBmodStatArray[ich].period = fsin->GetParameter(1);
	fBmodStatArray[ich].phase = fsin->GetParameter(2);
	fBmodStatArray[ich].offset =fsin->GetParameter(3);
	fBmodStatArray[ich].chisq =fsin->GetChisquare();
	fBmodStatArray[ich].ndf =(Double_t)fsin->GetNDF();
	fBmodStatArray[ich].ndata =npt;
			
	c1->Print("residual.pdf");
      } // ich loop
      
      coilID = fCoilArray[icoil];
      bmod->Fill();
    } // icoil loop
  } // icyc loop
  c1->Print("residual.pdf]");
  bmod->Write();
  output->Close();
  input->Close();
}

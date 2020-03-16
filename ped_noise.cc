#include "src/TaAccumulator.cc"
#include "src/TaEventRing.cc"
void ped_noise(Int_t slug);
void ped_noise(){
  for(int i=101;i<=167;i++)
    ped_noise(i);
}
void ped_noise(Int_t slug){
  ROOT::EnableImplicitMT(8);
  TFile* slug_file = TFile::Open(Form("ped_rootfiles/slug%d.root",slug),"RECREATE");
  TTree *ped_tree = new TTree("ped","ped noise");
  ped_tree->SetMarkerStyle(20);
  vector< TString > fDetList={"sam1","sam2","sam3","sam4","sam5","sam6","sam7","sam8"};
  Double_t *fSamYield_ptr = new Double_t[8];
  Double_t *fSamNoise_ptr = new Double_t[8];
  for(int isam=0;isam<8;isam++){
    ped_tree->Branch(Form("sam%d_yield",isam+1),&fSamYield_ptr[isam]);
    ped_tree->Branch(Form("sam%d_noise",isam+1),&fSamNoise_ptr[isam]);
  }
  TString runlist_filename = Form("./prex-runlist/simple_list/slug%d.list",slug);
  FILE *runlist = fopen(runlist_filename.Data(),"r");
  Int_t run_number;
  ped_tree->Branch("run",&run_number);
  while(!feof(runlist)){
    fscanf(runlist,"%d\n",&run_number);
    TStopwatch *tsw = new TStopwatch();
    TString filename = Form("prexPrompt_pass1_%d.000.root",run_number);
    TFile* input = TFile::Open("$QW_ROOTFILEs/"+filename);
    if(input==NULL)
      continue;
    TTreeReader myReader("mul",input);
    TTreeReader myReader_evt("evt",input);
    TTreeReaderValue<Double_t> fBCMValue(myReader,"yield_bcm_an_us.hw_sum");
    TTreeReaderValue<Double_t> fErrorFlag(myReader,"ErrorFlag");
    vector< TTreeReaderValue<Double_t> > fAsym;
    vector< TTreeReaderValue<Double_t> > fYield;
    Int_t nDet=fDetList.size();
    for(int idet=0;idet<nDet;idet++){
      TString chname=Form("yield_%s.hw_sum",
			  fDetList[idet].Data());
      TTreeReaderValue<Double_t> fbuff(myReader,chname);

      TString asym_name=Form("asym_%s.hw_sum",
			     fDetList[idet].Data());
      TTreeReaderValue<Double_t> fasym_buff(myReader,asym_name);

      fYield.push_back(fbuff);
      fAsym.push_back(fasym_buff);
    }

    Int_t pat_counter=0;
    TaEventRing* fEventRing = new TaEventRing();
    vector<Double_t> fDetData(nDet);
    TaAccumulator aProtoAccumualtor;
    vector<TaAccumulator> fGoodAvg(nDet,aProtoAccumualtor);
    vector<TaAccumulator> fPedestalAvg(nDet,aProtoAccumualtor);
    vector<TH1D*>  fh1dArray(nDet);
    gStyle->SetStatH(0.3);
    gStyle->SetStatW(0.3);
    for(int idet=0;idet<nDet;idet++){
      if(idet==1)
	fh1dArray[idet] = new TH1D("",Form("SAM%d",idet+1),100,-800,800);
      else
	fh1dArray[idet] = new TH1D("",Form("SAM%d",idet+1),100,-300,300);
    }

    // ================ Determine beam off readout in uA
    TTreeReaderValue<Double_t> fBCMraw(myReader_evt,"bcm_an_us.hw_sum_raw");
    TTreeReaderValue<Double_t> fBCMnsamp(myReader_evt,"bcm_an_us.num_samples");
    TTreeReaderValue<Double_t> fBCMuA(myReader_evt,"bcm_an_us.hw_sum");
    double bcm_adc[2];
    double bcm_uA[2];
    for(int i=0;i<20;i++) myReader_evt.Next();
    myReader_evt.Next();
    bcm_adc[0]=(*fBCMraw)/(*fBCMnsamp);
    bcm_uA[0]=(*fBCMuA);
    myReader_evt.Next();
    bcm_adc[1]=(*fBCMraw)/(*fBCMnsamp);
    bcm_uA[1]=(*fBCMuA);
  
    double scale = (bcm_uA[1]-bcm_uA[0])/(bcm_adc[1]- bcm_adc[0]);
    double pedestal = bcm_adc[1]-bcm_uA[1]/scale ;  
    double beam_off_adc = -346.7;
    double beam_off_uA = (beam_off_adc-pedestal)*scale;
    cout << "beam off uA " << beam_off_uA << endl;
    fEventRing->SetBeamOffLimit(beam_off_uA+0.10);
    // ================ 

    while(myReader.Next()){
      fEventRing->PushBeamCurrent(*fBCMValue);
      for(int idet=0;idet<nDet;idet++){
	fDetData[idet]= (*fYield[idet])*(*fAsym[idet]);
      }
      fEventRing->PushDetector(fDetData);
      if( fEventRing->isReady() ){
	// cout << " Ring is Ready " << endl;
	// cout << pat_counter << ":" << *fBCMValue << endl;
	vector<Double_t> fPopback ;
	if(fEventRing->Pop(fPopback)){
	  for(int idet=0;idet<nDet;idet++){
	    fPedestalAvg[idet].Update(fPopback[idet]);
	    fh1dArray[idet]->Fill(fPopback[idet]*1e6);
	  }
	}
      }
      if(*fErrorFlag==0){
	for(int idet=0;idet<nDet;idet++){
	  fGoodAvg[idet].Update((*fYield[idet])*(*fBCMValue));
	}
      }
      pat_counter++;
    }
    if(fGoodAvg[0].GetN()==0)
      continue;
    if(fPedestalAvg[0].GetN()==0)
      continue;

    for(int idet=0;idet<nDet;idet++){
      fSamYield_ptr[idet] = fGoodAvg[idet].GetMean1();
      // fSamNoise_ptr[idet] = fPedestalAvg[idet].GetRMS()/fGoodAvg[idet].GetMean1()*1e6; // in ppm;
      fSamNoise_ptr[idet] = fh1dArray[idet]->GetRMS()/fGoodAvg[idet].GetMean1(); // in ppm;

    }
    ped_tree->Fill();
    cout << " **** run " << run_number << endl;
    cout << " ======== Beam-ON Yield ========== " << endl;
    for(int idet=0;idet<nDet;idet++){
      cout << fGoodAvg[idet].GetMean1() << " Volt ";
      // cout << fPedestalAvg[idet].GetRMS()*1e6 << " uV ";
      // cout << fPedestalAvg[idet].GetRMS()/fGoodAvg[idet].GetMean1()*1e6 << " ppm "
      // 	   << endl;
      cout << fh1dArray[idet]->GetRMS() << " uV ";
      cout << fh1dArray[idet]->GetRMS()/fGoodAvg[idet].GetMean1() << " ppm "
	   << endl;;

    }
    TCanvas *c1 = new TCanvas("c1","c1",1000,600);
    c1->Divide(4,2);
    for(int idet=0;idet<nDet;idet++){
      c1->cd(idet+1);
      fh1dArray[idet]->Draw();
    }
    cout << pat_counter << " patterns in " ;
    tsw->Print();
  }
  slug_file->cd();
  ped_tree->Write();
  slug_file->Close();
}

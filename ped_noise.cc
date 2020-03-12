#include "src/TaAccumulator.cc"
// #include "src/TaEventRing.cc"
void ped_noise(Int_t run_number){
  TStopwatch *tsw = new TStopwatch();
  TString filename = Form("prexPrompt_pass1_%d.000.root",run_number);
  TFile* input = TFile::Open("$QW_ROOTFILEs/"+filename);
  TTreeReader myReader("mul",input);
  TTreeReaderValue<Double_t> fBCMValue(myReader,"yield_bcm_an_us.hw_sum");
  TTreeReaderValue<Double_t> fErrorFlag(myReader,"ErrorFlag");
  vector< TTreeReaderValue<Double_t> > fAsym;
  vector< TTreeReaderValue<Double_t> > fYield;
  vector< TString > fDetList={"sam1","sam2","sam3","sam4","sam5","sam6","sam7","sam8"};
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
  TaEventRing fEventRing;
  vector<TaAccumulator> fGoodAvg(nDet);
  vector<TaAccumulator> fPedestalAvg(nDet);
  while(myReader.Next()){
    if(*fBCMValue<5.0){
      fEventRing.PushBeamCurrent(*fBCMValue);
      fEventRing.Push(fYield);
      if( fEventRing.isReady() ){
	vector<Double_t> fPopback = fEventRing.Pop();
	for(int idet=0;idet<nDet;idet++){
	  fPedestalAvg[idet].Update(fPopback(idet));
	}
      }else
	fEventRing.Pop();

    }else if(*fErrorFlag==0){
      for(int idet=0;idet<nDet;idet++){
	fRunningAvg[idet].Update((*fYield[idet])*(*fBCMValue));
      }
    }
    pat_counter++;
  }

  // if(input==NULL)
  //   break;
  // TTree *mul_tree = input->Get("mul");
  // TEventList *elist = new TEventList("elist");
  // mul_tree->Draw(">>elist","yield_bcm_an_us<5.0");
  cout << " ======== Beam-ON Yield ========== " << endl;
  for(int idet=0;idet<nDet;idet++){
    cout << fRunningAvg[idet].GetMean1() << endl;
  }
  
  cout << pat_counter << " patterns in " ;
  tsw->Print();
}

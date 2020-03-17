TGraph* GraphDitheringSlope2D(Int_t slug,
			      TString ch1,TString ch2){

  TFile *slope_file= TFile::Open(Form("./averaged_slopes/slug%d_dit_slope_cyclewise_average.root",slug));
  TTree *dit_tree = (TTree*)slope_file->Get("dit");

  Double_t slope1;
  Double_t slope2;
  Int_t run;
  Int_t range;
  Double_t prev_range=-1;
  vector<Double_t> fSlope1 ;
  vector<Double_t> fSlope2 ;
  dit_tree->SetBranchAddress(ch1,&slope1);
  dit_tree->SetBranchAddress(ch2,&slope2);
  dit_tree->SetBranchAddress("run",&run);
  dit_tree->SetBranchAddress("range",&range);
  Int_t nevt = dit_tree->GetEntries();
  
  for(int ievt=0;ievt<nevt;ievt++){
    dit_tree->GetEntry(ievt);
    if(range!=prev_range){
      prev_range = range;
      fSlope1.push_back(slope1*1e3);
      fSlope2.push_back(slope2*1e3);
    }
  }// end of event loop
  Int_t nrange  = fSlope1.size();
  Double_t *fy = new Double_t[nrange];
  Double_t *fx = new Double_t[nrange];
  for(int i=0;i<nrange;i++){
    fy[i]=fSlope1[i];
    fx[i]=fSlope2[i];
  }
    
  TGraph *fret = new TGraph(nrange,fx,fy);
  fret->SetMarkerStyle(29);
  fret->SetMarkerSize(3);
  fret->SetMarkerColor(kBlue);
  slope_file->Close();
  return fret;
}

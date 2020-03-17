TMultiGraph* GraphDitheringSlope(Double_t *fRun,Int_t npt, TString chname){

  TFile *slope_file= TFile::Open("dit_averaged_slope.root");
  TTree *dit_tree = (TTree*)slope_file->Get("dit");
  Double_t *fSlopeVal  = new Double_t[npt];
  Double_t *fXcord  = new Double_t[npt];
  Double_t slope;
  Double_t run;
  Double_t range;
  Double_t prev_range=-1;
  vector<Int_t> start_pt;
  dit_tree->SetBranchAddress(chname,&slope);
  dit_tree->SetBranchAddress("run",&run);
  dit_tree->SetBranchAddress("range",&range);
  Int_t nevt = dit_tree->GetEntries();
  for(int ievt=0;ievt<nevt;ievt++){
    dit_tree->GetEntry(ievt);
    for(int ipt=0;ipt<npt;ipt++){
      fXcord[ipt]=ipt;
      if(run ==fRun[ipt] ) {
	fSlopeVal[ipt] = -slope*1e3; // old convention ,a minus sign
	if(range!=prev_range){
	  prev_range = range;
	  start_pt.push_back(ipt);
	}
      }
    } 
  }// end of event loop
  TMultiGraph *fmg = new TMultiGraph();
  Int_t nrange = start_pt.size();
  for(int i=0;i<nrange;i++){
    TGraph *g1;
    if(i!=nrange-1)
      g1 = new TGraph(start_pt[i+1]-start_pt[i],
		      fXcord+start_pt[i],
		      fSlopeVal+start_pt[i]);
    else
      g1 = new TGraph(npt-start_pt[i],
		      fXcord+start_pt[i],
		      fSlopeVal+start_pt[i]);
    
    g1->SetLineColor(kRed);
    g1->SetLineWidth(2);
    fmg->Add(g1,"l");
  }
  
  slope_file->Close();
  return fmg;
}

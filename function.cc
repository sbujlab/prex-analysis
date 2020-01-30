TString get_base(TString input){
  input.ReplaceAll("diff_","");
  input.ReplaceAll("asym_","");
  return input;
}

map<Int_t, vector<Double_t> > LoadSlopeMap(TTree *slope_tree, vector<TString> detlist, vector<TString> ivlist){
  
  map<Int_t, vector<Double_t> > fSlopeMap;
  Int_t ndet = detlist.size();
  Int_t nbpm = ivlist.size();
  vector<Double_t> fSlope(ndet*nbpm);
  Double_t fRun;
  slope_tree->SetBranchAddress("run",&fRun);

  for(int idet=0;idet<ndet;idet++){
    for(int ibpm=0;ibpm<nbpm;ibpm++){
      TString det_base = get_base(detlist[idet]);
      TString bpm_base = get_base(ivlist[ibpm]);
      TString chname = Form("%s_%s",det_base.Data(),bpm_base.Data());
      slope_tree->SetBranchAddress(chname,&fSlope[idet*nbpm+ibpm]);
    }
  }
  
  Int_t nevt = slope_tree->GetEntries();
  for(int i=0;i<nevt;i++){
    slope_tree->GetEntry(i);
    Int_t myRun = (Int_t)fRun;
    vector<Double_t> fMySlope;
    for(int idet=0;idet<ndet;idet++){
      for(int ibpm=0;ibpm<nbpm;ibpm++){
	Double_t slope_val  = fSlope[idet*nbpm+ibpm] ;
	if(detlist[idet].Contains("us_avg"))
	  slope_val = 0.5*(fSlope[ibpm]+fSlope[nbpm+ibpm]);
	else
	  fMySlope.push_back(slope_val);
	
	if(ivlist[ibpm].Contains("11X12X"))
	  fMySlope.push_back(0.4*slope_val);
      }
    }
    fSlopeMap[myRun]=fMySlope;
  }
  return fSlopeMap;
}



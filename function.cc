TString get_base(TString input){
  input.ReplaceAll("diff_","");
  input.ReplaceAll("asym_","");
  return input;
}

vector<Int_t> ParseRunList(TString list_name){
  vector<Int_t> fRet;
  FILE *runlist = fopen(list_name.Data(),"r");
  if(runlist==NULL)
    return fRet;
  
  while(!feof(runlist)){
    Int_t run_number=0;
    fscanf(runlist,"%d\n",&run_number);
    if(run_number!=0){
      fRet.push_back(run_number);
    }
  }
  fclose(runlist);
  return fRet;
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

	fMySlope.push_back(slope_val);
	
	if(ivlist[ibpm].Contains("11X12X"))
	  fMySlope.push_back(0.4*slope_val);
      }
    }
    fSlopeMap[myRun]=fMySlope;
  }
  return fSlopeMap;
}

map<Int_t, vector<Double_t> > LoadAvgSlopeMap(Int_t slug_number, vector<TString> detlist, vector<TString> ivlist){
  
  map<Int_t, vector<Double_t> > fSlopeMap;
  Int_t ndet = detlist.size();
  Int_t nbpm = ivlist.size();
  vector<Double_t> fSlope(ndet*nbpm);

  TFile *slope_rf;
  if(slug_number>=4)
    slope_rf = TFile::Open("./slopes/ditcoeff_slope_averaged.root");
  else
    slope_rf = TFile::Open("./slopes/ditcoeff_slope_averaged12X.root");
  if(slope_rf==NULL){
    cerr << " -- Dit Slope rootfile not found " << endl;
    return fSlopeMap;
  }

  TString runlist = Form("./prex-runlist/simple_list/slug%d.list",slug_number);
  vector<Int_t> fRunList = ParseRunList(runlist);

  TTree *slope_tree = (TTree*)slope_rf->Get("dit");
  Double_t fRun;
  slope_tree->SetBranchAddress("run",&fRun);
  for(int idet=0;idet<ndet;idet++){
    for(int ibpm=0;ibpm<nbpm;ibpm++){
      TString det_base = get_base(detlist[idet]);
      TString bpm_base = get_base(ivlist[ibpm]);
      bpm_base.ReplaceAll("bpm","");
      TString chname = Form("%s_%s",det_base.Data(),bpm_base.Data());
      slope_tree->SetBranchAddress(chname,&fSlope[idet*nbpm+ibpm]);
    }
  }

  Int_t nevt = slope_tree->GetEntries();
  for(int i=0;i<nevt;i++){
    slope_tree->GetEntry(i);
    Int_t myRun = (Int_t)fRun;
    if(find(fRunList.begin(),fRunList.end(),myRun)==fRunList.end())
      continue;
    vector<Double_t> fMySlope;
    for(int idet=0;idet<ndet;idet++){
      for(int ibpm=0;ibpm<nbpm;ibpm++){
	Double_t slope_val  = fSlope[idet*nbpm+ibpm] ;
	if(detlist[idet].Contains("us_avg"))
	  slope_val = 0.5*(fSlope[ibpm]+fSlope[nbpm+ibpm]);

	fMySlope.push_back(slope_val);
	
	if(ivlist[ibpm].Contains("11X12X"))
	  fMySlope.push_back(0.4*slope_val);
      }
    }
    fSlopeMap[myRun]=fMySlope;
  }
  slope_rf->Close();
  return fSlopeMap;

}


#include "device_list.h"
// For 120 Hz raw data
// nPattern=4 : (1+2)-(3+4) 30 Hz Equivalent
// nPattern=8 : (1+2+3+4)-(1+2+3+4) 15 Hz Equivalent

void GetPseudoAsym(TString filename, Int_t nPattern = 4){

  TFile *japanOutput = TFile::Open(filename);
  TTree* evt_tree = (TTree*)japanOutput->Get("evt");
  
  const Int_t ndet = sizeof(detector_channel)/sizeof(*detector_channel);
  const Int_t nmon = sizeof(monitor_channel)/sizeof(*monitor_channel);

  std::vector<TLeaf*> leaf_det;
  std::vector<TLeaf*> leaf_mon;
  for(int ich=0 ; ich<ndet;ich++){
    TBranch *branch_buff =(TBranch*)evt_tree->GetBranch(detector_channel[ich]);
    TLeaf *leaf_buff = branch_buff->GetLeaf("hw_sum");
    leaf_det.push_back(leaf_buff);
  }
  for(int ich=0 ; ich<nmon;ich++){
    TBranch *branch_buff =(TBranch*)evt_tree->GetBranch(monitor_channel[ich]);
    TLeaf *leaf_buff = branch_buff->GetLeaf("hw_sum");
    leaf_mon.push_back(leaf_buff);
  }
  TLeaf *leaf_ErrorFlag = evt_tree->GetLeaf("ErrorFlag");
  
  TString outputName;
  outputName = filename.ReplaceAll("prex","pseudo");
    
  TFile *peusdoOutput = TFile::Open(filename,"RECREATE");
  
  TTree *mul_tree = new TTree("mul","pseudo Helicity Tree");
  TTree *mulc_tree = new TTree("mulc","pseudo Helicity Combiner");
  const Int_t ncombo = sizeof(combiner_channel)/sizeof(*combiner_channel);

  Double_t kErrorFlag;
  Double_t asym_value[ndet];
  Double_t yield_value[ndet];
  Double_t diff_value[nmon];
  Double_t combo_value[ncombo];

  mul_tree->Branch("ErrorFlag",&kErrorFlag,"ErrorFlag/D");
  for(int idet=0;idet<ndet;idet++)
    mul_tree->Branch("asym_"+detector_channel[idet],
		     &(asym_value[idet]),"hw_sum/D");
  
  for(int imon=0;imon<nmon;imon++)
    mul_tree->Branch("diff_"+monitor_channel[imon],
		     &(diff_value[imon]),"hw_sum/D");

  for(int icom=0;icom<ncombo;icom++)
    mulc_tree->Branch(combiner_channel[icom],
		     &(combo_value[icom]),"hw_sum/D");
  
  Int_t nEntries = evt_tree->GetEntries();
  nEntries = nEntries - (nEntries%nPattern);
  Int_t pattern_counter=0;
  for(int ievt=0; ievt<nEntries; ){
    // Zero Initialization
    kErrorFlag = 0;
    for(int idet=0;idet<ndet;idet++){
      asym_value[idet] = 0.0;
      yield_value[idet] = 0.0;
    }
    for(int imon=0;imon<nmon;imon++)
      diff_value[imon] =0.0;

    for(int icom=0;icom<ncombo;icom++)
      combo_value[icom] =0.0;

    // Filling Pattern
    for(int ipattern=0;ipattern<nPattern;ipattern++){
      evt_tree->GetEntry(ievt++); // Not sure if it starts with 1 or zero
      kErrorFlag = (UInt_t)(kErrorFlag)|(UInt_t)(leaf_ErrorFlag->GetValue());
      Int_t polarity = (ipattern < (nPattern/2) ? 1: -1);
      
      for(int idet=0;idet<ndet;idet++){
	asym_value[idet]+= (polarity*(leaf_det[idet]->GetValue()));
	yield_value[idet]+= (leaf_det[idet]->GetValue());
      } // end of detector loop

      for(int imon=0;imon<nmon;imon++){
	diff_value[imon]+= (polarity*(leaf_mon[imon]->GetValue()));
      } // end of monitor loop

    } // end of pattern loop

    // Calculate diff and asym  then fill the tree
    for(int idet=0;idet<ndet;idet++)
      asym_value[idet] = asym_value[idet]/yield_value[idet];
    for(int imon=0;imon<nmon;imon++)
      diff_value[imon] = diff_value[imon]/nPattern;

    // Combination of detectors
    for(int icom=0;icom<ncombo;icom++){
      
      double weight_a = weight[icom][0];
      double weight_b = weight[icom][1];
      int det_a = det_pair[icom][0];
      int det_b = det_pair[icom][1];
      combo_value[icom] = weight_a*asym_value[det_a]+ weight_b*asym_value[det_b];
    }

    mul_tree->Fill();    
    mulc_tree->Fill();
    pattern_counter++;
    if(pattern_counter%1000==0)
      cout << pattern_counter << endl;
  } // end of raw entry look
  
  mul_tree->Write();
  mulc_tree->Write();

  japanOutput->Close();
  peusdoOutput->Close();
  
}

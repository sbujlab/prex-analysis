#include "device_list.h"
// For 120 Hz raw data
// over-integrating block0 and block1 to get equivalent 240 Hz Data

// FIXME : though I may have a better idea to do this, 
// Now I just pick up a way that is easy and simple

void GetPseudo240Hz(TString filename){
  TFile *japanOutput = TFile::Open(filename);
  TTree* evt_tree = (TTree*)japanOutput->Get("evt");

  Int_t nPattern = 4; // Needs 4 evt to build octet
  Int_t helicity[8] = { 1, -1, -1, 1,
			-1,  1,  1,-1}; // pseudo helicity, to cancel 60 Hz

  const Int_t ndet = sizeof(detector_channel)/sizeof(*detector_channel);
  const Int_t nmon = sizeof(monitor_channel)/sizeof(*monitor_channel);

  std::vector< vector<TLeaf*> > leaf_det; // [idet][iblock]
  std::vector< vector<TLeaf*> > leaf_mon; // [imon][iblock]
  for(int ich=0 ; ich<ndet;ich++){
    TBranch *branch_buff =(TBranch*)evt_tree->GetBranch(detector_channel[ich]);
    vector<TLeaf*> vec_buff;
    for(int iblock=0;iblock<4;iblock++){
      TLeaf *leaf_buff = branch_buff->GetLeaf(Form("block%d",iblock));
      vec_buff.push_back(leaf_buff);
    }
    leaf_det.push_back(vec_buff);
  }
  for(int ich=0 ; ich<nmon;ich++){
    TBranch *branch_buff =(TBranch*)evt_tree->GetBranch(monitor_channel[ich]);
    vector<TLeaf*> vec_buff;
    for(int iblock=0;iblock<4;iblock++){
      TLeaf *leaf_buff = branch_buff->GetLeaf(Form("block%d",iblock));
      vec_buff.push_back(leaf_buff);
    }
    leaf_mon.push_back(vec_buff);
  }
  TLeaf *leaf_ErrorFlag = evt_tree->GetLeaf("ErrorFlag");
  
  TString outputName;
  outputName = filename.ReplaceAll("prex","pseudo240Hz");
    
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

    Int_t index_hel = 0;
    // Filling Pattern
    for(int ipattern=0;ipattern<nPattern;ipattern++){
      evt_tree->GetEntry(ievt++); // Not sure if it starts with 1 or zero
      kErrorFlag = (UInt_t)(kErrorFlag)|(UInt_t)(leaf_ErrorFlag->GetValue());
      for(int ipair=0;ipair<2;ipair++){
	Int_t polarity = helicity[index_hel];
	index_hel ++;
	for(int idet=0;idet<ndet;idet++){
	  double pair_sum = (leaf_det[idet][ipair*2])->GetValue();
	  pair_sum += (leaf_det[idet][ipair*2+1])->GetValue();
	  asym_value[idet]+= (polarity*pair_sum);
	  yield_value[idet]+= pair_sum;
	} // end of detector loop
	for(int imon=0;imon<nmon;imon++){
	  double pair_sum = (leaf_mon[imon][ipair*2])->GetValue();
	  pair_sum += (leaf_mon[imon][ipair*2+1])->GetValue();
	  diff_value[imon]+= (polarity*pair_sum);
	} // end of monitor loop
      } // end of block-pair loop
    } // end of pattern loop

    // Calculate diff and asym  then fill the tree
    for(int idet=0;idet<ndet;idet++)
      asym_value[idet] = asym_value[idet]/yield_value[idet];
    for(int imon=0;imon<nmon;imon++)
      diff_value[imon] = diff_value[imon]/(nPattern*2);

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

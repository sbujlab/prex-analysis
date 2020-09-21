#include "LoadRunInfoMap.cc"
#include "TaRunInfo.cc"
// Find the lowest 7 correction from eigenvectors
// Use it as a cross-check

void MinirunTruncBySlug(){
  TString output_filename = "prex_minidiff_trunc_lagrall.root";
  TFile *output_file = TFile::Open(output_filename,"RECREATE");
  TTree *output_tree = new TTree("slug","");
  Int_t fSlug;
  Int_t fSign,fArm_flag;
  output_tree->Branch("slug",&fSlug,"slug/I");
  output_tree->Branch("sign",&fSign,"sign/I");
  output_tree->Branch("arm_flag",&fArm_flag);
  vector<Double_t> delta(4,0);
  vector<Double_t> distance(4,0);
  vector<Double_t> delta_error(4,0);
  TString det_array[]={"us_avg","usr","usl","us_dd"};
  Int_t nDet=4;
  for(int idet=0;idet<nDet;idet++){
    output_tree->Branch("delta_"+det_array[idet],&delta[idet]);
    output_tree->Branch("dist_"+det_array[idet],&distance[idet]);
    output_tree->Branch("delta_"+det_array[idet]+"_err",&delta_error[idet]);
  }
  for(int islug=1;islug<=94;islug++){
    // Load Input Trees
    TFile *burst_file = TFile::Open(Form("./treeMergeOutput/MergedSum_slug%d.root",islug));
    TTree *burst_tree  = (TTree*)burst_file->Get("burst");
    TFile *input_file1 = TFile::Open(Form("./treeMergeOutput/MergedLagrange_slug%d.root",islug));
    TFile *input_file2 = TFile::Open(Form("./treeMergeOutput/dit_eigslopes_allbpm_slug%d.root",islug));
    TTree *eig_tree = (TTree*)input_file1->Get("mini_eigall");
    TTree *slope_tree = (TTree*)input_file2->Get("slope");
    burst_tree->AddFriend(eig_tree);
    burst_tree->AddFriend(slope_tree);
    Int_t run_number;
    Int_t nIV;
    if(islug>=3)
      nIV = 12;
    else
      nIV = 10;
    vector<Double_t> fSlopes(4*nIV);
    vector<Double_t> fDiff_evMon_rms(nIV);
    vector<Double_t> fDiff_evMon_mean(nIV);
    
    burst_tree->SetBranchAddress("run",&run_number);
    for(int iiv=0;iiv<nIV;iiv++){
      TBranch *branch_ptr = eig_tree->GetBranch(Form("diff_evMon%d",iiv));
      branch_ptr->GetLeaf("mean")->SetAddress( &fDiff_evMon_mean[iiv]);
      branch_ptr->GetLeaf("rms")->SetAddress( &fDiff_evMon_rms[iiv]);
    }
    for(int idv=0;idv<4;idv++)
      for(int iiv=0;iiv<nIV;iiv++)
	slope_tree->SetBranchAddress(Form("%s_evMon%d",det_array[idv].Data(),iiv),
				     &fSlopes[idv*nIV+iiv]);
    
    for(int arm_switch=0;arm_switch<3;arm_switch++){
      map< Int_t, TaRunInfo > fRunInfoMap = LoadRunInfoMapBySlug(islug,arm_switch);
      if(fRunInfoMap.size()==0)
	continue;
      Int_t nEntries = burst_tree->GetEntries();

      
      for(int idet=0;idet<4;idet++){
	delta[idet]=0;
	distance[idet]=0;
	delta_error[idet]=0;
      }
      Int_t nCounts = 0;
      for(int ievt=0;ievt<nEntries;ievt++){
	burst_tree->GetEntry(ievt);
	if( fRunInfoMap.find(run_number) == fRunInfoMap.end())
	  continue;
	
	for(int idet=0;idet<nDet;idet++){
	  vector<Double_t> mySlope(nIV);
	  vector<Double_t> fProduct(nIV);
	  vector<Int_t> fRank(nIV);
	  for(int iiv=0;iiv<nIV;iiv++){
	    mySlope[iiv] = fSlopes[ idet*nIV +iiv];
	    fProduct[iiv] = fabs(mySlope[iiv]*fDiff_evMon_rms[iiv]);
	    fRank[iiv] = iiv;
	  }
	  // Ranking 
	  for(int i=1;i<nIV;i++){
	    for(int j=i-1;j>=0;j--){
	      if(fProduct[j]< fProduct[j+1] ){
		Double_t swap = fProduct[j];
		fProduct[j] = fProduct[j+1];
		fProduct[j+1] = swap;
		Int_t idx_swap = fRank[j];
		fRank[j] = fRank[j+1];
		fRank[j+1] = idx_swap;
	      }
	    }
	  }
	  double correction_sum =0;
	  for(int iiv=5;iiv<nIV;iiv++){
	    correction_sum += mySlope[ fRank[iiv] ]  * fDiff_evMon_mean[ fRank[iiv] ] ;
	  }
	  delta[idet] += correction_sum;
	  distance[idet] += TMath::Power(correction_sum,2);
	} // loop over detector channels
	nCounts++;
      } // end of minirun loop
      for(int idet=0;idet<4;idet++){
	delta[idet]  = delta[idet]/nCounts*1e9;
	distance[idet] = TMath::Sqrt(distance[idet]/nCounts)*1e9;
	delta_error[idet] = distance[idet]/sqrt(nCounts);
      }
      // Fill output tree with slug info
      fSlug = islug;
      fArm_flag = arm_switch;
      fSign = (*(fRunInfoMap.begin())).second.GetSign();
      output_tree->Fill();
    } // end of arm switch loop;

    input_file1->Close();
    input_file2->Close();
    burst_file->Close();
  } // end of loop over slugs

  output_file->cd();
  output_tree->Write();
  output_file->Close();
}

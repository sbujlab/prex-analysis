#include "utility.cc"
#include "TaAccumulator.cc"

void cov(Int_t slug_number=0){
  
  TString output_filename = "output.root";
  TFile* output_file = TFile::Open(output_filename,"RECREATE");
  TTree *cov_tree = new TTree("covv","Covariance / Variance Tree");
  
  vector<TString> dv_list ={"asym_usl","asym_usr","asym_us_avg"};
  vector<TString> iv_list ={"diff_bpm4aX","diff_bpm4eX","diff_bpm4aY","diff_bpm4eY",
  			    "diff_bpm11X","diff_bpm12X"};

  vector<TString> all_list = dv_list;
  all_list.insert(all_list.end(),iv_list.begin(),iv_list.end());
  
  Int_t ndv = dv_list.size();
  Int_t niv = iv_list.size();
  
  vector<STAT> fIVCovStat(niv*niv);
  vector<STAT> fDVIVCovStat(ndv*niv);
  vector<STAT> fVarStat(ndv+niv);

  Int_t fRun;
  Double_t fPatterns;
  for(int i=0;i<niv;i++){
    for(int j=0;j<niv;j++){
      TString base1 = get_basename(iv_list[i]);
      TString base2 = get_basename(iv_list[j]);
      cov_tree->Branch(base1+"_"+base2,&fIVCovStat[i*niv+j]);
    }
  }
  for(int i=0;i<ndv;i++){
    for(int j=0;j<niv;j++){
      TString base1 = get_basename(dv_list[i]);
      TString base2 = get_basename(iv_list[j]);
      cov_tree->Branch(base1+"_"+base2,&fDVIVCovStat[i*niv+j]);
    }
  }
  for(int i=0;i<niv+ndv;i++)
    cov_tree->Branch(all_list[i],&fVarStat[i]);

  
  cov_tree->Branch("run",&fRun,"run/I");
  cov_tree->Branch("nsamples",&fPatterns);
  vector<TaAccumulator> fIVCovAcccumulator(niv*niv);
  vector<TaAccumulator> fDVIVCovAcccumulator(ndv*niv);
  vector<TaAccumulator> fVarAcccumulator(ndv+niv);
  
  vector<Int_t> fRunlist=ParseRunList("test.list");
  auto iter = fRunlist.begin();
  while(iter!=fRunlist.end()){
    fRun = (*iter);
    TString input_path = "/home/yetao/workarea/japanOutput/";
    TString input_filename = Form("prexPrompt_pass1_%d.000.root",fRun);
    TFile *input_file=TFile::Open(input_path+input_filename);
    if(input_file==NULL){
      iter++;
      cout << " --Warning " << input_filename << " is not found " << endl;
      continue;
    }
    TTree *mul_tree = (TTree*)input_file->Get("mul");
    TTree *mulc_tree = (TTree*)input_file->Get("mulc");
    mul_tree->AddFriend(mulc_tree);

    vector<Double_t> fdv_val(ndv);
    vector<Double_t> fiv_val(niv);

    for(int idv=0;idv<ndv;idv++){
      TBranch *branch_ptr = mul_tree->GetBranch(dv_list[idv]);
      if(branch_ptr!=NULL)
	branch_ptr->GetLeaf("hw_sum")->SetAddress(&fdv_val[idv]);
      else
	fdv_val[idv] = 1e6;
    }
    for(int iiv=0;iiv<niv;iiv++){
      TBranch *branch_ptr = mul_tree->GetBranch(iv_list[iiv]);
      if(branch_ptr!=NULL)
	branch_ptr->GetLeaf("hw_sum")->SetAddress(&fiv_val[iiv]);
      else
	fiv_val[iiv] = 1e6;
    }
  
    TEventList *elist = new TEventList("elist");
    Int_t nevt = mul_tree->Draw(">>+elist","ErrorFlag==0");
    cout << " -- " << nevt << " good patterns. " << endl;
    fPatterns = nevt;
    for(int ievt=0;ievt<nevt;ievt++){
      Int_t index = elist->GetEntry(ievt);
      mul_tree->GetEntry(ievt);
      for(int i=0;i<niv;i++){
      	fVarAcccumulator[ndv+i].Update(fiv_val[i]);
      	for(int j=0;j<niv;j++){
      	  fIVCovAcccumulator[i*niv+j].Update(fiv_val[i],fiv_val[j]);
      	}
      }
    
      for(int i=0;i<ndv;i++){
      	fVarAcccumulator[i].Update(fdv_val[i]);
      	for(int j=0;j<niv;j++){
      	  fDVIVCovAcccumulator[i*niv+j].Update(fdv_val[i],fiv_val[j]);
      	}
      }
    } // end of event loop
    
    for(int i=0;i<niv;i++)
      for(int j=0;j<niv;j++)
    	fIVCovAcccumulator[i*niv+j].UpdateStat(fIVCovStat[i*niv+j]);
  
    for(int i=0;i<ndv;i++)
      for(int j=0;j<niv;j++)
    	fDVIVCovAcccumulator[i*niv+j].UpdateStat(fDVIVCovStat[i*niv+j]);
  
    for(int i=0;i<ndv+niv;i++)
      fVarAcccumulator[i].UpdateStat(fVarStat[i]);

    cov_tree->Fill();
    input_file->Close();
    iter++;
  }// end of loop over run;
  
  output_file->cd();
  cov_tree->Write();
  output_file->Close();
}

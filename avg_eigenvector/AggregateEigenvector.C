void AggregateEigenvector(){
  TFile *input = TFile::Open("rootfiles/all_slugs.root");
  TTree *eig_tree = (TTree*)input->Get("eig");
  eig_tree->AddFriend("reg");
  TString cut_text = "arm_flag==0";
  int nIV = 12;
  vector<double> fMean;
  vector<double> fRMS;
  TF1 *pol0 =new TF1("pol0","pol0",-1e6,1e6);
  for(int i=0; i<nIV;i++){
    eig_tree->Draw(Form("spin*diff_evMon%d.mean*1e6:reg_asym_us_avg.err*1e9:Entry$",i),
		   cut_text,"goff");

    TGraphErrors* ger_mean = new TGraphErrors(eig_tree->GetSelectedRows(),
					      eig_tree->GetV3(),eig_tree->GetV1(),
					      0, eig_tree->GetV2());
    ger_mean->Fit(pol0,"Q0");
    fMean.push_back(pol0->GetParameter(0));
    eig_tree->Draw(Form("diff_evMon%d.rms*1e3:reg_asym_us_avg.err*1e9:Entry$",i),
		   cut_text,"goff");

    TGraphErrors* ger_rms = new TGraphErrors(eig_tree->GetSelectedRows(),
					     eig_tree->GetV3(),eig_tree->GetV1(),
					      0, eig_tree->GetV2());
    ger_rms->Fit(pol0,"Q0");
    fRMS.push_back(pol0->GetParameter(0));
  }

  for(int iev=0;iev<nIV;iev++){
    cout << "evMon"<<iev <<  "\t &";
    printf("%.2f \t &", fMean[iev]);
    printf("%.2f \t \\\\ \n", fRMS[iev]);
  }
}

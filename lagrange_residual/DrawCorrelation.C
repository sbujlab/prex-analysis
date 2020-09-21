void DrawCorrelation(){
  gStyle->SetOptStat(111111);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.35);
  TChain *res = new TChain("res");
  for(int i=1;i<=94;i++){
    if(i==7)
      continue;
    res->Add(Form("rootfiles/slug%d_lagrange_residual.root",i));
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Print("lagrange_all_run_residual_correlation.pdf[");
  res->SetMarkerStyle(20);
  for(int icoil=1;icoil<=7;icoil++){
    res->Draw(Form("usl_coil%d_res*1e6 : usr_coil%d_res*1e6 " ,icoil,icoil),
	      Form("arm_flag==0 && usl_coil%d_flag && usr_coil%d_flag",icoil,icoil),"");
    c1->Print("lagrange_all_run_residual_correlation.pdf");
  }


  c1->Print("lagrange_all_run_residual_correlation.pdf]");
}

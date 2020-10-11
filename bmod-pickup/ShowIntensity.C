void ShowIntensity(Int_t run_number=3426,
			     TString user_cut=""){
  TString qw_path = "/lustre19/expphy/volatile/halla/parity/prex-respin2/japanOutput/";
  TString postpan_path = "/lustre19/expphy/volatile/halla/parity/prex-respin2/postpan_respin/";
  TString lagrange_path = "/lustre/expphy/volatile/halla/parity/prex-respin2/LagrangeOutput/rootfiles/";
  TFile* japan = TFile::Open(qw_path+Form("prexPrompt_pass2_%d.000.root",run_number));
  TFile* postpan = TFile::Open(postpan_path+Form("prexPrompt_%d_000_regress_overload.root",
						 run_number));
  TFile* lagrange = TFile::Open(lagrange_path+Form("prexRespin2_lagrange_eigen_%d.000.root",
						   run_number));
  
  TTree *mul = (TTree*)japan->Get("mul");
  TTree *mulc = (TTree*)japan->Get("mulc");
  TTree *reg = (TTree*)postpan->Get("reg");
  TTree *lagr = (TTree*)lagrange->Get("lagrall");
  mul->AddFriend(mulc);
  mul->AddFriend(reg);
  mul->AddFriend(lagr);

  mul->SetMarkerStyle(7);

  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  c3->cd();
  mul->Draw("reg_asym_sam_37_avg:Entry$","ErrorFlag==0x4019000 && BurstCounter==0","");
  
}

void ShowResidualCorrelation(Int_t run_number=3426,
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
  TTree *reg = (TTree*)lagrange->Get("regall");
  TTree *lagr = (TTree*)lagrange->Get("lagrall");
  mul->AddFriend(mulc);
  mul->AddFriend(reg);
  mul->AddFriend(lagr);

  mul->SetMarkerStyle(7);

  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  c3->cd();
  mul->Draw("reg_asym_us_dd/ppm:Entry$","ErrorFlag==0x4019000 && BurstCounter==0","*");
  TH2D *h2d_lagr_cycles  = (TH2D*)gPad->FindObject("htemp");
  h2d_lagr_cycles->SetTitle("Corrected Asymmetry(ppm) in BMW data");
  TGraph *g_lagr_cycles = (TGraph*)gPad->FindObject("Graph");
  g_lagr_cycles->SetName("all_cycles");
  g_lagr_cycles->SetMarkerColor(kWhite);
  g_lagr_cycles->SetMarkerStyle(7);

  mul->Draw("reg_asym_us_dd/ppm:Entry$","ErrorFlag==0x4019000 && BurstCounter==0 && actual_pattern_polarity==1","*same");
  TGraph *g_lagr_cycles_pos = (TGraph*)gPad->FindObject("Graph");
  g_lagr_cycles_pos->SetName("all_cycles_pos");
  g_lagr_cycles_pos->SetMarkerColor(kBlue);
  g_lagr_cycles_pos->SetMarkerStyle(7);

  mul->Draw("reg_asym_us_dd/ppm:Entry$","ErrorFlag==0x4019000 && BurstCounter==0 && actual_pattern_polarity==0","*same");
  TGraph *g_lagr_cycles_neg = (TGraph*)gPad->FindObject("Graph");
  g_lagr_cycles_neg->SetName("all_cycles_neg");
  g_lagr_cycles_neg->SetMarkerColor(kRed);
  g_lagr_cycles_neg->SetMarkerStyle(7);
  TLegend *leg_lagr_cycles = new TLegend(0.9,0.9,0.7,0.7);
  leg_lagr_cycles->AddEntry(g_lagr_cycles_pos,"Pattern Polarity Pos.","p");
  leg_lagr_cycles->AddEntry(g_lagr_cycles_neg,"Pattern Polarity Neg.","p");
  leg_lagr_cycles->Draw("same");


  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  mul->Draw("lagr_asym_us_avg/ppm:diff_bpm4eX/um","ErrorFlag==0x4019000 && yield_bmwobj==5 && BurstCounter==0 ","");
  TH2D *h2dCorCoil5 = (TH2D*)gPad->FindObject("htemp");
  h2dCorCoil5->SetTitle("Coil5: Corrected AVG Asymmetry(ppm)");
  h2dCorCoil5->GetYaxis()->Set(100000,-2e3,2.2e3);
  c1->cd(2);
  mul->Draw("(actual_pattern_polarity-0.5)*2*lagr_asym_us_avg/ppm:(actual_pattern_polarity-0.5)*2*diff_bpm4eX/um","ErrorFlag==0x4019000 && yield_bmwobj==5 && BurstCounter==0 ","");
  TH2D *h2dCorUnsigned = (TH2D*)gPad->FindObject("htemp");
  h2dCorUnsigned->SetTitle("Coil5: Corrected AVG Asymmetry(ppm), polarity unsigned");
  h2dCorUnsigned->GetYaxis()->Set(100000,-2e3,2.2e3);

  TCanvas *c4 = new TCanvas("c4","c4",1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  mul->Draw("lagr_asym_us_dd/ppm:diff_bpm4eX/um",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && BurstCounter==0 ","");
  TH2D *h2dCorCoil5_dd = (TH2D*)gPad->FindObject("htemp");
  h2dCorCoil5_dd->SetTitle("Coil5: Corrected DD Asymmetry(ppm)");
  c4->cd(2);
  mul->Draw("(actual_pattern_polarity-0.5)*2*lagr_asym_us_dd/ppm:(actual_pattern_polarity-0.5)*2*diff_bpm4eX/um",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && BurstCounter==0 ","");
  TH2D *h2dCorUnsigned_dd = (TH2D*)gPad->FindObject("htemp");
  h2dCorUnsigned_dd->SetTitle("Coil5: Corrected DD Asymmetry(ppm), polarity unsigned");

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  mul->Draw("asym_us_avg/ppm:diff_bpm4eX/um",
	    "(ErrorFlag&0xfbfe6fff)==0 && BurstCounter==0 && (yield_bmwobj==5 || yield_bmwobj<=0) ","");
  TH2D *h2dRaw = (TH2D*)gPad->FindObject("htemp");
  h2dRaw->SetTitle("Raw AVG Asymmetry(ppm)");

  c2->cd(2);
  mul->Draw("lagr_asym_us_avg/ppm:diff_bpm4eX/um",
	    "(ErrorFlag&0xfbfe6fff)==0 && BurstCounter==0&&(yield_bmwobj==5||yield_bmwobj<=0)","");
  TGraph *gcor = (TGraph*)gPad->FindObject("Graph");
  TH2D *h2dCor = (TH2D*)gPad->FindObject("htemp");
  h2dCor->SetTitle("Corrected AVG Asymmetry(ppm)");
  h2dCor->GetYaxis()->Set(100000,-2e3,2.2e3);

}

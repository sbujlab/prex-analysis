void ShowPhaseDependent(Int_t run_number=3426,
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

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  mul->Draw("lagr_asym_usl/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5","*");
  TH2D *huslUnsigned = (TH2D*)gPad->FindObject("htemp");
  huslUnsigned->SetTitle("Coil5: Corrected USL asymmetry(ppm)");
  TGraph *gusl = (TGraph*)gPad->FindObject("Graph");
  gusl->SetName("all");
  gusl->SetMarkerColor(kWhite);
  mul->Draw("lagr_asym_usl/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && actual_pattern_polarity==1","*same");
  TGraph *gusl_pos = (TGraph*)gPad->FindObject("Graph");
  gusl_pos->SetName("all");
  gusl_pos->SetMarkerColor(kBlue);
  gusl_pos->SetMarkerStyle(7);

  mul->Draw("lagr_asym_usl/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && actual_pattern_polarity==0","*same");
  TGraph *gusl_neg = (TGraph*)gPad->FindObject("Graph");
  gusl_neg->SetName("all");
  gusl_neg->SetMarkerColor(kRed);
  gusl_neg->SetMarkerStyle(7);
  TLegend *leg_usl = new TLegend(0.1,0.9,0.35,0.7);
  leg_usl->AddEntry(gusl_pos,"Pattern Polarity Pos.","p");
  leg_usl->AddEntry(gusl_neg,"Pattern Polarity Neg.","p");
  leg_usl->Draw("same");

  c1->cd(2);
  mul->Draw("(actual_pattern_polarity-0.5)*2*lagr_asym_usl/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5","");
  TH2D* h2d_usl_signed = (TH2D*)gPad->FindObject("htemp");
  h2d_usl_signed->SetTitle("Coil5: Corrected usl asymmetry(ppm), polarity unsigned");

  // USL
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  mul->Draw("lagr_asym_usr/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5","*");
  TH2D *husrUnsigned = (TH2D*)gPad->FindObject("htemp");
  husrUnsigned->SetTitle("Coil5: Corrected USR asymmetry(ppm)");
  TGraph *gusr = (TGraph*)gPad->FindObject("Graph");
  gusr->SetName("all");
  gusr->SetMarkerColor(kWhite);
  mul->Draw("lagr_asym_usr/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && actual_pattern_polarity==1","*same");
  TGraph *gusr_pos = (TGraph*)gPad->FindObject("Graph");
  gusr_pos->SetName("all");
  gusr_pos->SetMarkerColor(kBlue);
  gusr_pos->SetMarkerStyle(7);

  mul->Draw("lagr_asym_usr/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && actual_pattern_polarity==0","*same");
  TGraph *gusr_neg = (TGraph*)gPad->FindObject("Graph");
  gusr_neg->SetName("all");
  gusr_neg->SetMarkerColor(kRed);
  gusr_neg->SetMarkerStyle(7);
  TLegend *leg_usr = new TLegend(0.1,0.9,0.35,0.7);
  leg_usr->AddEntry(gusr_pos,"Pattern Polarity Pos.","p");
  leg_usr->AddEntry(gusr_neg,"Pattern Polarity Neg.","p");
  leg_usr->Draw("same");

  c2->cd(2);
  mul->Draw("(actual_pattern_polarity-0.5)*2*lagr_asym_usr/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5","");
  TH2D* h2d_usr_signed = (TH2D*)gPad->FindObject("htemp");
  h2d_usr_signed->SetTitle("Coil5: Corrected usr asymmetry(ppm), polarity unsigned");

  // US DD
  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  mul->Draw("lagr_asym_us_dd/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5","*");
  TH2D *hus_ddUnsigned = (TH2D*)gPad->FindObject("htemp");
  hus_ddUnsigned->SetTitle("Coil5: Corrected US_DD asymmetry(ppm)");
  TGraph *gus_dd = (TGraph*)gPad->FindObject("Graph");
  gus_dd->SetName("all");
  gus_dd->SetMarkerColor(kWhite);
  mul->Draw("lagr_asym_us_dd/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && actual_pattern_polarity==1","*same");
  TGraph *gus_dd_pos = (TGraph*)gPad->FindObject("Graph");
  gus_dd_pos->SetName("all");
  gus_dd_pos->SetMarkerColor(kBlue);
  gus_dd_pos->SetMarkerStyle(7);

  mul->Draw("lagr_asym_us_dd/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5 && actual_pattern_polarity==0","*same");
  TGraph *gus_dd_neg = (TGraph*)gPad->FindObject("Graph");
  gus_dd_neg->SetName("all");
  gus_dd_neg->SetMarkerColor(kRed);
  gus_dd_neg->SetMarkerStyle(7);
  TLegend *leg_us_dd = new TLegend(0.1,0.9,0.35,0.7);
  leg_us_dd->AddEntry(gus_dd_pos,"Pattern Polarity Pos.","p");
  leg_us_dd->AddEntry(gus_dd_neg,"Pattern Polarity Neg.","p");
  leg_us_dd->Draw("same");

  c3->cd(2);
  mul->Draw("(actual_pattern_polarity-0.5)*2*lagr_asym_us_dd/ppm:yield_bmod_ramp",
	    "ErrorFlag==0x4019000 && yield_bmwobj==5","");
  TH2D* h2d_us_dd_signed = (TH2D*)gPad->FindObject("htemp");
  h2d_us_dd_signed->SetTitle("Coil5: Corrected usdd asymmetry(ppm), polarity unsigned");

  TCanvas *c4 = new TCanvas
    
}

void ShowBmodPattern(Int_t run_number=3664,
		     Int_t bmwcycnum=872){
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

  TTree *evt = (TTree*)japan->Get("evt");
  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->Divide(3,1);
  c1->cd(1);
  evt->Draw("usl:bmod_ramp",
	    Form("ErrorFlag==0x4019000 && bmwobj==5 && bmwcycnum==%d",bmwcycnum),"*");
  TH2D *husl = (TH2D*)gPad->FindObject("htemp");
  husl->SetTitle("usl yield vs bmod ramp, group by pattern ");
  TGraph *gallusl = (TGraph*)gPad->FindObject("Graph");
  gallusl->SetName("gall");
  gallusl->SetMarkerColor(kWhite);
  evt->Draw("usl:bmod_ramp",
	    Form("ErrorFlag==0x4019000 && bmwobj==5 && bmwcycnum==%d && pattern_number%%2==0",bmwcycnum),"*same");
  TGraph *gevenusl = (TGraph*)gPad->FindObject("Graph");
  gevenusl->SetName("geven");
  gevenusl->SetMarkerColor(kBlue);
  gevenusl->SetMarkerStyle(20);
  evt->Draw("usl:bmod_ramp",
	    Form("ErrorFlag==0x4019000 && bmwobj==5 && bmwcycnum==%d && pattern_number%%2==1",bmwcycnum),"*same");
  TGraph *goddusl = (TGraph*)gPad->FindObject("Graph");
  goddusl->SetName("godd");
  goddusl->SetMarkerColor(kRed);
  goddusl->SetMarkerStyle(20);

  // USR
  c1->cd(2);
  evt->Draw("usr:bmod_ramp",
	    Form("ErrorFlag==0x4019000 && bmwobj==5 && bmwcycnum==%d",bmwcycnum),"*");
  TH2D *husr = (TH2D*)gPad->FindObject("htemp");
  husr->SetTitle("usr yield vs bmod ramp, group by pattern ");
  TGraph *gallusr = (TGraph*)gPad->FindObject("Graph");
  gallusr->SetName("gall");
  gallusr->SetMarkerColor(kWhite);

  evt->Draw("usr:bmod_ramp",
	    Form("ErrorFlag==0x4019000 && bmwobj==5 && bmwcycnum==%d && pattern_number%%2==0",bmwcycnum),"*same");
  TGraph *gevenusr = (TGraph*)gPad->FindObject("Graph");
  gevenusr->SetName("geven");
  gevenusr->SetMarkerColor(kBlue);
  gevenusr->SetMarkerStyle(20);
  evt->Draw("usr:bmod_ramp",
	    Form("ErrorFlag==0x4019000 && bmwobj==5 && bmwcycnum==%d && pattern_number%%2==1",bmwcycnum),"*same");
  TGraph *goddusr = (TGraph*)gPad->FindObject("Graph");
  goddusr->SetName("godd");
  goddusr->SetMarkerColor(kRed);
  goddusr->SetMarkerStyle(20);

  c1->cd(3);
  mul->Draw("lagr_asym_us_avg/ppm",
	    Form("ErrorFlag==0x4019000 && yield_bmwobj==5 && yield_bmwcycnum==%d",bmwcycnum));
  
}  

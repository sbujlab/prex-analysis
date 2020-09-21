void DrawEigVecFromDance(Int_t slug=94){
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  
  TFile *input = TFile::Open(Form("./rootfiles/MergedLagrange_slug%d.root",slug));
  TTree *eig_tree = (TTree*)input->Get("mini_eig");
  eig_tree->AddFriend("mini");
  eig_tree->AddFriend("slope",Form("./rootfiles/dit_eigslopes_slug%d.root",slug));
  TString  leaflist =  eig_tree->GetBranch("eigvec")->GetTitle() ;

  vector< TString > IVlist1 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm8X","bpm12X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm8Y","bpm12Y"};
  
  vector< TString > IVlist2 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm11X","bpm12X","bpm16X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm11Y","bpm12Y","bpm16Y"};
  vector<TString> fBPMList;
  if(slug>=3)
    fBPMList = IVlist2;
  else
    fBPMList = IVlist1;

  Int_t nBPM = fBPMList.size();
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Print(Form("plots/slug%d_vector.pdf[",slug));
  // Eigen Values
  TMultiGraph *mgval = new TMultiGraph();
  TLegend *legval = new TLegend(0.9,0.6,1,0.9);
  Int_t color_offset=0;
  for(int j=0;j<nBPM;j++){
    Int_t npt =  eig_tree->Draw(Form("diff_evMon%d.rms/um:Entry$",j),
				"","goff");
    Double_t *y_val = eig_tree->GetV1();
    Double_t *fEntries = eig_tree->GetV2();
    TGraph *g1 = new TGraph(npt,fEntries,y_val);
    if(j+1==5 || j+1==9)
      color_offset ++;
    g1->SetLineColor(j+1+color_offset);
    g1->SetMarkerColor(j+1+color_offset);
    g1->SetMarkerStyle(20);
    legval->AddEntry(g1,Form("Eigenvector %d",j),"lp");
    mgval->Add(g1,"lp");
  } // component loop
  c2->cd();
  c2->SetFillStyle(4000);
  c2->SetFrameFillStyle(4000);
  gPad->SetFillStyle(4000);
  mgval->Draw("ALP");
  mgval->SetTitle("#sqrt{Eigenvalues} (um);minirun ; um");\
  legval->SetFillStyle(4000);
  legval->Draw("same");
  c2->Print("eigen_values.pdf");
  // c2->Print(Form("plots/slug%d_vector.pdf",slug));
  c2->Clear("D");
  //Eigen vectors
  for(int i=0;i<nBPM;i++){
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.9,0.6,0.95,0.9);
    Int_t color_offset=0;
    for(int j=0;j<nBPM;j++){
      Int_t npt =  eig_tree->Draw(Form("sign%d*evMon%d_%s:Entry$",i,i,fBPMList[j].Data()),
				  "","goff");
      Double_t *y_val = eig_tree->GetV1();
      Double_t *fEntries = eig_tree->GetV2();
      TGraph *g1 = new TGraph(npt,fEntries,y_val);
      if(j+1==5 || j+1==9)
	color_offset ++;
      g1->SetLineColor(j+1+color_offset);
      g1->SetMarkerColor(j+1+color_offset);
      g1->SetMarkerStyle(20);
      leg->AddEntry(g1,fBPMList[j],"lp");
      mg->Add(g1,"lp");
    } // component loop
    c2->cd();
    mg->Draw("ALP");
    mg->SetTitle(Form("EigenVector #%d;minirun ;Magnitude",i));
    leg->Draw("same");
    c2->Print(Form("plots/slug%d_vector.pdf",slug));
  } // eigenvalue loop

  c2->Print(Form("plots/slug%d_vector.pdf]",slug));
  
  input->Close();
  
}

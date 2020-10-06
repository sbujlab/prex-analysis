void DrawEigenValuesBySlug( Int_t slug=94){
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  vector<TString> IVlist;  
  TString filename_tag;
  TString slug_title = Form("Slug%d",slug);
  vector< TString > IVlist1 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm8X","bpm12X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm8Y","bpm12Y"};
  vector< TString > IVlist2 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm11X","bpm12X","bpm16X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm11Y","bpm12Y","bpm16Y"};
  if(slug>=3)
    IVlist = IVlist2;
  else
    IVlist = IVlist1;
  filename_tag = "allbpm";

  TFile *input;
  if(slug>=3)
    input= TFile::Open(Form("./rootfiles/slug%d_sorted_eigenvector_allbpm.root",slug));
  else
    input= TFile::Open(Form("./rootfiles/slug%d_sorted_eigenvector_10bpm.root",slug));
  TTree *eig_tree;
  eig_tree = (TTree*)input->Get("eig");
  Int_t nBPM = IVlist.size();
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  // Eigen Values
  c2->cd();
  Int_t color_offset=0;
  vector<Int_t> fSlug;
  vector<Int_t> fRun;
  eig_tree->Draw("slug:run","","goff");
  Double_t *yval = eig_tree->GetV1();
  Double_t *xval = eig_tree->GetV2();
  Int_t nmini = eig_tree->GetSelectedRows();
  for(int i =0;i<nmini;i++){
    fSlug.push_back(yval[i]);
    fRun.push_back(xval[i]);
  }

  TMultiGraph *mgval = new TMultiGraph();
  TLegend *legval = new TLegend(0.9,0.6,0.99,0.9);

  for(int j=0;j<nBPM;j++){
    Int_t npt =  eig_tree->Draw(Form("diff_evMon%d.rms*1e3:Entry$",j),
  				"","goff");
    Double_t *y_val = eig_tree->GetV1();
    Double_t *fEntries = eig_tree->GetV2();
    TGraph *g1 = new TGraph(npt,fEntries,y_val);
    int bpm_index = j % (nBPM/2)+1;
    if(bpm_index>=5)
      color_offset=1;
    else
      color_offset =0;
    g1->SetLineColor(bpm_index+color_offset);
    g1->SetMarkerColor(bpm_index+color_offset);
    g1->SetMarkerStyle(20);
    legval->AddEntry(g1,Form("Eigenvector #%d",j),"lp");
    mgval->Add(g1,"lp");
  } // component loop
  mgval->Draw("ALP");
  TH1D *hmgval = (TH1D*)mgval->GetHistogram();
  hmgval->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
  Int_t last_run_label=0;
  double ymax = hmgval->GetYaxis()->GetXmax();
  double ymin = hmgval->GetYaxis()->GetXmin();
  for(int i=0;i<nmini;i++){
    if(fRun[i]!=last_run_label){
      last_run_label = fRun[i];
      hmgval->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
      TLine *line_buff = new TLine(i,ymin,i,ymax);
      line_buff->SetLineWidth(1);
      line_buff->SetLineStyle(4);
      line_buff->SetLineColor(14);
      line_buff->Draw("same");
    }else
      hmgval->GetXaxis()->SetBinLabel(i+1,"");
  }
  hmgval->SetTitle(slug_title+": #sqrt{Eigenvalues} (um); Run Number; um");
  hmgval->GetXaxis()->SetTitleOffset(1.25);
  legval->Draw("same");
  c2->Print(Form("plots/slug%d_eigenvalue.pdf",slug));

  input->Close();
}

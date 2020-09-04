void DrawEigenVectorsBySlug( Int_t slug=94){
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
  c2->Print(Form("plots/slug%d_eigenvector.pdf[",slug));
  // Eigen Values
  c2->cd();
  int nPad =4;
  vector<TPad*> fPad(nPad);
  Double_t header_height =0.05;
  Double_t footer_height =0.04;
  Double_t pad_height = (1 - header_height -footer_height)/(nPad);
  for(int i=0;i<nPad;i++){
    double y1 = 1.0-pad_height *(i+1) - header_height;
    double y2 = y1+pad_height;
    if(i==0)
      y2+=header_height;
    if(i+1>=nPad)
      y1 = 0.0;
    fPad[i] =  new TPad("","",0.0,y1,1.0,y2);
    if(i!=0)
      fPad[i]->SetTopMargin(0.01);
    if(i+1<nPad)
      fPad[i]->SetBottomMargin(0.01);
    else
      fPad[i]->SetBottomMargin(0.2);
    fPad[i]->Draw();
  }

  TMultiGraph *mgval = new TMultiGraph();
  TLegend *legval = new TLegend(0.9,0.6,0.99,0.9);
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

  vector<TMultiGraph *> fmgval(nPad);
  vector<TLegend *> fLegval(nPad);
  for(int ip=0;ip<nPad;ip++){
    fmgval[ip]  = new TMultiGraph();
    fLegval[ip] =  new TLegend(0.9,0.6,0.99,0.9);
  }

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
    g1->SetMarkerStyle(6);
    
    fLegval[j/3]->AddEntry(g1,Form("Eigenvector #%d",j),"lp");
    fmgval[j/3]->Add(g1,"lp");
  } // component loop
  for(int ip=0;ip<nPad;ip++){
    fPad[ip]->cd();
    fmgval[ip]->Draw("A");
    TH1D *hmgval = (TH1D*)fmgval[ip]->GetHistogram();
    hmgval->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
    Int_t last_slug_label=0;
    double ymax = hmgval->GetYaxis()->GetXmax();
    double ymin = hmgval->GetYaxis()->GetXmin();
    for(int i=0;i<nmini;i++){
      if(fSlug[i]!=last_slug_label){
	last_slug_label = fSlug[i];
	hmgval->GetXaxis()->SetBinLabel(i+1, Form("%d",fSlug[i]));
      }else
	hmgval->GetXaxis()->SetBinLabel(i+1,"");
    }
    if(ip==0)
      hmgval->SetTitle(" #sqrt{Eigenvalues} (um); Slug Number; um");
    if(ip+1>=nBPM+3)
      hmgval->SetTitle(";Slug Number ;Magnitude");
    	
    hmgval->GetXaxis()->SetTitleOffset(1.0);
    hmgval->GetYaxis()->SetTitleOffset(0.5);
		     
    hmgval->GetXaxis()->SetTitleSize(0.07);
    hmgval->GetYaxis()->SetTitleSize(0.07);
		     
    hmgval->GetXaxis()->SetLabelSize(0.09);
    hmgval->GetYaxis()->SetLabelSize(0.07);
    fLegval[ip]->Draw("same");
  }

  // c2->cd();
  // mgval->Draw("ALP");
  // TH1D *hmgval = (TH1D*)mgval->GetHistogram();
  // hmgval->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
  Int_t last_run_label=0;
  double ymax;// = hmgval->GetYaxis()->GetXmax();
  double ymin;// = hmgval->GetYaxis()->GetXmin();
  // for(int i=0;i<nmini;i++){
  //   if(fRun[i]!=last_run_label){
  //     last_run_label = fRun[i];
  //     hmgval->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
  //     TLine *line_buff = new TLine(i,ymin,i,ymax);
  //     line_buff->SetLineWidth(1);
  //     line_buff->SetLineStyle(4);
  //     line_buff->SetLineColor(14);
  //     line_buff->Draw("same");
  //   }else
  //     hmgval->GetXaxis()->SetBinLabel(i+1,"");
  // }
  // hmgval->SetTitle(slug_title+": #sqrt{Eigenvalues} (um); Run Number; um");
  // hmgval->GetXaxis()->SetTitleOffset(1.25);
  // legval->Draw("same");
  c2->Print(Form("plots/slug%d_eigenvector.pdf",slug));
  c2->Clear();
  for(int i=0;i<nBPM;i++){
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.9,0.6,0.99,0.9);
    Int_t color_offset=0;
    for(int j=0;j<nBPM;j++){
      Int_t npt =  eig_tree->Draw(Form("evMon%d_%s:Entry$",i,IVlist[j].Data()),
				  "","goff");
      Double_t *y_val = eig_tree->GetV1();
      Double_t *fEntries = eig_tree->GetV2();
      TGraph *g1 = new TGraph(npt,fEntries,y_val);
      if(j+1==5 || j+1==9)
	color_offset ++;
      g1->SetLineColor(j+1+color_offset);
      g1->SetMarkerColor(j+1+color_offset);
      g1->SetMarkerStyle(20);
      leg->AddEntry(g1,IVlist[j],"lp");
      mg->Add(g1,"lp");
    } // component loop
    c2->cd();
    mg->Draw("ALP");
    TH1D *hmg = (TH1D*)mg->GetHistogram();
    hmg->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
    last_run_label = 0;
    ymax  = hmg->GetYaxis()->GetXmax();
    ymin  = hmg->GetYaxis()->GetXmin();
    for(int i=0;i<nmini;i++){
      if(fRun[i]!=last_run_label){
	last_run_label = fRun[i];
	hmg->GetXaxis()->SetBinLabel(i+1, Form("%d",fRun[i]));
	TLine *line_buff= new TLine(i,ymin,i,ymax);
	line_buff->SetLineWidth(1);
	line_buff->SetLineStyle(4);
	line_buff->SetLineColor(14);
	line_buff->Draw("same");
      }else
	hmg->GetXaxis()->SetBinLabel(i+1,"");
    }
    hmg->SetTitle(slug_title+Form(": EigenVector #%d;Run Number ;Magnitude",i));
    hmg->GetXaxis()->SetTitleOffset(1.25);
    leg->Draw("same");
    c2->Print(Form("plots/slug%d_eigenvector.pdf",slug));

  } // eigenvalue loop
  c2->Print(Form("plots/slug%d_eigenvector.pdf]",slug));

  input->Close();
}

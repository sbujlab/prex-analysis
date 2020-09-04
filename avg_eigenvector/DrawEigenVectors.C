void DrawEigenVectors(){
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  vector<TString> IVlist;  
  vector< TString > IVlist1 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm8X","bpm12X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm8Y","bpm12Y"};
  vector< TString > IVlist2 =  {"bpm4aX","bpm4eX","bpm1X",
				"bpm11X","bpm12X","bpm16X",
				"bpm4aY","bpm4eY","bpm1Y",
				"bpm11Y","bpm12Y","bpm16Y"};
  IVlist = IVlist2;

  TFile *input = TFile::Open("./rootfiles/sorted_eigenvector_allbpm.root");
  TTree *eig_tree;
  eig_tree = (TTree*)input->Get("eig");
  Int_t nBPM = IVlist.size();
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,1000);
  c2->cd();
  vector<TPad*> fPad(nBPM/3);
  Double_t header_height =0.05;
  Double_t footer_height =0.04;
  Double_t pad_height = (1 - header_height -footer_height)/(nBPM/3);
  for(int i=0;i<nBPM/3;i++){
    double y1 = 1.0-pad_height *(i+1) - header_height;
    double y2 = y1+pad_height;
    if(i==0)
      y2+=header_height;
    if(i+1>=nBPM/3)
      y1 = 0.0;
    fPad[i] =  new TPad("","",0.0,y1,1.0,y2);
    if(i!=0)
      fPad[i]->SetTopMargin(0.01);
    if(i+1<nBPM/3)
      fPad[i]->SetBottomMargin(0.01);
    else
      fPad[i]->SetBottomMargin(0.2);
    fPad[i]->Draw();
  }

  c2->Print("eigenvector.pdf[");
  
  // Eigen Values
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

  vector<TMultiGraph *> fmgval(nBPM/3);
  vector<TLegend *> fLegval(nBPM/3);
  for(int ip=0;ip<nBPM/3;ip++){
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

  for(int ip=0;ip<nBPM/3;ip++){
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
  c2->Print("eigenvector.pdf");
  
  //Eigen vectors
  for(int i=0;i<nBPM;i++){ 
    vector<TMultiGraph *> fmg(nBPM/3);
    vector<TLegend *> fLeg(nBPM/3);
    for(int ip=0;ip<nBPM/3;ip++){
      fmg[ip]  = new TMultiGraph();
      fLeg[ip] =  new TLegend(0.9,0.6,0.99,0.9);
    }
    int color_offset=0;
    for(int j=0;j<nBPM;j++){
      Int_t npt =  eig_tree->Draw(Form("evMon%d_%s:Entry$",i,IVlist[j].Data()),
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

      fLeg[j/3]->AddEntry(g1,IVlist[j],"lp");
      fmg[j/3]->Add(g1,"lp");
    } // component loop
    for(int ip=0;ip<nBPM/3;ip++){
      fPad[ip]->cd();
      fmg[ip]->Draw("A");
      TH1D *hmg = (TH1D*)fmg[ip]->GetHistogram();
      hmg->GetXaxis()->Set(nmini,-0.5,nmini-0.5);
      int last_slug_label = 0;
      double ymax  = hmg->GetYaxis()->GetXmax();
      double ymin  = hmg->GetYaxis()->GetXmin();
      for(int i=0;i<nmini;i++){
	if(fSlug[i]!=last_slug_label){
	  last_slug_label = fSlug[i];
	  hmg->GetXaxis()->SetBinLabel(i+1, Form("%d",fSlug[i]));
	  // TLine *line_buff= new TLine(i,ymin,i,ymax);
	  // line_buff->SetLineWidth(1);
	  // line_buff->SetLineStyle(4);
	  // line_buff->SetLineColor(14);
	  // line_buff->Draw("same");
	}else
	  hmg->GetXaxis()->SetBinLabel(i+1,"");
      }
      if(ip==0){
	hmg->SetTitle(Form("EigenVector #%d;Slug Number ;Magnitude",i));
      }
      if(ip+1>=nBPM/3)
	hmg->SetTitle(";Slug Number ;Magnitude");
	
      hmg->GetXaxis()->SetTitleOffset(1.0);
      hmg->GetYaxis()->SetTitleOffset(0.5);
		     
      hmg->GetXaxis()->SetTitleSize(0.07);
      hmg->GetYaxis()->SetTitleSize(0.07);
		     
      hmg->GetXaxis()->SetLabelSize(0.09);
      hmg->GetYaxis()->SetLabelSize(0.07);

      fLeg[ip]->Draw("same");
    }

    c2->Print("eigenvector.pdf");
  } // eigenvalue loop
  c2->Print("eigenvector.pdf]");
  input->Close();
}

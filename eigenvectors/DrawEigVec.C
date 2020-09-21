void DrawEigVec(Int_t slug=94){
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  
  TFile *input = TFile::Open(Form("./rootfiles/slug%d.root",slug));
  TTree *eig_tree = (TTree*)input->Get("eig");
  TString  leaflist =  eig_tree->GetBranch("eigvec")->GetTitle() ;
  vector<TString> fBPMList;
  while(leaflist.Contains(":")){
    Ssiz_t pos = leaflist.First(':');
    TString extracted = TString( leaflist(0,pos) );
    extracted.ReplaceAll("/D","");
    fBPMList.push_back(extracted);
    leaflist.Remove(0,pos+1);
    cout << extracted << endl;
  }
  leaflist.ReplaceAll("/D","");
  fBPMList.push_back(leaflist);
  cout << leaflist << endl;
  Int_t nBPM = fBPMList.size();
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Print(Form("plots/slug%d_2d.pdf[",slug));
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Print(Form("plots/slug%d_chart.pdf[",slug));
  for(int i=0;i<nBPM;i++){
    c1->cd();
    Int_t nEntries  = eig_tree->GetEntries("MyStat>8000 && eigID==0");
    TH2D h2d("h2d","",nEntries,0,nEntries,nBPM,0,nBPM);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.9,0.6,0.95,0.9);
    Int_t color_offset=0;
    for(int j=0;j<nBPM;j++){
      Int_t npt =  eig_tree->Draw("eigvec."+fBPMList[j]+":run:minirun:Entry$",
				  Form("MyStat>8000 && eigID==%d",i),"goff");
      Double_t *y_val = eig_tree->GetV1();
      Double_t *frun = eig_tree->GetV2();
      Double_t *fMini = eig_tree->GetV3();
      Double_t *fEntries = eig_tree->GetV4();
      for(int ipt=0;ipt<npt;ipt++){
	h2d.Fill( Form("%.0f:%.0f",frun[ipt],fMini[ipt]),
		  j, y_val[ipt]);
      }
      TGraph *g1 = new TGraph(npt,fEntries,y_val);
      if(j+1==5 || j+1==9)
	color_offset ++;
      g1->SetLineColor(j+1+color_offset);
      g1->SetMarkerColor(j+1+color_offset);
      g1->SetMarkerStyle(20);
      leg->AddEntry(g1,fBPMList[j],"lp");
      mg->Add(g1,"lp");
    } // component loop
    h2d.Draw("COLZ");
    h2d.GetZaxis()->SetRangeUser(-1,1);
    h2d.SetTitle(Form(" EigenVector # %d vs minirun ",i));
    h2d.GetYaxis()->SetNdivisions(nBPM);
    h2d.GetYaxis()->SetLabelSize(0.05);
    for(int j=0;j<nBPM;j++){
      h2d.GetYaxis()->SetBinLabel(j+1, fBPMList[j]);
    }
    
    h2d.GetXaxis()->SetTitle("Run: Minirun");
    c1->Print(Form("plots/slug%d_2d.pdf",slug));
    c2->cd();
    mg->Draw("ALP");
    mg->SetTitle(Form("EigenVector #%d;minirun ;Magnitude",i));
    leg->Draw("same");
    c2->Print(Form("plots/slug%d_chart.pdf",slug));
  } // eigenvalue loop
  c1->Print(Form("plots/slug%d_2d.pdf]",slug));
  c2->Print(Form("plots/slug%d_chart.pdf]",slug));
  c1->cd();
  TMultiGraph *mg =new TMultiGraph();
  TLegend *leg = new TLegend(0.9,0.7,1.0,0.9);
  Int_t color_offset=0;
  for(int j=0;j<nBPM;j++){
    eig_tree->Draw("sqrt(eigval)*1e3:Entry$",
		   Form("MyStat>8000 && eigID==%d",j),"*");
    TGraph *g1 = new TGraph(eig_tree->GetSelectedRows(),
			    eig_tree->GetV2(),eig_tree->GetV1());
    g1->SetMarkerStyle(20);
    if(j+1==5 || j+1==9)
      color_offset ++;
    g1->SetLineColor(j+1+color_offset);
    g1->SetMarkerColor(j+1+color_offset);
    g1->SetMarkerStyle(20);
    leg->AddEntry(g1,Form("Eigenvector # %d",j));
    mg->Add(g1,"lp");
  }
  mg->Draw("ALP");
  mg->SetTitle("#sqrt{ EigenValue} (um) vs minrun ");
  mg->GetYaxis()->SetTitle(" um ");
  mg->GetXaxis()->SetTitle(" minirun # ");
  leg->Draw("same");
  c1->Print(Form("plots/slug%d_eigval.pdf",slug));
  
  input->Close();
  
}

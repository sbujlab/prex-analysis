void DrawEigVecFromDance(Int_t kSwitch, Int_t run){
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  vector<TString> IVlist;  
  TString filename_tag;
  TString run_title = Form("Run%d",run);
  if(kSwitch==0){  // 5 bpm
    vector< TString > IVlist2 =  {"bpm4aX","bpm4eX","bpm11X12X",
				  "bpm4aY","bpm4eY"};
      IVlist = IVlist2;
    filename_tag = "5bpm";
  }else if(kSwitch==1){  // all bpm
    vector< TString > IVlist2 =  {"bpm4aX","bpm4eX","bpm1X",
				  "bpm11X","bpm12X","bpm16X",
				  "bpm4aY","bpm4eY","bpm1Y",
				  "bpm11Y","bpm12Y","bpm16Y"};
    IVlist = IVlist2;
    filename_tag = "allbpm";
  }

  TFile *input = TFile::Open(Form("./rootfiles/mc_lagrange_%d.000.root",run));
  TTree *eig_tree;
  if(kSwitch==0){
    eig_tree = (TTree*)input->Get("mini_eig");
  }else if(kSwitch==1){
    eig_tree = (TTree*)input->Get("mini_eigall");
  }

  eig_tree->AddFriend("mini");
  eig_tree->AddFriend("slope",Form("./rootfiles/dit_eigslopes_%s_run%d.root",
  				   filename_tag.Data(),run));
  
  Int_t nBPM = IVlist.size();
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Print(Form("plots/run%d_vector_%s.pdf[",run,filename_tag.Data()));
  // Eigen Values
  TMultiGraph *mgval = new TMultiGraph();
  TLegend *legval = new TLegend(0.9,0.6,0.99,0.9);
  Int_t color_offset=0;
  vector<Int_t> fRun;
  vector<Int_t> fMini;
  eig_tree->Draw("run:mini","","goff");
  Double_t *yval = eig_tree->GetV1();
  Double_t *xval = eig_tree->GetV2();
  Int_t nmini = eig_tree->GetSelectedRows();
  for(int i =0;i<nmini;i++){
    fRun.push_back(yval[i]);
    fMini.push_back(xval[i]);
  }

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
    legval->AddEntry(g1,Form("Eigenvector #%d",j),"lp");
    mgval->Add(g1,"lp");
  } // component loop
  c2->cd();
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
  hmgval->SetTitle(run_title+": #sqrt{Eigenvalues} (um); Run Number; um");
  hmgval->GetXaxis()->SetTitleOffset(1.25);
  legval->Draw("same");
  c2->Print(Form("plots/run%d_vector_%s.pdf",run,filename_tag.Data()));
  c2->Clear("D");
  //Eigen vectors
  for(int i=0;i<nBPM;i++){
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.9,0.6,0.99,0.9);
    Int_t color_offset=0;
    for(int j=0;j<nBPM;j++){
      Int_t npt =  eig_tree->Draw(Form("sign%d*evMon%d_%s:Entry$",i,i,IVlist[j].Data()),
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
    hmg->SetTitle(run_title+Form(": EigenVector #%d;Run Number ;Magnitude",i));
    hmg->GetXaxis()->SetTitleOffset(1.25);
    leg->Draw("same");
    c2->Print(Form("plots/run%d_vector_%s.pdf",run,filename_tag.Data()));
  } // eigenvalue loop
  c2->Print(Form("plots/run%d_vector_%s.pdf]",run,filename_tag.Data()));
  input->Close();
}

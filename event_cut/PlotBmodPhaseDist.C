void PlotBmodPhaseDist(){
  TFile *file2 = TFile::Open("bmod_merged_all.root");
  TTree* bmod_tree = (TTree*)file2->Get("bmod");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,700);
  c1->cd();
  // c1->SetLogy();
  c1->Print("bmod_phase_hist_all_run.pdf[");
  vector<TString> fChannelX={"bpm4eX","bpm4aX","bpm11X","bpm12X"};
  vector<TString> fChannelY={"bpm4aY","bpm4eY"};
  Int_t count=0;
  for(int icoil=1;icoil<=7;icoil++){
    vector<TString> fChannel;
    if(icoil%2==1)
      fChannel = fChannelX;
    else
      fChannel = fChannelY;
    
    Int_t nch  = fChannel.size();
    TLegend *leg = new TLegend(1.0,0.9,0.9,0.7);
    TLegend *leg2 = new TLegend(1.0,0.9,0.9,0.7);
    TMultiGraph *mg = new TMultiGraph();
    THStack *hs = new THStack("hs",Form("coil%d Modulation phase",icoil));
    for(int ich=0;ich<nch;ich++){
      TString device_name = fChannel[ich];
      TH1D *hstd = new TH1D(Form("h%d",count),"",360,-110,250);
      bmod_tree->Draw( Form("(%s.phase/deg>250 ? %s.phase-2.0*TMath::Pi() : %s.phase)/deg>>h%d",
      			    device_name.Data(),device_name.Data(),device_name.Data(),count),
      		       Form("%s.ndata>100 && coil==%d ",device_name.Data(),icoil));
      hstd->SetLineColor(kBlack+ich);
      hs->Add(hstd);
      // hstd->Draw();
      count ++;
      leg->AddEntry(hstd,device_name,"l");
      
      bmod_tree->Draw( Form("(%s.phase/deg>250 ? %s.phase-2.0*TMath::Pi() : %s.phase)/deg:run",
      			    device_name.Data(),device_name.Data(),device_name.Data()),
      		       Form("%s.ndata>100 && coil==%d",device_name.Data(),icoil),"*");
      TGraph* g1 = new TGraph(bmod_tree->GetSelectedRows(),bmod_tree->GetV2(),bmod_tree->GetV1());
      g1->SetMarkerStyle(7);
      g1->SetMarkerColor(kBlack+ich);
      g1->SetLineColor(kBlack+ich);
      mg->Add(g1,"P");
      leg2->AddEntry(g1,device_name,"p");

    } // end of ich loop
    hs->Draw("nostack");
    leg->Draw("same");
    c1->Print("bmod_phase_hist_all_run.pdf");

    mg->Draw("AP");
    mg->SetTitle(Form("coil%d Modulation Phase;run number; Phase (degree)",icoil));
    leg2->Draw("same");
    c1->Print("bmod_phase_hist_all_run.pdf");
  }
  
  c1->Print("bmod_phase_hist_all_run.pdf]");

}

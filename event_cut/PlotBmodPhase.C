void PlotBmodPhase(){
  TFile *file2 = TFile::Open("bmod_merged_all.root");
  TTree* bmod_tree = (TTree*)file2->Get("bmod");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,700);
  c1->cd();
  c1->Print("bmod_phase_all_run.pdf[");
  c1->SetGridy();
  vector<Int_t> coilH= { 1,3,5,7};
  vector<Int_t> coilV= { 2,4,6};
  vector<Int_t> coilID=coilH;;

  vector<TString> fChannelArray={"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm11X","bpm12X","usl","usr"};
  Int_t nch = fChannelArray.size();
  
  for(int ich=0;ich<nch;ich++){
    TString device_name = fChannelArray[ich];
    if(device_name.Contains("X"))
      coilID= coilH;
    else if(device_name.Contains("Y"))
      coilID= coilV;


    TString unit="deg";
    TString yaxis_title=" Phase (degree)";
    // if(device_name.Contains("bpm")){
    //   unit = "um";
    //   yaxis_title = "   Modulation Phase (um)";
    // }
    // else
    //   unit = device_name+".offset";
    TString det_cut = "";
    if(device_name.Contains("s")){
      det_cut = "&&" +device_name+".offset>0.02";
    }
  
    TMultiGraph* mg =new TMultiGraph();
    TLegend *leg = new TLegend(1.0,0.9,0.9,0.7);
    for(int i=0;i<coilID.size();i++){
      bmod_tree->Draw(Form("(%s.phase/deg>250 ? %s.phase-2*TMath::Pi() : %s.phase ) /%s:run",
			   device_name.Data(),device_name.Data(),device_name.Data(),unit.Data()),
		      Form("%s.ndata>100  && coil==%d",device_name.Data(),coilID[i])+det_cut,"*");
      TGraph* g3 = new TGraph(bmod_tree->GetSelectedRows(),bmod_tree->GetV2(),bmod_tree->GetV1());
      g3->SetMarkerStyle(7);
      g3->SetMarkerColor(kBlack+i);
      g3->SetLineColor(kBlack+i);
      mg->Add(g3,"P");
      leg->AddEntry(g3,Form("coil %d",coilID[i]));
    }// end of coil loop
    mg->Draw("AP");
    mg->SetTitle(device_name+yaxis_title+";run number;"+yaxis_title);
    leg->Draw("same");
    c1->Print("bmod_phase_all_run.pdf");

  } // end of ich loop
  
  c1->Print("bmod_phase_all_run.pdf]");

}

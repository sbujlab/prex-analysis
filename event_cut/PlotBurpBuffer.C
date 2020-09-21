void PlotBurpBuffer(){
  TFile *file2 = TFile::Open("bmod_merged_all_new.root");
  TTree* bmod_tree = (TTree*)file2->Get("bmod");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,700);
  TPad pad1("pad1","",0.0,1.0,0.7,0.5);
  TPad pad2("pad2","",0.0,0.5,0.7,0.0);
  TPad pad3("pad3","",0.7,1.0,1.0,0.5);
  TPad pad4("pad4","",0.7,0.5,1.0,0.0);
  c1->cd();
  pad1.Draw("same");
  pad2.Draw("same");
  pad3.Draw("same");
  pad4.Draw("same");
  pad1.SetGridy();
  pad2.SetGridy();
  c1->Print("bmod_burpbuff_all_cycle.pdf[");
  
  vector<Int_t> coilH= { 1,3,5,7};
  vector<Int_t> coilV= { 2,4,6,};
  vector<Int_t> coilID=coilH;

  vector<TString> fChannelArray={"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm11X","bpm12X"};
  Int_t nch = fChannelArray.size();
  
  for(int ich=0;ich<nch;ich++){
    TString device_name = fChannelArray[ich];
    TString coil_select;
    if(device_name.Contains("X")){
      coilID= coilH;
      coil_select = " && coil%2==1";
    }
    else if(device_name.Contains("Y")){
      coilID= coilV;
      coil_select = " && coil%2==0";
    }
    TString unit="1";
    TString yaxis_title="";
    if(device_name.Contains("bpm")){
      unit = "um";
      yaxis_title = " : Scaled Amplitude + 3-sigma burp limit (um)";
    }
    else
      unit = device_name+".offset";
    TString det_cut = "";
    if(device_name.Contains("s")){
      det_cut = "&&" +device_name+".offset>0.02";
      yaxis_title = "   Modulation Amplitude (fraction yield)";
    }

    pad4.cd();
    bmod_tree->Draw(Form("%s.bkg/%s",
			 device_name.Data(),unit.Data()),
		    Form("%s.ndata>200 && run!=3763 && run!=4303",device_name.Data())+coil_select,"");
    TH1D *h1dptr = (TH1D*)gPad->FindObject("htemp");
    
    double rms_hi =10*(h1dptr->GetRMS());
    det_cut += Form("&& %s.bkg/%s<=%f",device_name.Data(),unit.Data(),rms_hi);
    TH1D hjitter("hjitter","",rms_hi,0,rms_hi);
    bmod_tree->Draw(Form("%s.bkg/%s>>hjitter",
			 device_name.Data(),unit.Data()),
		    Form("%s.ndata>200 && run!=3763 && run!=4303",device_name.Data())+det_cut+coil_select,"");

    TMultiGraph* mg =new TMultiGraph();
    TLegend *leg = new TLegend(1.0,0.9,0.9,0.7);
    pad1.cd();
    for(int i=0;i<coilID.size();i++){
      bmod_tree->Draw(Form("((run> 3877 ? 1.31 : 0.9)*%s + 3*%s.bkg)/%s:run",
			   device_name.Data(),device_name.Data(),unit.Data()),
		      Form("%s.ndata>200  && coil==%d && run!=3763 && run!=4303",device_name.Data(),coilID[i])+det_cut,"*");
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

    pad3.cd();
    bmod_tree->Draw(Form("((run> 3877 ? 1.31 : 0.9)*%s + 3*%s.bkg)/%s",
			 device_name.Data(),device_name.Data(),unit.Data()),
		    Form("%s.ndata>200 && run!=3763 && run!=4303",device_name.Data())+det_cut+coil_select,"");
    h1dptr = (TH1D*)gPad->FindObject("htemp");
    h1dptr->SetTitle(" ");

    Double_t nBins = h1dptr->GetNbinsX();
    Double_t BinMax = h1dptr->GetBinCenter(nBins);
    Double_t bin_size = BinMax / nBins;
    Double_t npass =h1dptr->Integral(0,350/bin_size);
    Double_t npass2 =h1dptr->Integral(0,600/bin_size);
    Double_t rate = npass/h1dptr->GetEntries();
    Double_t rate2 = npass2/h1dptr->GetEntries();
    TText label (0.3,0.65, Form( "Expect %.0f %% within 0.35 mm" , rate*100));
    TText label2 (0.3,0.7, Form( "Expect %.0f %% within 0.60 mm" , rate2*100));
    label.SetNDC();
    label.Draw("same");
    label2.SetNDC();
    label2.Draw("same");

    TMultiGraph* mgrms =new TMultiGraph();
    TLegend *legrms = new TLegend(1.0,0.9,0.9,0.7);
    pad2.cd();
    for(int i=0;i<coilID.size();i++){
      bmod_tree->Draw(Form("(%s.bkg)/%s:run",
			   device_name.Data(),unit.Data()),
		      Form("%s.ndata>200  && coil==%d && run!=3763 && run!=4303",device_name.Data(),coilID[i])+det_cut,"*");
      TGraph* g3 = new TGraph(bmod_tree->GetSelectedRows(),bmod_tree->GetV2(),bmod_tree->GetV1());
      g3->SetMarkerStyle(7);
      g3->SetMarkerColor(kBlack+i);
      g3->SetLineColor(kBlack+i);
      mgrms->Add(g3,"P");
      legrms->AddEntry(g3,Form("coil %d",coilID[i]));
    }// end of coil loop

    mgrms->Draw("AP");
    yaxis_title = " : beam random jitter (um)";
    mgrms->SetTitle(device_name+yaxis_title+";run number;"+yaxis_title);
    legrms->Draw("same");


    c1->Print("bmod_burpbuff_all_cycle.pdf");
    c1->Clear("D");
  } // end of ich loop
  c1->Print("bmod_burpbuff_all_cycle.pdf]");

}

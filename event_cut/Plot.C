void Plot(TString user_cut=""){
  TFile *file1 = TFile::Open("BPMCheck_shit.root");
  TTree* T1 = (TTree*)file1->Get("T");
  TCanvas *c1 = new TCanvas("c1","c1",1200,700);
  c1->cd();
  c1->Print("event_cut_all_run.pdf[");
  c1->SetGridy();
  vector<TString> fChannelArray={"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm12X","usl","usr"};
  Int_t nch = fChannelArray.size();
  
  for(int ich=0;ich<nch;ich++){
    TString device_name = fChannelArray[ich];
    TMultiGraph* mg = new TMultiGraph();

    TString unit="1";
    TString yaxis_title="";
    if(device_name.Contains("bpm")){
      unit = "um";
      yaxis_title = " yield (um)";
    }
    else
      unit = "yield_"+device_name;
    TString det_cut = "";
    if(device_name.Contains("s")){
      det_cut = "&& yield_"+device_name+">0.02";
      yaxis_title = " (fraction yield)";
    }
  
    // T1->Draw(Form("(yield_%s.ul-yield_%s.ll)/2/%s:run",device_name.Data(),device_name.Data(),unit.Data()),
    // 	   Form("yield_%s.num_samples>0",device_name.Data())+det_cut+user_cut,"*");
    // TGraph* g1 = new TGraph(T1->GetSelectedRows(),T1->GetV2(),T1->GetV1());
    // g1->SetName("buff");
    // g1->SetMarkerStyle(7);
    // g1->SetMarkerColor(kMagenta);
    // g1->SetLineColor(kMagenta);
    // mg->Add(g1,"P");

    // T1->Draw(Form("yield_%s/%s:run:yield_%s.sigma/%s",
    // 		device_name.Data(),unit.Data(),device_name.Data(),unit.Data()),
    // 	   Form("yield_%s.num_samples>0",device_name.Data())+det_cut+user_cut,"*");
  
    T1->Draw(Form("0:run:3*yield_%s.fwhm/2.335/%s",
		  device_name.Data(),unit.Data()),
	     Form("yield_%s.num_samples>0",device_name.Data())+det_cut+user_cut,"*");

    TGraphErrors* g1 = new TGraphErrors(T1->GetSelectedRows(),T1->GetV2(),T1->GetV1(),0,T1->GetV3());
    g1->SetName("buff");
    // g1->SetMarkerStyle(7);
    g1->SetMarkerColor(kBlack);
    g1->SetLineColor(kOrange);
    g1->SetFillColor(kOrange);
    g1->SetFillStyle(3001);
    mg->Add(g1,"P");
  
    T1->Draw(Form("(yield_%s.ll-yield_%s)/%s:run",
		  device_name.Data(),device_name.Data(),unit.Data()),
	     Form("yield_%s.num_samples>0",device_name.Data())+det_cut+user_cut,"*");
    TGraph* gll = new TGraph(T1->GetSelectedRows(),T1->GetV2(),T1->GetV1());
    gll->SetName("buff");
    gll->SetMarkerStyle(6);
    gll->SetMarkerColor(kRed);
    gll->SetLineColor(kRed);
    mg->Add(gll,"P");

    T1->Draw(Form("(yield_%s.ul - yield_%s)/%s:run",
		  device_name.Data(),device_name.Data(),unit.Data()),
	     Form("yield_%s.num_samples>0",device_name.Data())+det_cut+user_cut,"*");
    TGraph* gul = new TGraph(T1->GetSelectedRows(),T1->GetV2(),T1->GetV1());
    gul->SetName("buff");
    gul->SetMarkerStyle(6);
    gul->SetMarkerColor(kBlue);
    gul->SetLineColor(kBlue);
    mg->Add(gul,"P");

    // T1->Draw(Form("yield_%s/%s:run",device_name.Data(),unit.Data()),
    // 	   Form("yield_%s.num_samples>0",device_name.Data())+det_cut+user_cut,"*");
    // TGraph* g1 = new TGraph(T1->GetSelectedRows(),T1->GetV2(),T1->GetV1());
    // g1->SetName("buff");
    // g1->SetMarkerStyle(7);
    // g1->SetMarkerColor(kBlack);
    // g1->SetLineColor(kBlack);
    // mg->Add(g1,"P");

    // T1->Draw(Form("3*(yield_%s.sigma)/%s:run",device_name.Data(),unit.Data()),
    // 	   Form("yield_%s.num_samples>0",device_name.Data())+det_cut+user_cut,"*");
    // TGraph* g2 = new TGraph(T1->GetSelectedRows(),T1->GetV2(),T1->GetV1());
    // g2->SetMarkerStyle(7);
    // g2->SetMarkerColor(kOrange);
    // mg->Add(g2,"P");
    mg->Draw("AP");
    mg->SetTitle(fChannelArray[ich]+yaxis_title+";run number;"+yaxis_title);
    // if(kZoom){
    //   if(fChannelArray[ich].Contains("bpm4aX"))
    // 	mg->GetYaxis()->SetRangeUser(-100,100);
    //   else if(fChannelArray[ich].Contains("bpm"))
    // 	mg->GetYaxis()->SetRangeUser(-500,500);
    //   else
    // 	mg->GetYaxis()->SetRangeUser(-0.05,0.05);
    // }
    TLegend *leg = new TLegend(1.0,0.9,0.9,0.7);
    leg->AddEntry(g1," 3 sigma ");
    leg->AddEntry(gul," upper edge ");
    leg->AddEntry(gll," lower edge ");
    leg->Draw("same");
    c1->Print("event_cut_all_run.pdf");
  } // ich loop
  c1->Print("event_cut_all_run.pdf]");
}

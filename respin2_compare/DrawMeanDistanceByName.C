  
vector< pair<TString,TString> > fSets ={ {"dit","reg_5bpm"},
					 {"dit","lagr_6bpm"},
					 {"dit","lagr_all"},
					 {"dit","reg_all"},
					 {"reg_5bpm","reg_6bpm"},
					 {"reg_5bpm","reg_all"},
					 {"lagr_6bpm","reg_6bpm"},
					 {"lagr_all","reg_all"},
					 {"lagr_6bpm","lagr_all"},
					 {"lagr_all_trunc","lagr_all"}};


void DrawMeanDistanceByName(TString det_name="us_dd"){
  vector<Double_t> fDeltaMu_slug;
  vector<Double_t> fSigmaDeltaMu_slug;
  vector<Double_t> fChiSquare_slug;
  vector<Int_t> fNDF_slug;
  vector<Double_t> fProb_slug;
  
  vector<Double_t> fDeltaMu_slugfix;
  vector<Double_t> fSigmaDeltaMu_slugfix;
  vector<Double_t> fChiSquare_slugfix;
  vector<Int_t> fNDF_slugfix;
  vector<Double_t> fProb_slugfix;
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);
  TCanvas *c1 = new TCanvas("c1","c1",1400,600);
  // c1->SetGridx();
  // c1->SetGridy();
  c1->cd();
  TPad *pad1 = new TPad("","",0,0,0.66,1.0);
  TPad *pad2 = new TPad("","",0.66,0,1.0,1.0);
  pad1->Draw();
  pad2->Draw();
  pad1->SetGridx();
  pad1->SetGridy();
  TCanvas *c2 = new TCanvas("c2","c2",1400,600);
  c2->SetGridy();
  c2->SetGridx();
  auto iter= fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    TF1 *fp0fix = new TF1("fp0fix","[0]",-100,100);
    TF1 *fp0free = new TF1("fp0free","[0]",-100,100);
    fp0fix->FixParameter(0,0);
    
    TFile *input =  TFile::Open(Form("./rootfiles/prex_minidiff_%s_%s.root",set1.Data(),set2.Data()));
    TTree *slug_tree = (TTree*)input->Get("slug");
    c2->cd();
    TMultiGraph *mg_dist = new TMultiGraph();
    TLegend *leg_dist = new TLegend(0.99,0.9,0.9,0.7);
    slug_tree->Draw(Form("dist_%s:slug+arm_flag/6",det_name.Data()),"arm_flag==0","goff");
    TGraph* gdist = new TGraph(slug_tree->GetSelectedRows(),slug_tree->GetV2(),slug_tree->GetV1());
    gdist->SetLineColor(kBlack);
    gdist->SetMarkerColor(kBlack);
    gdist->SetMarkerStyle(20);
    mg_dist->Add(gdist,"lp");
    leg_dist->AddEntry(gdist,"both-arm","p");
    
    mg_dist->Draw("A");
    mg_dist->SetTitle(Form(" %s vs %s : #sqrt{ #LT(#Delta A )^{2}#GT } (ppb) vs Slug; Slug; #sqrt{ #LT(#Delta A )^{2}#GT } (ppb) ",
			   set1.Data(),set2.Data()));
    mg_dist->GetXaxis()->SetTitleSize(0.07);
    mg_dist->GetXaxis()->SetTitleOffset(0.5);
    mg_dist->GetYaxis()->SetTitleSize(0.07);
    mg_dist->GetYaxis()->SetTitleOffset(0.5);
    leg_dist->Draw("same");
    c2->SaveAs(Form("./plots/delta_A_rms_%s_%s_by_slug_%s.pdf",
		    set1.Data(),set2.Data(),det_name.Data()));


    pad1->cd();
    TMultiGraph *mg_er = new TMultiGraph();
    TLegend *leg_er = new TLegend(0.99,0.9,0.9,0.7);
    slug_tree->Draw(Form("sign*delta_%s:delta_%s_err:slug+arm_flag/6",det_name.Data(),det_name.Data()),
		    "arm_flag==0","goff");
    // pretending it is both-only...
    TGraphErrors* ger_both = new TGraphErrors(slug_tree->GetSelectedRows(),
					      slug_tree->GetV3(),slug_tree->GetV1(),
					      0,slug_tree->GetV2());
    
    ger_both->Fit("fp0free","Q");  
    double* yval = slug_tree->GetV1();
    double* yerr = slug_tree->GetV2();
    int ynpt = slug_tree->GetSelectedRows();
    TH1D *hpull_free = new TH1D("hpull_free","Pull Fit",40,-5,8);
    TH1D *hpull_fixed = new TH1D("hpull_fixed","Pull Fit",40,-5,8);
    double mean = fp0free->GetParameter(0);
    for(int i=0;i<ynpt;i++){
      hpull_free->Fill((yval[i]-mean)/yerr[i]);
      hpull_fixed->Fill(yval[i]/yerr[i]);
    }
    hpull_free->Fit("gaus","Q");
    hpull_fixed->Fit("gaus","Q");
    hpull_free->SetName("Pull");
    hpull_fixed->SetName("Pull");
    
    ger_both->SetLineColor(kBlack);
    ger_both->SetMarkerColor(kBlack);
    ger_both->SetMarkerStyle(20);
    
    mg_er->Add(ger_both,"p");
    leg_er->AddEntry(ger_both,"both-arm","p");
    pad1->cd();
    mg_er->Draw("A");
    // ger_both->Fit("fp0free","Q");  
    leg_er->Draw("same");
    mg_er->SetTitle(Form("%s -  %s vs %s :Sign Corrected  #Delta A #pm #sigma(#Delta A) (ppb) vs Slug; Slug; #Delta A (ppb) ",
			 det_name.Data(),set1.Data(),set2.Data()));
    mg_er->GetXaxis()->SetTitleSize(0.07);
    mg_er->GetXaxis()->SetTitleOffset(0.5);
    mg_er->GetYaxis()->SetTitleSize(0.07);
    mg_er->GetYaxis()->SetTitleOffset(0.5);
    double ymax = mg_er->GetYaxis()->GetXmax();
    double ymin = mg_er->GetYaxis()->GetXmin();
    mg_er->GetYaxis()->SetRangeUser(ymin,ymax+0.3*(ymax-ymin));
    gPad->Update();
    TPaveStats *ps = (TPaveStats*)ger_both->FindObject("stats");
    ps->SetX2NDC(0.9);
    pad2->cd();
    hpull_free->Draw();
    pad2->Update();
    TPaveStats *pspull_free = (TPaveStats*)hpull_free->FindObject("stats");
    pspull_free->SetOptStat(111111);
    pspull_free->SetY1NDC(0.5);
    pspull_free->SetX1NDC(0.5);
    fDeltaMu_slug.push_back(fp0free->GetParameter(0));
    fSigmaDeltaMu_slug.push_back(fp0free->GetParError(0));
    fChiSquare_slug.push_back(fp0free->GetChisquare());
    fNDF_slug.push_back(fp0free->GetNDF());
    fProb_slug.push_back(fp0free->GetProb());
    c1->SaveAs(Form("./plots/delta_A_freepar_fit_%s_%s_by_slug_%s.pdf",
		    set1.Data(),set2.Data(),det_name.Data()));
    
    pad1->cd();
    ger_both->Fit("fp0fix","Q");
    fDeltaMu_slugfix.push_back(fp0fix->GetParameter(0));
    fSigmaDeltaMu_slugfix.push_back(fp0fix->GetParError(0));
    fChiSquare_slugfix.push_back(fp0fix->GetChisquare());
    fNDF_slugfix.push_back(fp0fix->GetNDF());
    fProb_slugfix.push_back(fp0fix->GetProb());
    pad2->cd();
    hpull_fixed->Draw();
    pad2->Update();
    TPaveStats *pspull_fixed = (TPaveStats*)hpull_fixed->FindObject("stats");
    pspull_fixed->SetOptStat(111111);
    pspull_fixed->SetY1NDC(0.5);
    pspull_fixed->SetX1NDC(0.5);

    c1->SaveAs(Form("./plots/delta_A_fixpar_fit_%s_%s_by_slug_%s.pdf",
		    set1.Data(),set2.Data(),det_name.Data()));
    
    input->Close();
    iter++;
  } // end of iteration comparison sets

  cout << " ********* Slug Summary ************ " << endl;
  Int_t icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");
    cout << set1 <<  " vs " << set2 << "&";
    printf("%.1f & %.1f & %.1f & %d & %.2f\\\\ \n",
	   fDeltaMu_slug[icount],fSigmaDeltaMu_slug[icount],
	   fChiSquare_slug[icount],fNDF_slug[icount],fProb_slug[icount]);
    iter++;
    icount++;
  }
  cout << " ********* fix para " << endl;
  icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");

    cout << set1 <<  " vs " << set2 << "&";
    printf("%.0f & %.0f & %.1f & %d & %.2f \\\\ \n",
	   fDeltaMu_slugfix[icount],fSigmaDeltaMu_slugfix[icount],
	   fChiSquare_slugfix[icount],fNDF_slugfix[icount],fProb_slugfix[icount]);
    iter++;
    icount++;
  }
  
}

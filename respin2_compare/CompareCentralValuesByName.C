void CompareCentralValuesByName(TString det_name="us_dd"){

  map<TString,TString> fLinks={ {"reg_5bpm","prex_grand_average_regress_cbpm.root"},
				{"dit","prex_grand_average_dither_respin2.root"},
				{"reg_6bpm","prex_grand_average_reg_6bpm.root"},
				{"lagr_6bpm","prex_grand_average_lagr_6bpm.root"},
				{"reg_all","prex_grand_average_reg_allbpm.root"},
				{"lagr_all","prex_grand_average_lagr_allbpm.root"},
				{"lagr_all_trunc","prex_grand_average_lagr_allbpm_truncated.root"}};

  
  vector< pair<TString,TString> > fSets ={ {"dit","reg_5bpm"},
					   {"dit","lagr_6bpm"},
					   {"dit","lagr_all"},
					   {"dit","reg_all"},
					   {"reg_5bpm","reg_6bpm"},
					   {"reg_5bpm","reg_all"},
					   {"lagr_6bpm","reg_6bpm"},
					   {"lagr_all","reg_all"},
					   {"lagr_6bpm","lagr_all"},
					   {"lagr_all_trunc","lagr_all"}  };

  vector<Double_t> fDeltaMu_slug;
  vector<Double_t> fSigmaDeltaMu_slug;
  vector<Double_t> fChiSquare_slug;
  vector<Int_t> fNDF_slug;
  vector<Double_t> fProb_slug;
  
  vector<Double_t> fDeltaMu_pitt;
  vector<Double_t> fSigmaDeltaMu_pitt;
  vector<Double_t> fChiSquare_pitt;
  vector<Int_t> fNDF_pitt;
  vector<Double_t> fProb_pitt;
  
  vector<Double_t> fDeltaMu_slugfix;
  vector<Double_t> fSigmaDeltaMu_slugfix;
  vector<Double_t> fChiSquare_slugfix;
  vector<Int_t> fNDF_slugfix;
  vector<Double_t> fProb_slugfix;

  vector<Double_t> fDeltaMu_pittfix;
  vector<Double_t> fSigmaDeltaMu_pittfix;
  vector<Double_t> fChiSquare_pittfix;
  vector<Int_t> fNDF_pittfix;
  vector<Double_t> fProb_pittfix;
  
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->SetGridx();
  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->Divide(1,2);
  auto iter= fSets.begin();

  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;

    cout << set1 <<": " << fLinks[set1] << endl;
    cout << set2 <<": " << fLinks[set2] << endl;

    TString branch_name1, branch_name2;
    if(set1.Contains("dit"))
      branch_name1 = "dit_asym_"+det_name;
    else if(set1.Contains("reg"))
      branch_name1 = "reg_asym_"+det_name;
    else if(set1.Contains("lagr"))
      branch_name1 = "lagr_asym_"+det_name;

    if(set2.Contains("dit"))
      branch_name2 = "dit_asym_"+det_name;
    else if(set2.Contains("reg"))
      branch_name2 = "reg_asym_"+det_name;
    else if(set2.Contains("lagr"))
      branch_name2 = "lagr_asym_"+det_name;

    TFile* file1 = TFile::Open(fLinks[set1]);
    TTree *slug_tree1 = (TTree*)file1->Get("slug");
    slug_tree1->AddFriend(set2+"=slug",fLinks[set2]);
    TTree *pitt_tree1 = (TTree*)file1->Get("pitt");
    pitt_tree1->AddFriend(set2+"=pitt",fLinks[set2]);
    
    TF1 *ffix = new TF1("ffix","[0]",-100,100);
    ffix->FixParameter(0,0);
    TF1 *ffloat = new TF1("ffloat","[0]",-100,100);
    c1->cd();
    slug_tree1->Draw(Form("sign*(%s.mean-%s.%s.mean)*1e9:sqrt(abs(%s.error**2-%s.%s.error**2))*1e9:slug",
			  branch_name1.Data(),set2.Data(),branch_name2.Data(),
			  branch_name1.Data(),set2.Data(),branch_name2.Data()),
		     "arm_flag==0","goff");

    TGraphErrors *ger_slug = new TGraphErrors(slug_tree1->GetSelectedRows(),
					      slug_tree1->GetV3(),slug_tree1->GetV1(),
					      0,slug_tree1->GetV2());
    ger_slug->SetMarkerStyle(20);
    ger_slug->Draw("AP");
    ger_slug->SetTitle(set1 + " vs "+set2 + ":      sign corrected (#mu_{1} - #mu_{2}) #pm #sqrt{#sigma_{1}^{2} - #sigma_{2}^{2} }; Slug ; (ppb) "  );
    ger_slug->GetYaxis()->SetTitleSize(0.05);
    ger_slug->GetXaxis()->SetTitleSize(0.05);
    ger_slug->Fit("ffloat","Q");
    c1->SaveAs(Form("compare_%s_vs_%s_by_slug_freepar_%s.png",set1.Data(),set2.Data(),det_name.Data()) );

    fDeltaMu_slug.push_back(ffloat->GetParameter(0));
    fSigmaDeltaMu_slug.push_back(ffloat->GetParError(0));
    fChiSquare_slug.push_back(ffloat->GetChisquare());
    fNDF_slug.push_back(ffloat->GetNDF());
    fProb_slug.push_back(ffloat->GetProb());
    
    ger_slug->Fit("ffix","Q");
    fDeltaMu_slugfix.push_back(ffix->GetParameter(0));
    fSigmaDeltaMu_slugfix.push_back(ffix->GetParError(0));
    fChiSquare_slugfix.push_back(ffix->GetChisquare());
    fNDF_slugfix.push_back(ffix->GetNDF());
    fProb_slugfix.push_back(ffix->GetProb());
    
    c1->SaveAs(Form("compare_%s_vs_%s_by_slug_fixpar_%s.png",set1.Data(),set2.Data(),det_name.Data()) );

    c2->cd(1);
    gPad->SetGridx();
    slug_tree1->Draw(Form("%s.rms*1e6:slug",branch_name1.Data()),"arm_flag==0","goff");
    
    TGraph *grms1 = new TGraph(slug_tree1->GetSelectedRows(),
			       slug_tree1->GetV2(),slug_tree1->GetV1());
    grms1->SetMarkerStyle(20);
    grms1->SetMarkerColor(kBlue);
    grms1->SetLineColor(kBlue);

    slug_tree1->Draw(Form("%s.%s.rms*1e6:slug",set2.Data(),branch_name2.Data()),"arm_flag==0","goff");
    TGraph *grms2 = new TGraph(slug_tree1->GetSelectedRows(),
			       slug_tree1->GetV2(),slug_tree1->GetV1());
    grms2->SetMarkerStyle(20);
    grms2->SetMarkerColor(kRed);
    grms2->SetLineColor(kRed);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(grms1,"lp");
    mg->Add(grms2,"lp");
    mg->Draw("A");
    mg->SetTitle(set1 + " vs "+set2 + ":  RMS(ppm); Slug ; (ppm) "  );
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetTitleSize(0.05);
    TLegend* legrms = new TLegend(0.9,0.7,0.99,0.9);
    legrms->AddEntry(grms1,set1,"lp");
    legrms->AddEntry(grms2,set2,"lp");
    legrms->Draw("same");
    c2->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    slug_tree1->Draw(Form("sqrt(abs(%s.rms**2-%s.%s.rms**2))*1e6:slug",
			  branch_name1.Data(),set2.Data(),branch_name2.Data()),
		     "arm_flag==0","goff");
    TGraph *grms_diff = new TGraph(slug_tree1->GetSelectedRows(),
				   slug_tree1->GetV2(),slug_tree1->GetV1());

    grms_diff->SetMarkerStyle(20);
    grms_diff->Draw("ALP");
    grms_diff->SetTitle(set1 + " vs "+set2 + ": #sqrt{RMS_{1}^{2} - RMS_{2}^{2} }; Slug ; (ppm) "  );
    grms_diff->GetYaxis()->SetTitleSize(0.05);
    grms_diff->GetXaxis()->SetTitleSize(0.05);
    c2->SaveAs(Form("rms_%s_vs_%s_by_slug_%s.png",set1.Data(),set2.Data(),det_name.Data()) );
    
    /////////////////// PITT 
    c1->cd();
    pitt_tree1->Draw(Form("(%s.mean-%s.%s.mean)*1e9:sqrt(abs(%s.error**2-%s.%s.error**2))*1e9:pitt",
			  branch_name1.Data(),set2.Data(),branch_name2.Data(),
			  branch_name1.Data(),set2.Data(),branch_name2.Data()),"","goff");

    TGraphErrors *ger_pitt = new TGraphErrors(pitt_tree1->GetSelectedRows(),
					      pitt_tree1->GetV3(),pitt_tree1->GetV1(),
					      0,pitt_tree1->GetV2());
    ger_pitt->SetMarkerStyle(20);
    ger_pitt->Draw("AP");
    ger_pitt->SetTitle(set1 + " vs "+set2 + ":       (#mu_{1} - #mu_{2}) #pm #sqrt{#sigma_{1}^{2} - #sigma_{2}^{2} }; Pitt ; (ppb) "  );
    ger_pitt->GetYaxis()->SetTitleSize(0.05);
    ger_pitt->GetXaxis()->SetTitleSize(0.05);
    ger_pitt->Fit("ffloat","Q");
    
    fDeltaMu_pitt.push_back(ffloat->GetParameter(0));
    fSigmaDeltaMu_pitt.push_back(ffloat->GetParError(0));
    fChiSquare_pitt.push_back(ffloat->GetChisquare());
    fNDF_pitt.push_back(ffloat->GetNDF());
    fProb_pitt.push_back(ffloat->GetProb());
    
    c1->SaveAs(Form("compare_%s_vs_%s_by_pitt_freepar_%s.png",set1.Data(),set2.Data(),det_name.Data()) );

    ger_pitt->Fit("ffix","Q");
    fDeltaMu_pittfix.push_back(ffix->GetParameter(0));
    fSigmaDeltaMu_pittfix.push_back(ffix->GetParError(0));
    fChiSquare_pittfix.push_back(ffix->GetChisquare());
    fNDF_pittfix.push_back(ffix->GetNDF());
    fProb_pittfix.push_back(ffix->GetProb());
    c1->SaveAs(Form("compare_%s_vs_%s_by_pitt_fixpar_%s.png",set1.Data(),set2.Data(),det_name.Data()) );

    file1->Close();
    
    iter++;
  }

  cout << " ********* Slug Summary ************ " << endl;
  Int_t icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");

    cout << set1 <<  " vs " << set2 << "&";
    printf("%.1f & %.1f & %.1f & %d  & %.2f\\\\ \n",
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

  cout << " ********* Pitt Summary ************ " << endl;
  icount=0;
  iter = fSets.begin();
  while(iter!= fSets.end()){
    TString set1 = (*iter).first;
    TString set2 = (*iter).second;
    set1.ReplaceAll("_","\\_");
    set2.ReplaceAll("_","\\_");

    cout << set1 <<  " vs " << set2 << "&";
    printf("%.1f & %.1f & %.1f & %d & %.2f \\\\ \n",
	   fDeltaMu_pitt[icount],fSigmaDeltaMu_pitt[icount],
	   fChiSquare_pitt[icount],fNDF_pitt[icount],fProb_pitt[icount]);
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
    printf("%.0f & %.0f & %.1f & %d & %.2f\\\\ \n",
	   fDeltaMu_pittfix[icount],fSigmaDeltaMu_pittfix[icount],
	   fChiSquare_pittfix[icount],fNDF_pittfix[icount],fProb_pittfix[icount]);
    iter++;
    icount++;
  }

}

void Compare(Int_t run){
  TFile *input = TFile::Open(Form("rootfiles/mc_lagrange_%d.000.root",run));
  vector<TString> label={"raw","ideal","dither","reg5bpm",
			 "regall","lagr_all","truncated_lagr_all",
			 "wrong_lagr_all","wrong dither"};

  vector<TString> tree_name={"mini","mini","mini_dit","mini_reg",
			     "mini_regall","mini_lagrall","mini_eiglagrall_tr",
			     "mini_lagr_wrg","mini_dit_false"};

  vector<TString> channel_name={"asym_us_avg","asym_us_avg_ideal",
				"dit_asym_us_avg","reg_asym_us_avg",
				"reg_asym_us_avg","lagr_asym_us_avg",
				"lagr_asym_us_avg",
				"lagr_asym_us_avg","dit_asym_us_avg"};

  vector<Double_t> fMean;
  vector<Double_t> fError;
  vector<Double_t> fRMS;
  vector<Double_t> fChiSquare;
  vector<Int_t> fNDF;
  
  Int_t color_code[]={1,2,3,4,6,7,9,28,39};
  TMultiGraph* mg = new TMultiGraph();
  TLegend *leg = new TLegend(0.9,0.7,0.99,0.95);
  Int_t nplot = label.size();
  TCanvas* c1 = new TCanvas("c1","c1",1200,600);
  c1->Print(Form("mc_compare_run%d.pdf[",run));
  gStyle->SetOptFit(1);
  Int_t color_offset=0;
  for(int i=0;i<nplot;i++){
    TTree *tree = (TTree*)input->Get(tree_name[i]);
    tree->Draw(Form("%s.rms*1e6:Entry$+1",channel_name[i].Data()),"","goff");
    
    TGraph *grms = new TGraph(tree->GetSelectedRows(),
			      tree->GetV2(),tree->GetV1());
    
    grms->SetMarkerStyle(20);
    grms->SetMarkerColor(color_code[i]);
    grms->SetLineColor(color_code[i]);
    if(i!=0){
      mg->Add(grms,"lp");
      leg ->AddEntry(grms,label[i]);
    }
    tree->Draw(Form("%s.rms*1e6",channel_name[i].Data()),"","goff");
    TH1D *h1d = (TH1D*)gDirectory->FindObject("htemp");
    fRMS.push_back(h1d->GetMean());
  }
  mg->Draw("A");
  leg->Draw("same");
  mg->SetTitle("Monte Carlo Runs RMS width (ppm); minirun ; RMS(ppm)");
  c1->Print(Form("mc_compare_run%d.pdf",run));

  for(int i=0;i<nplot;i++){

    TTree *tree = (TTree*)input->Get(tree_name[i]);
    tree->Draw(Form("%s*1e9:%s.err*1e9:Entry$+1",channel_name[i].Data(),channel_name[i].Data())
	       ,"","goff");
    TGraphErrors *ger = new TGraphErrors(tree->GetSelectedRows(),
					 tree->GetV3(),tree->GetV1(),0,tree->GetV2());
    ger->SetMarkerStyle(20);
    ger->Draw("AP");
    ger->Fit("pol0","Q");
    ger->SetTitle(Form("%s: asym_us_avg(ppb); minirun ; (ppb)",label[i].Data()));
    TF1 *f1ptr = ger->GetFunction("pol0");
    fMean.push_back(f1ptr->GetParameter(0));
    fError.push_back(f1ptr->GetParError(0));
    fChiSquare.push_back(f1ptr->GetChisquare());
    fNDF.push_back(f1ptr->GetNDF());
    c1->Print(Form("mc_compare_run%d.pdf",run));
  }
  c1->Print(Form("mc_compare_run%d.pdf]",run));

  TString file_name = Form("report/run%d_summary.tex",run);
  FILE *tex_file  = fopen(file_name.Data(),"w");
  fprintf(tex_file," \\documentclass{article} \n " );
  fprintf(tex_file," \\Large\n " );
  fprintf(tex_file," \\usepackage[landscape,margin=0.5in]{geometry}\n");
  fprintf(tex_file," \\usepackage{float}\n");
  fprintf(tex_file," \\begin{document}\n " );

  fprintf(tex_file," \\begin{table}[H] \n " );
  fprintf(tex_file," \\centering \n " );
  fprintf(tex_file," \\begin{tabular}{c|c|c|c|c|c}  \n " );
  fprintf(tex_file,"  & Mean(ppb) & Error(ppb) & RMS(ppm) & $\\chi^2$ & NDof  \\\\  \n ");
  fprintf(tex_file," \\hline \n " );

  for(int i=0;i<nplot;i++){
    TString test_name = label[i];
    test_name.ReplaceAll("_"," ");
    fprintf(tex_file,"%s & %.0f & %.0f & %.1f  & %.1f &%d \\\\ \n ",
	    test_name.Data(), fMean[i], fError[i], fRMS[i], fChiSquare[i], fNDF[i]);
  }
  fprintf(tex_file," \\hline \n " );
  fprintf(tex_file," \\end{tabular}  \n " );
  fprintf(tex_file," \\caption{Run Averages}  \n " );
  fprintf(tex_file," \\end{table}  \n " );
  fclose(tex_file);
}

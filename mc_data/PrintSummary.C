void PrintSummary(Int_t run){
  TString file_name = Form("report/run%d_summary.tex",run);
  FILE *tex_file  = fopen(file_name.Data(),"a");

  TFile* input = TFile::Open(Form("rootfiles/mc_lagrange_%d.000.root",run));
  TTree* mini_tree = (TTree*)input->Get("mini");
  TTree* eig_tree = (TTree*)input->Get("mini_eigall");
  TTree* lagr_tree = (TTree*)input->Get("mini_lagrall");
  TTree* false_lagr_tree = (TTree*)input->Get("mini_lagr_wrg");
  TTree* reg_tree = (TTree*)input->Get("mini_regall");
  vector<TTree*> fTreeArray={reg_tree,lagr_tree,false_lagr_tree};
  
  TFile* input_slope1 = TFile::Open(Form("rootfiles/dit_eigslopes_allbpm_run%d.root",run));
  TTree* slope_true = (TTree*)input_slope1->Get("slope");
  slope_true->AddFriend(mini_tree);
  slope_true->AddFriend(eig_tree,"eig");
  
  TFile* input_slope2 = TFile::Open(Form("rootfiles/dit_eigslopes_allbpm_false_run%d.root",run));
  TTree* slope_false = (TTree*)input_slope2->Get("slope");
  slope_false->AddFriend(mini_tree);
  slope_false->AddFriend(eig_tree);
  eig_tree->AddFriend(slope_true,"true");
  eig_tree->AddFriend(slope_false,"false");
  
  vector<TString> fAliasArray={"mini_eigall","true","false"}
    ;
  Int_t nTree = fAliasArray.size();
  vector<TString> fBPMArray= {"bpm4aX","bpm4eX","bpm1X",
			      "bpm11X","bpm12X","bpm16X",
			      "bpm4aY","bpm4eY","bpm1Y",
			      "bpm11Y","bpm12Y","bpm16Y"};
  Int_t nBPM = fBPMArray.size();
  vector<Double_t> fBPM_Diff;
  vector<Double_t> fBPM_Diff_error;
  vector<Double_t> fBPM_RMS;
  vector<vector<Double_t> >  fBPM_slope_by_method;

  for(int i=0;i<nBPM;i++){
    TH1D *h1dptr;
    mini_tree->Draw(Form("diff_%s/(nm):diff_%s.err/(nm):Entry$",
			 fBPMArray[i].Data(),fBPMArray[i].Data()),"","goff");
    double *yval= mini_tree->GetV1();
    double *yerr= mini_tree->GetV2();
    double *xval= mini_tree->GetV3();
    TGraphErrors *mger = new TGraphErrors(mini_tree->GetSelectedRows(),
					  xval,yval,0,yerr);
    mger->Fit("pol0","Q");
    TF1* fit = mger->GetFunction("pol0");
    fBPM_Diff.push_back( fit->GetParameter(0));
    fBPM_Diff_error.push_back( fit->GetParError(0));
    mini_tree->Draw(Form("diff_%s.rms/(um)",fBPMArray[i].Data()),"","goff");
    h1dptr= (TH1D*)gDirectory->FindObject("htemp");
    fBPM_RMS.push_back(h1dptr->GetMean());
    
  }
  for(int it=0;it<nTree;it++){
    vector<Double_t> fBPM_slope;
    TH1D *h1dptr;
    for(int i=0;i<nBPM;i++){
      fTreeArray[it]->Draw(Form("us_avg_%s/(ppm/um)",fBPMArray[i].Data()),"","goff");
      h1dptr= (TH1D*)gDirectory->FindObject("htemp");
      fBPM_slope.push_back(h1dptr->GetMean());
    }
    fBPM_slope_by_method.push_back(fBPM_slope);
  }
  vector<Double_t> fEVM_Diff;
  vector<Double_t> fEVM_Diff_error;
  vector<Double_t> fEVM_RMS;
  vector<vector<Double_t> > fEVM_slope_by_method;
  vector<vector<Double_t> > fEVM_deltaA_by_method;
  vector<vector<Double_t> > fEVM_noise_by_method;

  for(int i=0;i<nBPM;i++){
    TH1D *h1dptr;
    eig_tree->Draw(Form("sign%d * diff_evMon%d/(nm):diff_evMon%d.err/(nm):Entry$",i,i,i),"","goff");
    double *yval= eig_tree->GetV1();
    double *yerr= eig_tree->GetV2();
    double *xval= eig_tree->GetV3();
    TGraphErrors *mger = new TGraphErrors(eig_tree->GetSelectedRows(),
					  xval,yval,0,yerr);
    mger->Fit("pol0","Q");
    TF1* fit = mger->GetFunction("pol0");
    fEVM_Diff.push_back(fit->GetParameter(0));
    fEVM_Diff_error.push_back(fit->GetParError(0));
    
    eig_tree->Draw(Form("diff_evMon%d.rms/(um)",i),"","goff");
    h1dptr= (TH1D*)gDirectory->FindObject("htemp");
    fEVM_RMS.push_back(h1dptr->GetMean());
  }
  
  for(int it=0;it<nTree;it++){
    vector<Double_t> fEVM_slope;
    vector<Double_t> fEVM_deltaA;
    vector<Double_t> fEVM_noise;
    TH1D *h1dptr;
    for(int i=0;i<nBPM;i++){
      eig_tree->Draw(Form("sign%d*%s.us_avg_evMon%d/(ppm/um)",
			  i,fAliasArray[it].Data(),i),"","goff");
      h1dptr= (TH1D*)gDirectory->FindObject("htemp");
      fEVM_slope.push_back(h1dptr->GetMean());
      
      eig_tree->Draw(Form("%s.us_avg_evMon%d * diff_evMon%d /ppb",
			  fAliasArray[it].Data(),i,i),"","goff");
      h1dptr= (TH1D*)gDirectory->FindObject("htemp");
      fEVM_deltaA.push_back(h1dptr->GetMean());
      
      eig_tree->Draw(Form("fabs(%s.us_avg_evMon%d * diff_evMon%d.rms)/ppm",
			  fAliasArray[it].Data(),i,i),"","goff");
      h1dptr= (TH1D*)gDirectory->FindObject("htemp");
      fEVM_noise.push_back(h1dptr->GetMean());
    }

    fEVM_slope_by_method.push_back(fEVM_slope);
    fEVM_deltaA_by_method.push_back(fEVM_deltaA);
    fEVM_noise_by_method.push_back(fEVM_noise);
  }

  
  fprintf(tex_file," \\begin{table}[H] \n " );
  fprintf(tex_file," \\centering \n " );
  fprintf(tex_file," \\begin{tabular}{c c c| c c c | c c c | c c c}  \n " );
  fprintf(tex_file," EV\\# & $\\Delta B$(nm) &  RMS(um)  & slope(ppm/um) & A$_{beam}$(ppb) & $\\sigma$(ppm) & slope(ppm/um) & A$_{beam}$(ppb) & $\\sigma$(ppm) & slope(ppm/um) & A$_{beam}$(ppb) & $\\sigma $(ppm) \\\\  \n ");
  fprintf(tex_file," \\hline \n " );
  for(int i=0;i<nBPM;i++){
    fprintf(tex_file, "%d " , i);
    fprintf(tex_file, " & %.1f $\\pm$ %.1f " , fEVM_Diff[i], fEVM_Diff_error[i]);
    fprintf(tex_file, " & %.1f " , fEVM_RMS[i]);
    for(int j=0;j<nTree;j++){
      fprintf(tex_file, "& %.2f " , fEVM_slope_by_method[j][i]);
      fprintf(tex_file, "& %.1f " , fEVM_slope_by_method[j][i] * fEVM_Diff[i]);
      fprintf(tex_file, "& %.1f" ,  fabs(fEVM_slope_by_method[j][i] * fEVM_RMS[i]) );
    }
    fprintf(tex_file, "\\\\ \n");
  }
  fprintf(tex_file," \\hline \n " );
  fprintf(tex_file," \\end{tabular}  \n " );
  fprintf(tex_file," \\caption{All BPM Corrections on US AVG: 1)Eigenvectors; 2)Regression; 3)Lagrange; 4)False Lagrange}  \n " );
  fprintf(tex_file," \\end{table} \n " );


  fprintf(tex_file," \\begin{table}[H] \n " );
  fprintf(tex_file," \\centering \n " );
  fprintf(tex_file," \\begin{tabular}{c c c| c | c | c}  \n " );
  fprintf(tex_file," BPM  & $\\Delta B$(nm) &  RMS(um)  & slope(ppm/um) & slope(ppm/um) & slope(ppm/um)\\\\  \n ");
  fprintf(tex_file," \\hline \n " );
  for(int i=0;i<nBPM;i++){
    fprintf(tex_file, "%s " , fBPMArray[i].Data());
    fprintf(tex_file, " & %.1f " , fBPM_Diff[i]);
    fprintf(tex_file, " & %.1f " , fBPM_RMS[i]);
    for(int j=0;j<nTree;j++){
      fprintf(tex_file, "& %.2f " , fBPM_slope_by_method[j][i]);
    }
    fprintf(tex_file, "\\\\ \n");
  }
  fprintf(tex_file," \\hline \n " );
  fprintf(tex_file," \\end{tabular}  \n " );
  fprintf(tex_file," \\caption{All BPM Corrections on US AVG: 1) BPMs Diff; 2)Regression; 2)Lagrange; 3)Lagrange(false)}  \n " );
  fprintf(tex_file," \\end{table} \n " );

  fprintf(tex_file," \\end{document}\n " );
  fclose(tex_file);
  

}

void GetForcedVectors(){
  TFile *input = TFile::Open("rootfiles/all_slugs.root");
  TTree *eig_tree = (TTree*)input->Get("eig");
  FILE *output = fopen("12bpm_forced_vector.txt","w");
  vector<TString> raw_bpm = {"4aX","4eX","1X","11X","12X","16X",
			     "4aY","4eY","1Y","11Y","12Y","16Y"};
  Int_t nBPM = raw_bpm.size();
  for(int iev=0;iev<nBPM;iev++){
    vector<Double_t> fWeights(nBPM);
    for(int ibpm=0;ibpm<nBPM;ibpm++){
      eig_tree->Draw(Form("evMon%d_bpm%s",iev,raw_bpm[ibpm].Data()),"","goff");
      TH1D *h1d = (TH1D*)gDirectory->FindObject("htemp");
      fWeights[ibpm] = h1d->GetMean();
    }
    Double_t norm_factor = 0.0;
    for(int ibpm=0;ibpm<nBPM;ibpm++)
      norm_factor += pow(fWeights[ibpm],2);
    norm_factor = sqrt(norm_factor);
    for(int ibpm=0;ibpm<nBPM;ibpm++)
      fWeights[ibpm] /= norm_factor;
    
    fprintf(output, "def diff_evMon%d=",iev);   
    for(int ibpm=0;ibpm<nBPM;ibpm++){
      if(fWeights[ibpm]>0)
	fprintf(output,"+");
      fprintf(output,"%.4f*diff_bpm%s",fWeights[ibpm],raw_bpm[ibpm].Data());
    }
    fprintf(output,"\n");
  }
  
  
  input->Close();
  fclose(output);

}

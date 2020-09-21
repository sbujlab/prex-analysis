void MakeMatrices(Int_t run){
  
  vector<TString> device_array={"usl","usr",
				"bpm4aX","bpm4aY","bpm4eX","bpm4eY",
				"bpm1X","bpm1Y","bpm11X","bpm11Y",
				"bpm12X","bpm12Y","bpm16X","bpm16Y","bpm11X12X"};
  vector<TString> coil_array(7);
  for(int i=0;i<7;i++)
    coil_array[i] = Form("bmod_trim%d",i+1);
  
  vector<TString> det_name={"usl","usr"};
  vector<TString> bpm_name={"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm11X12X"};
  Int_t ndev = device_array.size();


  TFile *input = TFile::Open(Form("dit-coeffs/mc_ditcoeffs_%d.000.root",run));
  TTree *sens_tree = (TTree*)input->Get("sens");
  TTree *slope_tree = (TTree*)input->Get("dit_slope1");

  TMatrixD sens_matrix(ndev,7);
  for(int idev=0;idev<ndev;idev++){
    for(int icoil=0;icoil<7;icoil++){
      sens_tree->Draw(Form("%s_coil%d",device_array[idev].Data(),icoil+1),"","goff");
      TH1D *h1ptr = (TH1D*)gDirectory->FindObject("htemp");
      sens_matrix[idev][icoil] = h1ptr->GetMean();
    }
  }
  
  TMatrixD slope_matrix(2,5);
  for(int idet=0;idet<2;idet++){
    for(int ibpm=0;ibpm<5;ibpm++){
      slope_tree->Draw(Form("%s_%s",det_name[idet].Data(),bpm_name[ibpm].Data()),"","goff");
      TH1D *h1ptr = (TH1D*)gDirectory->FindObject("htemp");
      slope_matrix[idet][ibpm] = h1ptr->GetMean();
    }
  }

  TFile *sens_output = TFile::Open(Form("matrices/mc_sens_matrix.%d.root",run),"RECREATE");
  sens_output->cd();
  sens_output->WriteObject(&device_array,"dv_array");
  sens_output->WriteObject(&coil_array,"coil_array");
  sens_output->WriteObject(&sens_matrix,"sens_matrix");
  sens_output->Close();

  sens_output = TFile::Open(Form("matrices/mc_sens_matrix_wrong.%d.root",run),"RECREATE");
  /// >>> Introducing offset to dithering sensitivities
  sens_matrix[0][6] = sens_matrix[0][6]*1.2;
  sens_output->cd();
  sens_output->WriteObject(&device_array,"dv_array");
  sens_output->WriteObject(&coil_array,"coil_array");
  sens_output->WriteObject(&sens_matrix,"sens_matrix");
  sens_output->Close();

  TFile *slope_output = TFile::Open(Form("matrices/mc_dit_slope.%d.root",run),"RECREATE");
  slope_output->cd();
  slope_output->WriteObject(&det_name,"dv_array");
  slope_output->WriteObject(&bpm_name,"iv_array");
  slope_output->WriteObject(&slope_matrix,"slope_matrix");
  slope_output->Close();


  /// >>> Introducing offset to dithering slopes

  TMatrixD slope_matrix_false(2,5); //  [idet][imon]
  TMatrixD dD_dC(2,5); // [idet][icoil]
  TMatrixD dM_dC(5,5); // [icoil][imon]
  Int_t coil_idx[]={1,3,4,6,7};
  for(int idet=0;idet<2;idet++){
    for(int ic=0;ic<5;ic++){
      sens_tree->Draw(Form("%s_coil%d",det_name[idet].Data(),coil_idx[ic]),"","goff");
      TH1D *h1ptr = (TH1D*)gDirectory->FindObject("htemp");
      dD_dC[idet][ic] = h1ptr->GetMean();
    }
  }
  dD_dC[0][4] =   dD_dC[0][4]*1.2;

  for(int imon=0;imon<5;imon++){
    for(int ic=0;ic<5;ic++){
      sens_tree->Draw(Form("%s_coil%d",bpm_name[imon].Data(),coil_idx[ic]),"","goff");
      TH1D *h1ptr = (TH1D*)gDirectory->FindObject("htemp");
      dM_dC[imon][ic] = h1ptr->GetMean();
    }
  }

  TMatrixD inv_dM_dC(dM_dC);
  inv_dM_dC.Invert();
  slope_matrix_false = dD_dC * inv_dM_dC;
  
  TFile *false_slope_output = TFile::Open(Form("matrices/mc_dit_false_slope.%d.root",run),"RECREATE");
  false_slope_output->cd();
  false_slope_output->WriteObject(&det_name,"dv_array");
  false_slope_output->WriteObject(&bpm_name,"iv_array");
  false_slope_output->WriteObject(&slope_matrix_false,"slope_matrix");
  false_slope_output->Close();
  
}

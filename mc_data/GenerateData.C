#include  "parameters.hh"
typedef struct {Double_t hw_sum, block0,block1,block2, block3, Device_Error_Code;} LEAFLIST;

void GenerateData(Int_t run, Int_t random_seed=0){
  TDatime *fDatime = new TDatime();
  if(random_seed==0)
    random_seed =  fDatime->Get();
  TRandom3 *fRNG = new TRandom3(random_seed);

  const Int_t nBPM =  12;
  const Int_t nPar = 5;
  const Int_t nDet = 3;
  Int_t nMiniRun  = 20;
  Int_t mini_run_size = 9000;
  
  TFile* output  = TFile::Open(Form("mc_%d.000.root",run),"RECREATE");
  TTree *evt_tree = new TTree("evt"," evt tree (Monte Carlo Simulated)");
  TTree *mul_tree = new TTree("mul"," mul tree (Monte Carlo Simulated)");
  TString fDetName[nDet]={"usl","usr","us_avg"};
  TString fBPMName[nBPM]={"bpm4aX","bpm4aY","bpm4eX","bpm4eY",
			"bpm1X","bpm1Y","bpm11X","bpm11Y",
			"bpm12X","bpm12Y","bpm16X","bpm16Y"};
  TString fParName[nPar]={"pos_x","pos_y","angle_x","angle_y","energy"};

  // event level outputs;
  LEAFLIST fDet_evt[nDet];
  LEAFLIST fBPM_evt[nBPM];
  LEAFLIST fAlpha[nPar];
  LEAFLIST fDet_yield[nDet],fDet_asym[nDet];
  TString leaflist_str = "hw_sum/D:block0:block1:block2:block3:Device_Error_Code";
  for(int idet=0;idet<nDet-1;idet++){  // At event level ; just individual detectors
    evt_tree->Branch(fDetName[idet],&fDet_evt[idet],leaflist_str);
    fDet_evt[idet].Device_Error_Code = 0;
    
    mul_tree->Branch("yield_"+fDetName[idet],&fDet_yield[idet],leaflist_str);
    fDet_yield[idet].Device_Error_Code=0;
  }
  for(int ipar=0;ipar<nPar;ipar++){
    evt_tree->Branch(fParName[ipar],&fAlpha[ipar],leaflist_str);
    fAlpha[ipar].Device_Error_Code = 0;
  }
  for(int ibpm=0;ibpm<nBPM;ibpm++){
    evt_tree->Branch(fBPMName[ibpm],&fBPM_evt[ibpm],leaflist_str);
    fBPM_evt[ibpm].Device_Error_Code = 0;
  }

  Double_t fCodaEventNumber  = 0;
  Double_t fBmwObj, fBmodRamp;
  Double_t fBmwCycNum = 1000;
  LEAFLIST fBmodTrimCard[7];
  evt_tree->Branch("bmwobj",&fBmwObj);
  evt_tree->Branch("bmwcycnum",&fBmwCycNum);
  evt_tree->Branch("bmod_ramp",&fBmodRamp);
  evt_tree->Branch("CodaEventNumber",&fCodaEventNumber);
  for(int icoil=0;icoil<7;icoil++){
    evt_tree->Branch(Form("bmod_trim%d",icoil+1), &fBmodTrimCard[icoil],"value/D:v0:v1:v2:v3:Device_Error_Code");
    fBmodTrimCard[icoil].Device_Error_Code=0;
  }

  // pair / pattern level outputs
  LEAFLIST fDet_asym_ideal[nDet], fDet_yield_ideal[nDet];

  LEAFLIST fBPM_yield[nBPM];
  LEAFLIST fBPM_diff[nBPM];
  LEAFLIST fAlpha_yield[nPar];
  LEAFLIST fAlpha_diff[nPar];

  for(int idet=0;idet<nDet;idet++){
    mul_tree->Branch("asym_"+fDetName[idet],&fDet_asym[idet],leaflist_str);
    fDet_asym[idet].Device_Error_Code=0;
    mul_tree->Branch("asym_"+fDetName[idet]+"_ideal",&fDet_asym_ideal[idet],leaflist_str);
    fDet_asym_ideal[idet].Device_Error_Code=0;
  }

  for(int ipar=0;ipar<nPar;ipar++){
    mul_tree->Branch("diff_"+fParName[ipar],&fAlpha_diff[ipar],leaflist_str);
    mul_tree->Branch("yield_"+fParName[ipar],&fAlpha_yield[ipar],leaflist_str);
    fAlpha_yield[ipar].Device_Error_Code=0;
    fAlpha_diff[ipar].Device_Error_Code=0;
  }
  for(int ibpm=0;ibpm<nBPM;ibpm++){
    mul_tree->Branch("diff_"+fBPMName[ibpm],&fBPM_diff[ibpm],leaflist_str);
    mul_tree->Branch("yield_"+fBPMName[ibpm],&fBPM_yield[ibpm],leaflist_str);
    fBPM_yield[ibpm].Device_Error_Code=0;
    fBPM_diff[ibpm].Device_Error_Code=0;
  }
  
  Double_t fErrorFlag;
  Double_t fErrorFlag_pattern;
  Double_t fBmodFlag = 0x8000 + 0x1000 + 0x10000 + 0x4000000;
  Short_t fBurstCounter=0;
  mul_tree->Branch("ErrorFlag",&fErrorFlag_pattern);
  mul_tree->Branch("BurstCounter",&fBurstCounter);
  evt_tree->Branch("ErrorFlag",&fErrorFlag);
  evt_tree->Branch("BurstCounter",&fBurstCounter);
  Double_t um_to_mm = 1e-3;
  Double_t ppm_to_frac = 1e-6;
  Double_t fBPM_common_mode;
  Double_t fBmod_pulse;
  Int_t nGoodPattern=0;
  Int_t kHelicitySign = 1; // flippling between +1 or -1
  Int_t last_burst_counter=-1;
  Bool_t kInBmwCycle;
  Int_t fBmodCounts=0;
  Int_t fBmodSteps = 400;
  Int_t fBmodPeriod = 8;
  double det_flux_random[3];
  const double two_pi  = TMath::Pi()*2.0;
  while(fBurstCounter<nMiniRun){
    fCodaEventNumber ++;

    if(fBurstCounter>last_burst_counter){
      cout << "Mini-run " <<  fBurstCounter <<endl;
      last_burst_counter  = fBurstCounter;
      fBmwCycNum ++;
      kInBmwCycle = kTRUE;
      cout << "Supercycle " <<  fBmwCycNum  << " - bmwobj: ";
      fBmwObj = 1;
      cout << fBmwObj << "  ";
    }
    
    if(fBmodCounts==fBmodSteps){
      fBmwObj++;
      if(fBmwObj<=7)
	cout << fBmwObj << "  ";
      else{
	fBmwObj = -1;
	kInBmwCycle = kFALSE;
	cout << endl;
      }
      fBmodCounts = 0;
    }

    if(kInBmwCycle){
      fBmodCounts++;
      fBmodRamp  = fBmodCounts%fBmodPeriod;
      fBmod_pulse = 300.0*TMath::Sin(fBmodRamp/fBmodPeriod* two_pi);
      for(int icoil=0;icoil<7;icoil++){
	fBmodTrimCard[icoil].hw_sum = 0.0;
	fBmodTrimCard[icoil].Device_Error_Code = 0.0;
      }
      fBmodTrimCard[(int)fBmwObj-1].hw_sum = fBmod_pulse;
    } // end of if in bmod cycle
    
    det_flux_random[0] = fRNG->PoissonD(det_flux);
    fDet_evt[0].hw_sum = det_flux_random[0]; 
    det_flux_random[1] = fRNG->PoissonD(det_flux);
    fDet_evt[1].hw_sum = det_flux_random[1];
    
    for(int ipar=0;ipar<nPar;ipar++){
      fAlpha[ipar].hw_sum = fRNG->Gaus( 0, sigma_alpha[ipar]) + kHelicitySign * delta_alpha[ipar];
      // unit: um;urad;ppm
      if(kInBmwCycle)
      	fAlpha[ipar].hw_sum  += fBmod_pulse * dAlpha_dCoil[ipar][(int)fBmwObj-1];

      for(int idet=0;idet<nDet-1;idet++)
	fDet_evt[idet].hw_sum += det_flux_random[idet]* dDet_dAlpha[idet][ipar]* fAlpha[ipar].hw_sum * ppm_to_frac;

    }
    
    fBPM_common_mode = fRNG->Gaus(0,1);
    for(int ibpm=0;ibpm<nBPM;ibpm++){
      fBPM_evt[ibpm].hw_sum = fRNG->Gaus(0,sigma_bpm_idn[ibpm]*um_to_mm );
      for(int jbpm=0;jbpm<nBPM;jbpm++)
	fBPM_evt[ibpm].hw_sum += beta[ibpm][jbpm] * sigma_bpm_cor[jbpm] * fBPM_common_mode * um_to_mm;
      for(int ipar=0;ipar<nPar;ipar++)
	fBPM_evt[ibpm].hw_sum  += (dBPM_dAlpha[ibpm][ipar] * fAlpha[ipar].hw_sum)*um_to_mm;
    }
    
    if(kInBmwCycle){
      fErrorFlag = fBmodFlag;
    }else
      fErrorFlag = 0;
    
    evt_tree->Fill();
    
    // Processing Helicity Info
    if(kHelicitySign>0){
      fErrorFlag_pattern = fErrorFlag;
      for(int idet=0;idet<nDet-1;idet++){
	fDet_yield[idet].hw_sum = fDet_evt[idet].hw_sum;
	fDet_asym[idet].hw_sum = fDet_evt[idet].hw_sum;
	fDet_asym_ideal[idet].hw_sum = det_flux_random[idet];
	fDet_yield_ideal[idet].hw_sum = det_flux_random[idet];
      }
      for(int ipar=0;ipar<nPar;ipar++){
	fAlpha_yield[ipar].hw_sum =  fAlpha[ipar].hw_sum;
	fAlpha_diff[ipar].hw_sum =  fAlpha[ipar].hw_sum;
      }
      for(int ibpm=0;ibpm<nBPM;ibpm++){
	fBPM_yield[ibpm].hw_sum =  fBPM_evt[ibpm].hw_sum;
	fBPM_diff[ibpm].hw_sum =  fBPM_evt[ibpm].hw_sum;
      }
    }else{ // in the second window
      for(int ipar=0;ipar<nPar;ipar++){
	fAlpha_yield[ipar].hw_sum +=  fAlpha[ipar].hw_sum;
	fAlpha_diff[ipar].hw_sum -=  fAlpha[ipar].hw_sum;
	fAlpha_yield[ipar].hw_sum= fAlpha_yield[ipar].hw_sum/2.0;
	fAlpha_diff[ipar].hw_sum= fAlpha_diff[ipar].hw_sum/2.0;
      }
      for(int ibpm=0;ibpm<nBPM;ibpm++){
	fBPM_yield[ibpm].hw_sum +=  fBPM_evt[ibpm].hw_sum;
	fBPM_diff[ibpm].hw_sum -=  fBPM_evt[ibpm].hw_sum;
	fBPM_yield[ibpm].hw_sum= fBPM_yield[ibpm].hw_sum/2.0;
	fBPM_diff[ibpm].hw_sum= fBPM_diff[ibpm].hw_sum/2.0;
      }
      double fbcm_noise = fRNG->Gaus(0,bcm_res)*ppm_to_frac;
      for(int idet=0;idet<nDet-1;idet++){
	fDet_yield[idet].hw_sum += fDet_evt[idet].hw_sum;
	fDet_asym[idet].hw_sum -= fDet_evt[idet].hw_sum;
	fDet_yield[idet].hw_sum = fDet_yield[idet].hw_sum/2.0;
	fDet_asym[idet].hw_sum = fDet_asym[idet].hw_sum/2.0;
	fDet_asym[idet].hw_sum = fDet_asym[idet].hw_sum/fDet_yield[idet].hw_sum;
	fDet_asym[idet].hw_sum += fbcm_noise;
	
	fDet_yield_ideal[idet].hw_sum += det_flux_random[idet];
	fDet_asym_ideal[idet].hw_sum -= det_flux_random[idet];
	fDet_yield_ideal[idet].hw_sum = fDet_yield_ideal[idet].hw_sum/2.0;
	fDet_asym_ideal[idet].hw_sum = fDet_asym_ideal[idet].hw_sum/2.0;
	fDet_asym_ideal[idet].hw_sum = fDet_asym_ideal[idet].hw_sum/fDet_yield_ideal[idet].hw_sum;
      }

      fDet_asym_ideal[2].hw_sum = (fDet_asym_ideal[1].hw_sum + fDet_asym_ideal[0].hw_sum)/2.0;
      fDet_asym[2].hw_sum = (fDet_asym[1].hw_sum + fDet_asym[0].hw_sum)/2.0;
      
      mul_tree->Fill();
      fErrorFlag_pattern = ((int)fErrorFlag_pattern|(int)fErrorFlag);
      if(fErrorFlag_pattern==0){
	nGoodPattern ++;
	if(nGoodPattern%mini_run_size==0)
	  fBurstCounter++;
      }
    }
    // Flip Helicity Sign for next window
    if(kHelicitySign>0)
      kHelicitySign = -1;
    else
      kHelicitySign = 1;
    
  } // end of event loop
  mul_tree->Write();
  evt_tree->Write();
  
  TDirectory *par_dir = output->GetDirectory("/")->mkdir("parameters");
  par_dir->cd();
  TObjString* parameter_text = new TObjString();
  TObjString* random_seed_text = new TObjString();
  TString string_text="";
  
  std::ifstream file("parameters.hh");
  std::string line_str;
  std::string file_contents;
  file_contents += "parameters.hh \n";
  while(std::getline(file,line_str)){
    file_contents += line_str;
    file_contents.push_back('\n');
  }
  string_text = file_contents;
  parameter_text->SetString( string_text.Data());
  random_seed_text->SetString( Form("seed:%d",random_seed));
  par_dir->WriteTObject(parameter_text);
  par_dir->WriteTObject(random_seed_text);
  output->Close();
}

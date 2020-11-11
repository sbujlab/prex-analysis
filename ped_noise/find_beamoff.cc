#include "src/TaAccumulator.cc"
#include "src/TaEventRing.cc"
#include "LoadNormalizationMap.C"
#include "utilities.cc"
void find_beamoff(Int_t slug);
void find_beamoff(){
  for(int i=1;i<=94;i++)
    find_beamoff(i);
}
void find_beamoff(Int_t slug){
  ROOT::EnableImplicitMT(2);
  gStyle->SetStatH(0.35);
  gStyle->SetStatW(0.35);
  gStyle->SetOptStat("Mer");
  gStyle->SetStatFormat("6.2g");

  TStopwatch *tsw = new TStopwatch();
  
  map< Int_t,TString> fBCMRunMap = LoadNormalizationMap();
  TString runlist_filename = Form("prex-runlist/good_list/slug%d.list",slug);
  FILE *runlist = fopen(runlist_filename.Data(),"r");
  TString output_filename = Form("beamoff_cut/slug%d.txt",slug);
  FILE *output_pedcut = fopen(output_filename.Data(),"w");
  while(!feof(runlist)){
    Int_t run_number=0.0;
    fscanf(runlist,"%d/n",&run_number);
    if(run_number==0)
      continue;
    cout << "run " << run_number << endl;
    TString filename = Form("prexPrompt_pass1_%d.000.root",run_number);
    // TFile* input = TFile::Open("$QW_ROOTFILEs/"+filename);
    // TFile* input = TFile::Open("/media/yetao/prex/PREXII-respin1/"+filename);
    if(input==NULL)
      continue;
    if(fBCMRunMap.find(run_number)==fBCMRunMap.end()){
      cerr << "-- normalizing BCM info not found for run  "
	   << run_number << endl;
      continue;  // FIXME
    }
  
    TTreeReader myReader("mul",input);
    TString bcm_name = "yield_"+fBCMRunMap[run_number]+".hw_sum";
    TTreeReaderValue<Double_t> fBCMValue(myReader,bcm_name);
    TTreeReaderValue<Double_t> fErrorFlag(myReader,"ErrorFlag");
    TTreeReaderValue<Double_t> fPatternCounter(myReader,"pattern_number");
    Int_t pat_counter=0;
    TaEventRing* fEventRing = new TaEventRing();
    // ================ Determine beam off readout in uA
    TH1D h_bcm("h_bcm","h_bcm",21000,-10.5,199.5);
    TTree *mul_tree = (TTree*)input->Get("mul");
    mul_tree->Draw(bcm_name+">>h_bcm","","goff");
    Int_t bin_counter = h_bcm.FindFirstBinAbove(0.0);
    Double_t bin_center = h_bcm.GetBinCenter(bin_counter) ;
    if(bin_center>2.5){
      cerr << "-- Error: lowest current readback is "
	   << bin_center << " higher than 2.5 uA  "<< endl;
      continue;
    }
    Double_t beam_off_uA = bin_center;
    cout << "beam off uA " << beam_off_uA << endl;
    Double_t bumper = 0.03;
    fEventRing->SetBeamOffLimit(beam_off_uA+bumper);
    cout << "Set beam-off Limit(uA) " << beam_off_uA+bumper << endl;
  
    // === Start of Searching 
    Double_t fPopback;
    vector<Double_t> fEventCounterArray;
    while(myReader.Next()){
      if(pat_counter>10){
	fEventRing->PushBeamCurrent(*fBCMValue);
	fEventRing->PushEventCounter(*fPatternCounter);
	if( fEventRing->isReady() ){
	  Int_t fPopErrorFlag = fEventRing->PopEventCounter(fPopback);
	  if(fPopErrorFlag==0)
	    fEventCounterArray.push_back(fPopback);
	}
      }
      pat_counter++;
    }
    TString beamoff_cut = generate_cut(fEventCounterArray);
    cout << beamoff_cut << endl;
    // ==== Done with Searching
    if(beamoff_cut=="1==0"){
      input->Close();
      continue;
    }
    fprintf(output_pedcut,"%d:(%s)\n",run_number,beamoff_cut.Data());
    
    vector<TString> device_list={"yield_bcm_an_us","yield_bcm_an_ds","yield_bcm_dg_us","yield_bcm_dg_ds",
				 "diff_bcm_an_us","diff_bcm_an_ds","diff_bcm_dg_us","diff_bcm_dg_ds",
				 "diff_usl","diff_dsl","diff_atl1","diff_atl2",
				 "diff_usr","diff_dsr","diff_atr1","diff_atr2",
				 "diff_bpm4aWS","diff_bpm4eWS",
				 "diff_bpm4acWS","diff_bpm4ecWS",
				 "diff_bpm11WS","diff_bpm12WS","diff_bpm16WS",
				 "diff_bpm1WS","diff_cav4cQ",
				 "diff_bpm1p02bWS","diff_bpm1p03aWS",
				 "yield_usl","yield_dsl","yield_atl1","yield_atl2",
				 "idot_check","something_I_put_here","for_fun",
				 "yield_usr","yield_dsr","yield_atr1","yield_atr2",
				 "yield_bpm4aWS","yield_bpm4eWS",
				 "yield_bpm4acWS","yield_bpm4ecWS",
				 "yield_bpm11WS","yield_bpm12WS","yield_bpm16WS",
				 "yield_bpm1WS","yield_cav4cQ",
				 "yield_bpm1p02bWS","yield_bpm1p03aWS"};
  
    TTree *mulc_tree = (TTree*)input->Get("mulc");
    mul_tree->AddFriend(mulc_tree);
  
    TCanvas *c1 = new TCanvas("c1","c1",1600,400);
    c1->Divide(4,1);

    auto iter_dev = device_list.begin();
    while(iter_dev!=device_list.end()){
      TString ch_name = *iter_dev;
      if(mul_tree->GetBranch(ch_name)==NULL){
	if(ch_name.Contains("yield")){
	  device_list.erase(iter_dev);
	}else if(ch_name.Contains("diff")){
	  TString alias_name = ch_name;
	  TString raw_name = ch_name.ReplaceAll("diff_","");
	  if(mul_tree->GetBranch("yield_"+raw_name)==NULL){
	    device_list.erase(iter_dev);
	    continue;
	  }
	  mul_tree->SetAlias(alias_name,Form("asym_%s * yield_%s",raw_name.Data(),raw_name.Data()));
	  iter_dev++;
	}else{
	  device_list.erase(iter_dev);
	}
      }else
	iter_dev++;
    }
  
    c1->Print(Form("./plots/run%d_beam_off.pdf[",run_number));
    Int_t plot_counter=0.0;
    Int_t ndev = device_list.size();
    for(int idev=0;idev<ndev;idev++){
      c1->cd(idev%4+1);
      TString unit="";
      Double_t rescale=1.0;
      if(device_list[idev].Contains("yield_bcm")){
	rescale = 1;
	unit = " (uA) ";
      }
      if(device_list[idev].Contains("diff_bcm")){
	rescale = 1e6;
	unit = " (pA) ";
      }
      if(device_list[idev].Contains("diff_") &&
	 (device_list[idev].Contains("l") || device_list[idev].Contains("r"))){
	rescale = 1e6;
	unit = " (uV) ";
      }
      if(device_list[idev].Contains("diff_bpm")
	 && device_list[idev].Contains("WS")){
	rescale = 76e-6*1e6;
	unit = " (uV) ";
      }
      TString asym_chname = device_list[idev];
      asym_chname.ReplaceAll("diff_","asym_");
      asym_chname.ReplaceAll("yield_","asym_");
      TString not_blinder = Form("&&((%s.Device_Error_Code&512)!=512)",asym_chname.Data());
      mul_tree->Draw(Form("%s*%f",device_list[idev].Data(),rescale),
		     "("+beamoff_cut+")"+not_blinder);
      
      TH1D *hptr = (TH1D*)gPad->FindObject("htemp");
      if(hptr!=NULL){
	hptr->SetTitle(device_list[idev]+unit);
	hptr->SetName(Form("h%d",plot_counter++));
      }
      if( idev%4==3 || idev==ndev){
	c1->Print(Form("./plots/run%d_beam_off.pdf",run_number));
	c1->Clear("D");
      }
    }
    c1->Print(Form("./plots/run%d_beam_off.pdf]",run_number));
    
    input->Close();
  } // end of runlist loop
  
  fclose(output_pedcut);
  tsw->Print();
}



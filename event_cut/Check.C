#include "Plot.C"
void Check(){
  TFile *file1 = TFile::Open("BPMCheck_shit.root");
  TTree* T1 = (TTree*)file1->Get("T");
 
  vector<pair<TString,TString> > fConditions=
    {
     {"3*yield_usl.fwhm/2.335/yield_usl>0.05 && yield_usl>0.02","usl 3-sigma width greater than 5\%" },
     {"(yield_usl.ul-yield_usl)/yield_usl>0.05 && yield_usl>0.02","usl upper deviation greater than 5\%"},
     {"(yield_usl-yield_usl.ll)/yield_usl>0.05 && yield_usl>0.02","usl lower deviation greater than 5\%"},
     // {"(yield_usl.ul-yield_usl.ll)<5.0/2.335*yield_usl.fwhm && yield_usl>0.02","usl likely double peaking"},
     {"3*yield_usr.fwhm/2.335/yield_usr>0.05 && yield_usr>0.02","usr 3-sigma width greater than 5\%" },
     {"(yield_usr.ul-yield_usr)/yield_usr>0.05 && yield_usr>0.02","usr upper deviation greater than 5\%"},
     {"(yield_usr-yield_usr.ll)/yield_usr>0.05 && yield_usr>0.02","usr lower deviation greater than 5\%"},
     // {"(yield_usr.ul-yield_usr.ll)<5.0/2.335*yield_usr.fwhm && yield_usr>0.02","usr likely double peaking"},
     // {"3*yield_bpm4aX.fwhm/2.335/um>100 && yield_bpm4aX.num_samples>0","bpm4aX 3-sigma width greater than 100 um"},
     // {"(yield_bpm4aX.ul-yield_bpm4aX)/um>100 && yield_bpm4aX.num_samples>0", "bpm4aX upper deviation deviated more than 100 um"},
     // {"(yield_bpm4aX-yield_bpm4aX.ll)/um>100 && yield_bpm4aX.num_samples>0", "bpm4aX lower deviation deviated more than 100 um"},
     // {"(yield_bpm4aX.ul-yield_bpm4aX.ll)<3*2.335*2*yield_bpm4aX.fwhm&&yield_bpm4aX.num_samples>0","bpm4aX likely double peaking "},
     // {"(yield_bpm4aX.ul-yield_bpm4aX.ll)>5*2.335*2*yield_bpm4aX.fwhm&&yield_bpm4aX.num_samples>0","bpm4aX likely has outliers "},
     // {"3*yield_bpm4eX.fwhm/2.335/um>300 && yield_bpm4eX.num_samples>0", "bpm4aX 3-sigma witdh greater than 300 um" },
     // {"(yield_bpm4eX.ul-yield_bpm4eX)/um>100 && yield_bpm4eX.num_samples>0", "bpm4eX upper deviation deviated more than 100 um"},
     // {"(yield_bpm4eX-yield_bpm4eX.ll)/um>100 && yield_bpm4eX.num_samples>0", "bpm4eX lower deviation deviated more than 100 um"},
     // {"(yield_bpm4eX.ul-yield_bpm4eX.ll)<3.0*2.335*2*yield_bpm4eX.fwhm&&yield_bpm4eX.num_samples>0","bpm4eX most likely to be double peaking "},
     // {"3*yield_bpm4aY.fwhm/2.335/um>300 &&yield_bpm4aY.num_samples>0","bpm4aY 3-sigma width greater than 300 um"},
     {"(yield_bpm4aY.ul-yield_bpm4aY)/um>2000&&yield_bpm4aY.num_samples>0", "bpm4aY upper deviation deviated more than 2000 um"},
     {"(yield_bpm4aY-yield_bpm4aY.ll)/um>2000&&yield_bpm4aY.num_samples>0", "bpm4aY lower deviation deviated more than 2000 um"},
     // {"(yield_bpm4aY.ul-yield_bpm4aY.ll)<3*2.335*2*yield_bpm4aY.fwhm&&yield_bpm4aY.num_samples>0","bpm4aY most likely to be double peaking "},
     // {"3*yield_bpm4eY.fwhm/2.335/um>200 &&yield_bpm4eY.num_samples>0","bpm4eY 3-sigma width greater than 200 um"},
     // {"(yield_bpm4eY.ul-yield_bpm4eY)/um>100&&yield_bpm4eY.num_samples>0", "bpm4eY upper deviation deviated more than 100 um"},
     // {"(yield_bpm4eY-yield_bpm4eY.ll)/um>100&&yield_bpm4eY.num_samples>0", "bpm4eY lower deviation deviated more than 100 um"},
     // {"(yield_bpm4eY.ul-yield_bpm4eY.ll)<3*2.335*2*yield_bpm4eY.fwhm&&yield_bpm4eY.num_samples>0","bpm4eY most likely to be double peaking "},
     // {"3*yield_bpm12X.fwhm/2.335/um> 200 &&yield_bpm12X.num_samples>0","bpm12X 3-sigma width greater than 200 um"},
     // {"(yield_bpm12X.ul-yield_bpm12X)/um>100&&yield_bpm12X.num_samples>0", "bpm12X upper deviation deviated more than 100 um"},
     // {"(yield_bpm12X-yield_bpm12X.ll)/um>100&&yield_bpm12X.num_samples>0", "bpm12X lower deviation deviated more than 100 um"},
     // {"(yield_bpm12X.ul-yield_bpm12X.ll)<3.0*2.335*2*yield_bpm12X.fwhm&&yield_bpm12X.num_samples>0","bpm12X most likely to be double peaking "}
    };

  map<Int_t, vector<TString> > fDescriptionMap; // description of problem
  TString cut_text;
  auto iter = fConditions.begin();
  while(iter!=fConditions.end()){
    cut_text += Form("&& (!(%s))",(*iter).first.Data());
    Int_t npt = T1->Draw("run", (*iter).first,"goff");
    Double_t *fRun = T1->GetV1();
    if(npt==0){
      cout << " no run " << (*iter).second << endl;
    }
    for(int i=0;i<npt;i++)
      fDescriptionMap[ fRun[i] ].push_back( (*iter).second );
    iter++;
  } // loop over conditions


  
  cout << " == Summary == " << endl;
  auto irun = fDescriptionMap.begin();
  while(irun!=fDescriptionMap.end()){
    cout << " ** run " << (*irun).first << endl;
    Int_t ndesc = (*irun).second.size();
    for(int i=0;i<ndesc;i++)
      cout << "\t " << ((*irun).second)[i] << endl;
    cout << "\n" << endl;
    irun++;
  }
  
  Plot(cut_text);
}  


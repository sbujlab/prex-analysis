void CycleCrawler(Int_t run_number){
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.1f");
  TString japanOutput_path=getenv("QW_ROOTFILES");
  TString japanOutput_prefix="/prexPrompt_pass1_";
  TString filename=japanOutput_path+japanOutput_prefix+Form("%d.000.root",run_number);
  TFile *japanOutput = TFile::Open(filename);
  if(japanOutput==NULL)
    return;
  TTree *evt_tree = (TTree*)japanOutput->Get("evt");
  if(evt_tree==NULL)
    return;
  FILE *report = fopen("./crawler.txt","a+");
  Int_t nCoils=7;
  Double_t bmwcycnum;
  evt_tree->SetBranchAddress("bmwcycnum",&bmwcycnum);
  TEventList *elist = new TEventList("elist");
  evt_tree->Draw(">>+elist","bmwcycnum>0 && bmod_ramp>0 && bmwobj>0");
  evt_tree->SetEventList(elist);
  TH1D *hcyc = new TH1D("hcyc","",5000,-0.5,4999.5);
  evt_tree->Draw("bmwcycnum>>hcyc","","goff");
  Int_t bin_first = hcyc->FindFirstBinAbove(0);
  Int_t bin_last = hcyc->FindLastBinAbove(0);
  Int_t cycID_first = hcyc->GetBinCenter(bin_first);
  Int_t cycID_last = hcyc->GetBinCenter(bin_last);
  Int_t nCycle = cycID_last-cycID_first+1;
  TString *cycle_label = new TString[7*nCycle];
  for(int icyc=0;icyc<nCycle;icyc++)
    for(int icoil=0;icoil<nCoils;icoil++)
      cycle_label[icoil+icyc*nCoils]=Form("%d:%d",cycID_first+icyc,icoil+1);

  TCanvas *c1 = new TCanvas("c1","c1",nCycle*500,500);
  c1->cd();
  c1->SetLogz();
  Bool_t kFound =kFALSE;
  TH2D *h2d = new TH2D("h2d","",7*nCycle,-0.5,7*nCycle-0.5,7,0.5,7+0.5);
  for(int icyc=0;icyc<nCycle;icyc++){
    for(int icoil=0;icoil<nCoils;icoil++){
      for(int iobj=0;iobj<nCoils;iobj++){
        evt_tree->Draw(Form("bmod_trim%d",icoil+1),
                       Form("bmwobj==%d && bmwcycnum==%d",iobj+1,cycID_first+icyc),
                       "goff");
        TH1D* htemp=(TH1D*)gDirectory->FindObject("htemp");
        double rms=htemp->GetRMS();
        if(rms>10 && icoil!=iobj){
	  if(!kFound)                                                                                           
            fprintf(report," == run %d ==\n",run_number);                                                       
          kFound=kTRUE;                                                                                         
          fprintf(report," -- cycleID %d : bmod_trimd%d at bmwobj %d \n" ,                                      
               cycID_first+icyc,icoil+1,iobj+1);                                                                
        }
        h2d->Fill(cycle_label[iobj+nCoils*icyc],icoil+1,rms);
      }
    }
  }
  h2d->SetMarkerSize(1.5);
  h2d->SetTitle(Form("run %d : bmod_trim vs bmwobj ;bmwcycnum:bmwobj ; trimcard #",run_number));
  h2d->GetYaxis()->SetNdivisions(7);
  h2d->GetXaxis()->SetNdivisions(700+nCycle);
  h2d->GetXaxis()->SetLabelSize(0.04);
  h2d->GetXaxis()->SetTitleSize(0.04);
  h2d->Draw("COL TEXT");
  if(kFound)
    c1->SaveAs(Form("./crawler_output/run%d_cycle.png",run_number));

  fclose(report);
}


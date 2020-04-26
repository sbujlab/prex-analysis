void BCMStatus(){
  TFile *input = TFile::Open("prex_bcm_status.root");
  TTree *T = (TTree*) input->Get("T");
  gStyle->SetOptStat(0);
  // gStyle->SetNumberContours(3);
  TH2D *h2d = new TH2D("h2d","",1700,3300-0.5,5000-0.5,4,-0.5,3.5);
  TH2D *hslug = new TH2D("hslug","",95,-0.5,95-0.5,4,-0.5,3.5);
  h2d->GetXaxis()->SetTitle("Run Number");
  hslug->GetXaxis()->SetTitle("Slug Number");
  Int_t nEntries = T->GetEntries();
  TString device_list[]={"bcm_an_us","bcm_an_ds",
			 "bcm_an_ds3","bcm_an_ds10",
			 "bcm_dg_us","bcm_dg_ds",
			 "unser"};
  Int_t slug,run_number,bcm_index;
  Int_t bcm_flag[7];
  T->SetBranchAddress("run",&run_number);
  T->SetBranchAddress("slug",&slug);
  T->SetBranchAddress("bcm_index",&bcm_index);
  for(int i=0;i<7;i++)
    T->SetBranchAddress(device_list[i]+"_flag",&bcm_flag[i]);

  for(int ievt=0;ievt<nEntries;ievt++){
    T->GetEntry(ievt);
    for(int i=0;i<7;i++){
      if(bcm_flag[i]!=0){
	h2d->Fill(run_number,device_list[i].Data(),bcm_flag[i]);
	hslug->Fill(slug,device_list[i].Data(),bcm_flag[i]);
      }
    }
  }
  TCanvas *c1 = new TCanvas("c1","c1",800,300);
  c1->cd();
  // h2d->GetZaxis()->SetRangeUser(-0.25,1.25);
  h2d->GetYaxis()->SetLabelSize(0.05);
  h2d->Draw("BOX");
  c1->SaveAs("normalizing_bcm_history_run.pdf");
  TCanvas *c2 = new TCanvas("c2","c2",800,300);
  c2->cd();
  hslug->SetContour(30);
  hslug->GetYaxis()->SetLabelSize(0.05);
  hslug->GetZaxis()->SetRangeUser(0,30);
  hslug->Draw("COLZ0");
  c2->SaveAs("normalizing_bcm_history_slug.pdf");
}

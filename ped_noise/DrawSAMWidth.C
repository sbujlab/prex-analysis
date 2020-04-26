void DrawSAMWidth(){
  gStyle->SetOptStat(0);
  TFile *red_file = TFile::Open("postpan.root");
  TFile *ped_file = TFile::Open("test.root");
  TTree *ped_tree = (TTree*)ped_file->Get("ped");
  TTree *mini_tree = (TTree*)red_file->Get("T");
  TString user_cut = "run!=5616 && run!=5619 && run!=5620 && run!=5649";
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->cd();
  c1->Print("sam_pedestal_noise.pdf[","pdf");
  for(int isam=1;isam<=8;isam++){
    mini_tree->Draw(Form("reg_asym_sam%d.rms/ppm:run",isam),
		    user_cut+Form("&&reg_asym_sam%d.rms>0",isam),"prof");
    TH1D *hreg = (TH1D*)gPad->FindObject("htemp");
    hreg->SetName("reg");
    hreg->SetTitle(Form("SAM %d",isam));
    hreg->SetMarkerColor(kMagenta);
    hreg->SetMarkerStyle(20);
    hreg->GetYaxis()->SetTitle("Width (ppm)");
    hreg->GetYaxis()->SetTitleSize(0.03);
    hreg->GetYaxis()->SetLabelSize(0.03);
    hreg->GetXaxis()->SetTitle("run number");
    hreg->GetXaxis()->SetTitleSize(0.03);
    hreg->GetXaxis()->SetLabelSize(0.03);
    double ymax = hreg->GetMaximum();
    double ymin =0;
    hreg->GetYaxis()->SetRangeUser(ymin,ymax+0.5*(ymax-ymin));
    
    ped_tree->Draw(Form("sam%d_noise:run",isam),"","same");
    TH1D *hped= new TH1D("hped","hped",100,0,10);
    hped->SetMarkerStyle(20);
    TLegend leg(0.9,0.9,0.7,0.7);
    leg.AddEntry(hreg,"regression width","p");
    leg.AddEntry(hped,"pedestal noise","p");
    leg.Draw("same");
    c1->Print("sam_pedestal_noise.pdf");
  }
  c1->Print("sam_pedestal_noise.pdf]");
  // g_ped->SetName("ped");
  // TMultiGraph *mg = new TMultiGraph();
  // mg->Add(g_reg,"p");
  
  // mg->Add(g_ped,"p");
  // mg->Draw("A");
  // g_buff = (TGraph*)gPad->FindObject("Graph");
}

void DrawSAMWidth(){

  TFile *red_file = TFile::Open("postpan.root");
  TFile *ped_file = TFile::Open("./ped_rootfiles/merged.root");
  TTree *ped_tree = (TTree*)ped_file->Get("ped");
  TTree *mini_tree = (TTree*)red_file->Get("T");
  TString user_cut="run!=5616 && run!=5619 && run!=5620";
  int isam=4;

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  mini_tree->Draw(Form("reg_asym_sam%d.rms/ppm:run",isam),
		  user_cut,"prof");
  TGraph* g_reg = (TGraph*)gPad->FindObject("Graph");
  g_reg->SetName("regress");
  ped_tree->SetMarkerStyle(20);
  ped_tree->SetMarkerColor(kMagenta);
  ped_tree->Draw(Form("sam%d_noise:run",isam),"","prof");
  TGraph* g_ped = (TGraph*)gPad->FindObject("Graph");
  g_ped->SetName("ped");
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(g_reg,"p");
  mg->Add(g_ped,"p");
  mg->Draw("A");
  // g_buff = (TGraph*)gPad->FindObject("Graph");
}

void bmwTestRawRMS(){

  TFile* reg_input = TFile::Open("regall_prod_only_grand_average_signed_weighted.root");
  TFile* reg_bmw_input = TFile::Open("regall_bmw_only_grand_average_signed_weighted.root");
  
  TTree *reg_tree = (TTree*)reg_input->Get("slug");
  TTree *reg_bmw_tree = (TTree*)reg_bmw_input->Get("slug");
  reg_tree->AddFriend(reg_bmw_tree,"reg_bmw");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->cd();

  reg_tree->Draw("reg_bmw.AdetRaw.rms*1e6*1.05 : slug","arm_flag!=0","goff");
  TGraph *mg_single = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
  mg_single->SetMarkerStyle(23);
  
  reg_tree->Draw("AdetRaw.rms*1e6 : slug","","goff");
  TGraph *mg_prod = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
  mg_prod->SetLineColor(kBlue);
  mg_prod->SetLineStyle(7);
  mg_prod->SetLineWidth(2);
  
  reg_tree->Draw("reg_bmw.AdetRaw.rms*1e6 : slug","reg_bmw.AdetRaw.error>0","goff");
  TGraph *mg_bmw = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
  mg_bmw->SetLineColor(kBlue);
  mg_bmw->SetLineWidth(2);

  TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
  leg->AddEntry(mg_prod, " Raw: Production Only");
  leg->AddEntry(mg_bmw, "  Raw: BMW Only");
  leg->AddEntry(mg_single, " Single Arm","p");
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(mg_prod,"l");
  mg->Add(mg_bmw,"l");
  mg->Add(mg_single,"p");

  mg->Draw("A");
  mg->GetXaxis()->SetTitle("Slug Number");
  mg->GetYaxis()->SetTitle(" RMS (ppm)");
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->SetTitle();
  leg->Draw("same");
  c1->SaveAs("raw_rms_by_slug.pdf");
}

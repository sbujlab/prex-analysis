void bmwOnlyRMSBySlug(){


  TFile* reg_input = TFile::Open("regall_prod_only_grand_average_signed_weighted.root");
  TFile* reg_bmw_input = TFile::Open("regall_bmw_only_grand_average_signed_weighted.root");

  TFile* lagr_input = TFile::Open("lagrall_prod_only_grand_average_signed_weighted.root");
  TFile* lagr_bmw_input = TFile::Open("lagrall_bmw_only_grand_average_signed_weighted.root");

  TFile* dit_input = TFile::Open("dit_prod_only_grand_average_signed_weighted.root");
  TFile* dit_bmw_input = TFile::Open("dit_bmw_only_grand_average_signed_weighted.root");

  TTree *lagr_tree = (TTree*)lagr_input->Get("slug");
  TTree *reg_tree = (TTree*)reg_input->Get("slug");
  TTree *dit_tree = (TTree*)dit_input->Get("slug");
  
  TTree *lagr_bmw_tree = (TTree*)lagr_bmw_input->Get("slug");
  TTree *reg_bmw_tree = (TTree*)reg_bmw_input->Get("slug");
  TTree *dit_bmw_tree = (TTree*)dit_bmw_input->Get("slug");
  
  lagr_tree->AddFriend(reg_tree,"reg");
  lagr_tree->AddFriend(dit_tree,"dit");
  lagr_tree->AddFriend(reg_bmw_tree,"reg_bmw");
  lagr_tree->AddFriend(lagr_bmw_tree,"lagr_bmw");
  lagr_tree->AddFriend(dit_bmw_tree,"dit_bmw");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->cd();
  // c1->SetGridy();
  // c1->SetGridx();

  lagr_tree->Draw("lagr_bmw.Adet.rms*1e6*1.05 : slug","arm_flag!=0","goff");
  TGraph *mg_single = new TGraph(lagr_tree->GetSelectedRows(),lagr_tree->GetV2(),lagr_tree->GetV1());
  mg_single->SetMarkerStyle(23);
  
  lagr_tree->Draw("Adet.rms*1e6 : slug","","goff");
  TGraph *mg_lagr = new TGraph(lagr_tree->GetSelectedRows(),lagr_tree->GetV2(),lagr_tree->GetV1());
  mg_lagr->SetLineColor(kBlue);
  mg_lagr->SetLineStyle(7);
  mg_lagr->SetLineWidth(1);
  
  lagr_tree->Draw("reg.Adet.rms*1e6 : slug","","goff");
  TGraph *mg_reg = new TGraph(lagr_tree->GetSelectedRows(),lagr_tree->GetV2(),lagr_tree->GetV1());
  mg_reg->SetLineColor(kRed);
  mg_reg->SetLineStyle(7);
  mg_reg->SetLineWidth(1);
    
  lagr_tree->Draw("lagr_bmw.Adet.rms*1e6 : slug","lagr_bmw.Adet.error>0","goff");
  TGraph *mg_lagr_bmw = new TGraph(lagr_tree->GetSelectedRows(),lagr_tree->GetV2(),lagr_tree->GetV1());
  mg_lagr_bmw->SetMarkerColor(kBlue);
  mg_lagr_bmw->SetLineColor(kBlue);
  mg_lagr_bmw->SetMarkerStyle(7);
    
  lagr_tree->Draw("reg_bmw.Adet.rms*1e6 : slug","reg_bmw.Adet.error>0","goff");
  TGraph *mg_reg_bmw = new TGraph(lagr_tree->GetSelectedRows(),lagr_tree->GetV2(),lagr_tree->GetV1());
  mg_reg_bmw->SetMarkerColor(kRed);
  mg_reg_bmw->SetLineColor(kRed);
  mg_reg_bmw->SetMarkerStyle(7);

  lagr_tree->Draw("dit_bmw.Adet.rms*1e6 : slug","dit_bmw.Adet.error>0","goff");
  TGraph *mg_dit_bmw = new TGraph(lagr_tree->GetSelectedRows(),lagr_tree->GetV2(),lagr_tree->GetV1());
  mg_dit_bmw->SetMarkerColor(kOrange);
  mg_dit_bmw->SetLineColor(kOrange);
  mg_dit_bmw->SetMarkerStyle(7);


  TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
  leg->AddEntry(mg_lagr, " Lagrange : Production Only");
  leg->AddEntry(mg_reg, " Regression : Production Only");
  leg->AddEntry(mg_lagr_bmw, " Lagrange : BMW Only");
  leg->AddEntry(mg_reg_bmw, " Regression : BMW Only");
  leg->AddEntry(mg_dit_bmw, " Dithering : BMW Only");
  leg->AddEntry(mg_single, " Single Arm","p");
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(mg_lagr,"l");
  mg->Add(mg_reg,"l");
  mg->Add(mg_lagr_bmw,"lp");
  mg->Add(mg_reg_bmw,"lp");
  mg->Add(mg_dit_bmw,"lp");
  mg->Add(mg_single,"p");

  mg->Draw("A");
  mg->GetXaxis()->SetTitle("Slug Number");
  mg->GetYaxis()->SetTitle(" RMS (ppm)");
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->SetTitle();
  leg->Draw("same");
  c1->SaveAs("rms_by_slug.pdf");
}

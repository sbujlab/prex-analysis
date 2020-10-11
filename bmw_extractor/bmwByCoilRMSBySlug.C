void bmwByCoilRMSBySlug(TString det_name="us_dd",
			TString user_cut="arm_flag==0"){
  
  TFile* reg_input = TFile::Open("regall_prod_only_grand_average_signed_weighted.root");
  TTree *reg_tree = (TTree*)reg_input->Get("slug");
  TFile* lagr_input = TFile::Open("lagrall_prod_only_grand_average_signed_weighted.root");
  TTree *lagr_tree = (TTree*)lagr_input->Get("slug");
  
  int coil_array[]={7,1,3,5,2,4,6};
  int ncoil = sizeof(coil_array)/sizeof(*coil_array);
  
  for(int icoil=0;icoil<ncoil;icoil++){
    TFile* reg_bmw_input = TFile::Open(Form("regall_bmw_coil%d_grand_average_signed_weighted.root",coil_array[icoil]));
    TTree *reg_bmw_tree = (TTree*)reg_bmw_input->Get("slug");
    reg_tree->AddFriend(reg_bmw_tree,Form("c%d",coil_array[icoil]));
    TFile* lagr_bmw_input = TFile::Open(Form("lagrall_bmw_coil%d_grand_average_signed_weighted.root",coil_array[icoil]));
    TTree *lagr_bmw_tree = (TTree*)lagr_bmw_input->Get("slug");
    lagr_tree->AddFriend(lagr_bmw_tree,Form("c%d",coil_array[icoil]));
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(1,2);
  c1->Print(Form("%s_coil_rms_by_slug.pdf[",det_name.Data()));
  for(int ic=0;ic<ncoil;ic++){
    c1->Clear("D");
    c1->cd(1);
    // reg_tree->Draw("reg_bmw.Adet.rms*1e6*1.05 : slug","arm_flag!=0","goff");
    // TGraph *mg_single = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
    // mg_single->SetMarkerStyle(23);
    TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
    TMultiGraph *mg = new TMultiGraph();
    Int_t color_pa[]={1,2,3,4,6,7,8};
  
    reg_tree->Draw("reg_asym_"+det_name+".rms*1e6 : slug",user_cut,"goff");
    TGraph *mg_reg = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
    // mg_reg->SetLineColor(kMagenta);
    // mg_reg->SetMarkerColor(kMagenta);
    mg_reg->SetLineStyle(7);
    mg_reg->SetLineWidth(1);
    mg_reg->SetMarkerStyle(7);
    mg->Add(mg_reg,"l");
    leg->AddEntry(mg_reg, " Regression : Production Only");
  

    reg_tree->Draw(Form("c%d.reg_asym_%s.rms*1e6 : slug",coil_array[ic],det_name.Data()),
		   user_cut+Form("&&c%d.reg_asym_%s.error>0",coil_array[ic],det_name.Data()),"goff");
    TGraph *mg_reg_bmw = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
    mg_reg_bmw->SetMarkerColor(kRed);
    mg_reg_bmw->SetLineColor(kRed);
    mg_reg_bmw->SetMarkerStyle(7);
    mg->Add(mg_reg_bmw,"lp");
    leg->AddEntry(mg_reg_bmw,Form("Regression: coil%d",coil_array[ic]));


    lagr_tree->Draw(Form("c%d.lagr_asym_%s.rms*1e6 : slug",coil_array[ic],det_name.Data()),
		   user_cut+Form("&& c%d.lagr_asym_%s.error>0",coil_array[ic],det_name.Data()),"goff");
    TGraph *mg_lagr_bmw = new TGraph(lagr_tree->GetSelectedRows(),lagr_tree->GetV2(),lagr_tree->GetV1());
    mg_lagr_bmw->SetMarkerColor(kBlue);
    mg_lagr_bmw->SetLineColor(kBlue);
    mg_lagr_bmw->SetMarkerStyle(7);
    mg->Add(mg_lagr_bmw,"lp");
    leg->AddEntry(mg_lagr_bmw,Form("Lagrange: coil%d",coil_array[ic]));

    mg->SetTitle(Form("%s Corrected Width(ppm): BMW Coil%d Test",det_name.Data(), coil_array[ic]));
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("Slug Number");
    mg->GetYaxis()->SetTitle(" RMS (ppm)");
    mg->GetYaxis()->SetTitleSize(0.05);
    leg->Draw("same");

    c1->cd(2);
    TLegend *legRaw = new TLegend(0.7,0.6,0.9,0.9);
    TMultiGraph *mgRaw = new TMultiGraph();
  
    reg_tree->Draw("asym_"+det_name+".rms*1e6 : slug",user_cut,"goff");
    TGraph *gRaw = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
    gRaw->SetLineColor(kBlack);
    gRaw->SetLineStyle(7);
    gRaw->SetLineWidth(1);
    gRaw->SetMarkerStyle(7);
    mgRaw->Add(gRaw,"lp");
    legRaw->AddEntry(gRaw, " Raw : Production Only");

    reg_tree->Draw(Form("c%d.asym_%s.rms*1e6 : slug",coil_array[ic],det_name.Data()),
		   user_cut+Form("&&c%d.asym_%s.error>0",coil_array[ic],det_name.Data()),"goff");
    TGraph *gRaw_bmw = new TGraph(reg_tree->GetSelectedRows(),reg_tree->GetV2(),reg_tree->GetV1());
    gRaw_bmw->SetMarkerColor(kBlue);
    gRaw_bmw->SetLineColor(kBlue);
    gRaw_bmw->SetMarkerStyle(7);
    mgRaw->Add(gRaw_bmw,"lp");
    legRaw->AddEntry(gRaw_bmw,Form("Raw: coil%d",coil_array[ic]));
    
    mgRaw->SetTitle(Form("%s Raw Width(ppm): BMW Coil%d Test",det_name.Data(), coil_array[ic]));
    mgRaw->Draw("A");
    mgRaw->GetXaxis()->SetTitle("Slug Number");
    mgRaw->GetYaxis()->SetTitle(" RMS (ppm)");
    mgRaw->GetYaxis()->SetTitleSize(0.05);
    legRaw->Draw("same");
    gPad->Update();
    c1->Print(Form("%s_coil_rms_by_slug.pdf",det_name.Data()));
  }// end of coil loop
  c1->Print(Form("%s_coil_rms_by_slug.pdf]",det_name.Data()));
}

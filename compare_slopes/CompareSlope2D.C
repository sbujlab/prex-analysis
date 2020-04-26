#include "GraphDitheringSlope2D.C"
// #include "GraphDitheringSlope12X2D.C"
// void CompareSlope2D(Int_t slug);
// void CompareSlope2D(){
//   for(int i=1;i<=94;i++)
//     CompareSlope2D(i);
// }
void CompareSlope2D(Int_t slug=90){
  TFile * input;
  if(slug>=4)
    input= TFile::Open(Form("merged_rootfiles/prexPrompt_slug%d_cbpm.root",slug));
  else
    input= TFile::Open(Form("merged_rootfiles/prexPrompt_slug%d.root",slug));
  TTree *reg = (TTree*)input->Get("T");
  reg->SetMarkerStyle(20);
  TString det_name[2]={"usl","usr"};
  vector<pair<TString,TString> > fBPMPair={ {"4aX","4eX"},
					    {"4aX","E"},
					    {"4eX","E"},
					    {"4aY","4eY"}};
				     
  TCanvas *c1 =new TCanvas("c1","c1",1300,600);
  c1->Divide(2,1);
  c1->Print(Form("./plots/slug%d_compare_slope2d.pdf[",slug),"pdf");
  for(int idet=0;idet<2;idet++){
    for(int iset=0;iset<4;iset++){
      TString bpm1 = fBPMPair[iset].first;
      TString bpm2 = fBPMPair[iset].second;
      if(slug<=3 && bpm2=="E")
	bpm2 = "12X";
      
      c1->cd(1);
      gPad->SetLeftMargin(0.15);
      reg->Draw(Form("diff_bpm%s/um:diff_bpm%s/um",
		     bpm1.Data(),bpm2.Data()));

      c1->cd(2);
      gPad->SetRightMargin(0.15);
      reg->Draw(Form("%s_bpm%s/(ppm/um):%s_bpm%s/(ppm/um):run",
		     det_name[idet].Data(),bpm1.Data(),
		     det_name[idet].Data(),bpm2.Data()),"","goff");
      TGraph *g_reg = new TGraph(reg->GetSelectedRows(),
				 reg->GetV2(),
				 reg->GetV1());
      Double_t *fRun = reg->GetV3();
      Int_t npt = reg->GetSelectedRows();
      g_reg->SetMarkerStyle(20);
      TLegend *leg = new TLegend(0.8,0.8,0.99,0.99);
      TMultiGraph *mg =new TMultiGraph();
      mg->Add(g_reg,"P");
      leg->AddEntry(g_reg,"regression","p");
      
      TString ch1 = Form("%s_bpm%s",det_name[idet].Data(),bpm1.Data());
      if(bpm2=="E")
	bpm2 = "11X12X";
      TString ch2 = Form("%s_bpm%s",det_name[idet].Data(),bpm2.Data());
      TGraph *g_dit = GraphDitheringSlope2D(slug,ch1,ch2);
      mg->Add(g_dit,"P");
      leg->AddEntry(g_dit,"dithering","p");
      mg->Draw("A");
      mg->SetTitle(Form("%s  vs %s; (ppm/um) ; (ppm/um)",
			ch1.Data(),ch2.Data()));
      leg->Draw("same");
      c1->Print(Form("./plots/slug%d_compare_slope2d.pdf",slug));
    }
  }
  c1->Print(Form("./plots/slug%d_compare_slope2d.pdf]",slug));
}



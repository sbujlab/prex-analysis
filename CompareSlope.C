#include "GraphDitheringSlope.C"
#include "GraphDitheringSlope12X.C"
void CompareSlope(Int_t slug);
void CompareSlope(){
  for(int i=1;i<=94;i++)
    CompareSlope(i);
}
void CompareSlope(Int_t slug){
  TFile * input;
  if(slug>=4)
    input= TFile::Open(Form("merged_rootfiles/prexPrompt_slug%d_cbpm.root",slug));
  else
    input= TFile::Open(Form("merged_rootfiles/prexPrompt_slug%d.root",slug));
  TTree *reg = (TTree*)input->Get("T");

  TString det_name[2]={"usl","usr"};
  TString bpm_name[5]={"4aX","4eX","4aY","4eY","E"};
  if(slug<=3)
    bpm_name[4] ="12X";
  int idet=0;
  int ibpm=0;
  Int_t npt;
  TCanvas *c1 = new TCanvas("c1","c1",1200,500);
  c1->SetGridx();
  c1->cd();
  c1->Print(Form("slug%d_compare_slope.pdf[",slug),"pdf");
  for(int idet=0;idet<2;idet++){
    for(int ibpm=0;ibpm<5;ibpm++){
      c1->Clear("D");

      npt = reg->Draw(Form("%s_bpm%s*1e3:Entry$:run",det_name[idet].Data(),bpm_name[ibpm].Data()),
		      "","goff");
      TGraph *g_reg = new TGraph(reg->GetSelectedRows(),
				 reg->GetV2(),reg->GetV1());
      g_reg->SetMarkerStyle(20);
      g_reg->SetMarkerColor(kBlue);
      g_reg->SetMarkerSize(2);
      TMultiGraph *mg = new TMultiGraph();
      Double_t *fRun = reg->GetV3();
      TString chname;
      if(bpm_name[ibpm]=="E")
	chname= Form("%s_bpm11X12X",det_name[idet].Data());
      else
	chname= Form("%s_bpm%s",det_name[idet].Data(),bpm_name[ibpm].Data());
      
      TMultiGraph *g_dit;
      g_dit= GraphDitheringSlope(fRun,npt,slug,chname);
      mg->Add(g_reg,"p");
      mg->Add(g_dit,"l");
      mg->Draw("AP");

      TH1D *htreg  = (TH1D*)mg->GetHistogram();
      htreg->GetXaxis()->Set(npt,-0.5,npt-0.5);
      for(int i=0;i<npt;i++){
	if(i==0)
	  htreg->GetXaxis()->SetBinLabel(i+1,Form("%0.f",fRun[i]));
	else if( fRun[i]!=fRun[i-1])
	  htreg->GetXaxis()->SetBinLabel(i+1,Form("%0.f",fRun[i]));
	else
	  htreg->GetXaxis()->SetBinLabel(i+1,"");
      }
      htreg->GetXaxis()->SetLabelSize(0.05);
      htreg->SetTitle(chname+";run number;slope(ppm/um)");
      double ymax = htreg->GetYaxis()->GetXmax();
      double ymin = htreg->GetYaxis()->GetXmin();
      htreg->GetYaxis()->SetRangeUser(ymin,ymax+0.5*(ymax-ymin));

      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry(g_reg,"Regression","p");
      TIter next(g_dit->GetListOfGraphs());
      leg->AddEntry((TGraph*)next(),"dithering 5x5 slug average","l");
      leg->Draw("same");

      c1->Print(Form("slug%d_compare_slope.pdf",slug));
    }
  }
  c1->Print(Form("slug%d_compare_slope.pdf]",slug));
}



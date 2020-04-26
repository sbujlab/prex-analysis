using namespace std;
bool compare_pair(const pair<int,int> &a, const pair<int,int> &b){
  return (a.first <  b.first);
}
void Compare(TString user_cut);
void Compare();

void Compare(){
  TCanvas *c1 = new TCanvas("c1","c1",900,600);
  c1->cd();
  for(int i=1;i<=94;i++){
    TString cut = Form("&& slug==%d",i);
    Compare(cut);
    c1->SaveAs(Form("us_avg_width_slug%d.png",i));
  }

  Compare(" && (slug==38  || slug==39)");
  c1->SaveAs("us_avg_width_slug38_39.png");
}

void Compare(TString user_cut){
  TFile *dit_file = TFile::Open("merged_dit_avg.root");
  // TFile *dit_file = TFile::Open("merged_dit.root");
  TFile *reg_file = TFile::Open("merged_reg.root");

  TTree* dit_tree = (TTree*)dit_file->Get("cor_avg");
  // TTree* dit_tree = (TTree*)dit_file->Get("cor");
  TTree* reg_tree = (TTree*)reg_file->Get("reg");

  TGraph* g_dit;
  TGraph* g_reg;
  Int_t npt;
  npt =dit_tree->Draw("dit_us_avg.rms/ppm:run","arm==0 && good"+user_cut,"goff");
  Double_t* y_ptr = dit_tree->GetV1();
  Double_t* x_ptr= dit_tree->GetV2();
  Double_t* y_dit = new Double_t[npt];
  Double_t* x_dit = new Double_t[npt];
  vector< pair<int,int> >  fDitIndex;
  for(int i=0;i<npt;i++){
    fDitIndex.push_back(make_pair(x_ptr[i],i));
  }
  sort(fDitIndex.begin(),fDitIndex.end(),compare_pair);
  int count=0;
  auto it_idx = fDitIndex.begin();
  while(it_idx!=fDitIndex.end()){
    int myindex = (*it_idx).second;
    int myrun = (*it_idx).first;
    x_dit[count]=myrun;
    y_dit[count]=y_ptr[myindex];
    count++;
    it_idx++;
  }
  g_dit = new TGraph(npt,x_dit,y_dit);
  g_dit->SetMarkerStyle(20);
  g_dit->SetMarkerColor(kBlue);
  g_dit->SetLineColor(kBlue);

  npt=reg_tree->Draw("reg_us_avg.rms/ppm:run","arm==0 && good"+user_cut,"goff");
  Double_t* y_ptr1 = (Double_t*)reg_tree->GetV1();
  Double_t* x_ptr1=(Double_t*)reg_tree->GetV2();
  Double_t* y_reg = new Double_t[npt];
  Double_t* x_reg = new Double_t[npt];
  vector< pair<int,int> >  fRegIndex;
  for(int i=0;i<npt;i++){
    fRegIndex.push_back(make_pair(x_ptr1[i],i));
  }
  sort(fRegIndex.begin(),fRegIndex.end(),compare_pair);

  auto it_reg = fRegIndex.begin();
  count=0;
  while(it_reg!=fRegIndex.end()){
    int myindex = (*it_reg).second;
    int myrun = (*it_reg).first;
    x_reg[count]=myrun;
    y_reg[count]=y_ptr1[myindex];
    it_reg++;
    count++;
  }
  g_reg = new TGraph(npt,x_reg,y_reg);
  g_reg->SetMarkerStyle(20);
  g_reg->SetMarkerSize(1.1);
  g_reg->SetMarkerColor(28);
  g_reg->SetLineColor(kOrange);
  g_reg->SetFillColor(kOrange);
  g_reg->SetFillStyle(3354);
  g_reg->SetLineWidth(-9903);
  g_reg->SetLineStyle(9);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(g_dit,"lp");
  mg->Add(g_reg,"lp");
  mg->Draw("A");
  mg->GetYaxis()->SetTitle("Us Avg Width(ppm)");
  mg->GetXaxis()->SetTitle("run number");
  TLegend *leg = new TLegend(0.8,0.85,1.0,1.0);
  leg->AddEntry(g_dit,"Dit. Corrected","lp");
  leg->AddEntry(g_reg,"Regression Floor","lpf");
  leg->Draw("same");
}



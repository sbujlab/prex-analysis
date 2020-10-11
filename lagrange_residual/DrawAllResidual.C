/*
   author: Tao Ye
   Last update: Aug 2019
 */

void DrawAllResidual(Int_t kSwitch){
  // kSwitch: 1-Lagrange, 2-Regression
  TString tag, full_tag;
  if(kSwitch==1){
    tag="lagr";
    full_tag = "lagrange";
  }else if(kSwitch==2){
    tag="reg";
    full_tag = "regression";
  }
    
  gStyle->SetOptStat(111111);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.35);
  TChain *res = new TChain("res");
  for(int i=1;i<=94;i++){
    if(i==7)
      continue;
    res->Add(Form("rootfiles/slug%d_%s_residual.root",i,full_tag.Data()));
  }
  FILE *table = fopen(Form("table_%s.txt",tag.Data()),"w");
  FILE *table_red = fopen(Form("table_red_%s.txt",tag.Data()),"w");
  // TString det_name[]={"usl","usr","us_avg","us_dd"};
  TString det_name[]={"us_avg","us_dd"};
  Int_t ndet = sizeof(det_name)/sizeof(*det_name);
  for(int idet =0;idet<ndet;idet++){
    TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
    c1->Divide(1,3);
    c1->cd(1);
    TPad *pad1l = new TPad("p1l","",0,0,0.5,1.0);
    TPad *pad1r = new TPad("p1r","",0.5,0,0.75,1.0);
    TPad *pad1rr = new TPad("p1rr","",0.75,0,1.0,1.0);
    pad1l->Draw();
    pad1r->Draw();
    pad1rr->Draw();

    c1->cd(2);
    TPad *pad2l = new TPad("p2l","",0,0,0.5,1.0);
    TPad *pad2r = new TPad("p2r","",0.5,0,0.75,1.0);
    TPad *pad2rr = new TPad("p2rr","",0.75,0,1.0,1.0);
    pad2l->Draw();
    pad2r->Draw();
    pad2rr->Draw();
    c1->cd(3);
    TPad *pad3l = new TPad("p3l","",0,0,0.5,1.0);
    TPad *pad3r = new TPad("p3r","",0.5,0,0.75,1.0);
    TPad *pad3rr = new TPad("p3rr","",0.75,0,1.0,1.0);
    pad3l->Draw();
    pad3r->Draw();
    pad3rr->Draw();

    c1->Print(Form("./plots/residual_%s_all_run_%s.pdf[",
		   full_tag.Data(),det_name[idet].Data()));
    
    TString arm_cut ="";
    if(det_name[idet].Contains("l"))
      arm_cut = "arm_flag != 1";
    if(det_name[idet].Contains("r"))
      arm_cut = "arm_flag != 2";
    if(det_name[idet].Contains("avg") || det_name[idet].Contains("dd"))
      arm_cut = "arm_flag==0";
    
    // Piggy-back runflag filter here
    arm_cut += "&&kGood==1";
    TVirtualPad *vpad;
    int npt, nptred;

    TString redundant_cut[7]={ "(run>=4735)","(run<4439 || run>4450)",
			       "(0)","(0)",
			       "(run<4745)","(run>=4439 && run<=4450)",
			       "(0)"};
    for(int icoil=1;icoil<=7;icoil++){
      
      pad1l->cd();
      TMultiGraph *mgsens = new TMultiGraph();
      TLegend *leg = new TLegend(0.85,0.8,1.0,0.95);
      npt = res->Draw(Form("%s_coil%d_sens*1e6:run",det_name[idet].Data(),icoil),
		arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&!"+redundant_cut[icoil-1],"*");
      double *y2 = res->GetV1();
      double *x2 = res->GetV2();
      TGraph *g2 = new TGraph(npt,x2,y2);
      g2->SetMarkerStyle(20);
      mgsens->Add(g2,"p");

      nptred = res->Draw(Form("%s_coil%d_sens*1e6:run",det_name[idet].Data(),icoil),
			 arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&"+redundant_cut[icoil-1],"*");
      double *y2red = res->GetV1();
      double *x2red = res->GetV2();
      if(nptred>0){
	TGraph *g2red = new TGraph(nptred,x2red,y2red);
	g2red->SetMarkerStyle(4);
	g2red->SetMarkerColor(kRed);
	mgsens->Add(g2red,"p");
	leg->AddEntry(g2red,"Redundant","p");
      }
      mgsens->Draw("AP");
      leg->Draw("same");
      mgsens->SetTitle( Form("%s_coil%d Sensitivity (ppm/count); Supercycle # ; (ppm/count)",
			 det_name[idet].Data(),icoil));
      

      // Residual Sensitivity
      pad2l->cd();
      TMultiGraph *mgres = new TMultiGraph();
      npt = res->Draw(Form("%s_coil%d_res*1e6:run",det_name[idet].Data(),icoil),
		      arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&!"+redundant_cut[icoil-1],"goff");
      double *y1 = res->GetV1();
      double *x1 = res->GetV2();
      TGraph *g1 = new TGraph(npt,x1,y1);
      g1->SetMarkerStyle(20);
      mgres->Add(g1,"p");

      // redundant coil
      nptred = res->Draw(Form("%s_coil%d_res*1e6:run",det_name[idet].Data(),icoil),
		      arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&"+redundant_cut[icoil-1],"goff");
      double *y1red = res->GetV1();
      double *x1red = res->GetV2();
      if(nptred>0){
	TGraph *g1red = new TGraph(nptred,x1red,y1red);
	g1red->SetMarkerStyle(4);
	g1red->SetMarkerColor(kRed);
	mgres->Add(g1red,"p");
      }
      mgres->Draw("AP");
      mgres->SetTitle( Form("%s_coil%d Residual (ppm/count); Supercycle # ; (ppm/count)",
			    det_name[idet].Data(),icoil));

      // Drawing histograms

      pad1r->cd();

      double y_sens_max= mgsens->GetYaxis()->GetXmax();
      double y_sens_min= mgsens->GetYaxis()->GetXmin();
      TH1D *hsens_dit = new TH1D("hsens_dit","",100,y_sens_min, y_sens_max+0.3*(y_sens_max-y_sens_min));
      res->Draw(Form("%s_coil%d_sens*1e6>>hsens_dit",det_name[idet].Data(),icoil),
		arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&!"+redundant_cut[icoil-1],"");
      if(hsens_dit!=NULL){
	hsens_dit->SetName("Dit Raw");
	hsens_dit->SetTitle("Raw Sensitivity (ppm/count)");
      }

      pad1rr->cd();
      TH1D *hsens_rd = new TH1D("hsens_rd","",100,y_sens_min, y_sens_max+0.3*(y_sens_max-y_sens_min));;
      res->Draw(Form("%s_coil%d_sens*1e6>>hsens_rd",det_name[idet].Data(),icoil),
		arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&"+redundant_cut[icoil-1],"");
      if(hsens_rd!=NULL){
	hsens_rd->SetName("Redundant Raw");
	hsens_rd->SetTitle("Raw Sensitivity (ppm/count)");
	hsens_rd->SetLineColor(kRed);
      }

      pad2r->cd();
      double ymax= mgres->GetYaxis()->GetXmax();
      double ymin= mgres->GetYaxis()->GetXmin();
      TH1D *hres = new TH1D("hres","",100,ymin,ymax+0.3*(ymax-ymin));
      res->Draw(Form("%s_coil%d_res*1e6>>hres",det_name[idet].Data(),icoil),
		arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&& !"+redundant_cut[icoil-1],"");
      hres->SetTitle("Residual Sensitivity (ppm/count)");
      hres->SetName("Dit");

      pad2rr->cd();
      TH1D *hred = new TH1D("hred","",100,ymin,ymax+0.3*(ymax-ymin));
      res->Draw(Form("%s_coil%d_res*1e6>>hred",det_name[idet].Data(),icoil),
		arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&"+redundant_cut[icoil-1],"");
      hred->SetTitle("Residual Sensitivity (ppm/count)");
      hred->SetName("Redundant");
      hred->SetLineColor(kRed);

      // Fractional Sensitivity
      pad3l->cd();
      TMultiGraph *mgfrac = new TMultiGraph();
      npt=res->Draw(Form("%s_coil%d_res/%s_coil%d_sens:run",
		     det_name[idet].Data(),icoil,det_name[idet].Data(),icoil),
		arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&!"+redundant_cut[icoil-1],"*");
      double *y3 = res->GetV1();
      double *x3 = res->GetV2();
      TGraph *g3 = new TGraph(npt,x3,y3);
      g3->SetMarkerStyle(20);
      mgfrac->Add(g3,"p");

      nptred=res->Draw(Form("%s_coil%d_res/%s_coil%d_sens:run",
			    det_name[idet].Data(),icoil,det_name[idet].Data(),icoil),
		       arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&"+redundant_cut[icoil-1],"*");
      double *y3red = res->GetV1();
      double *x3red = res->GetV2();

      if(nptred>0){
	TGraph *g3red = new TGraph(nptred,x3red,y3red);
	g3red->SetMarkerStyle(4);
	g3red->SetMarkerColor(kRed);
	mgfrac->Add(g3red,"p");
      }
      mgfrac->Draw("AP");
      mgfrac->SetTitle( Form("%s_coil%d Fractional Residual; Supercycle # ; Fraction",
			 det_name[idet].Data(),icoil));
      pad3r->cd();

      // TH1D *hgauge = new TH1D("hgauge","Fractional Residual",100,-1,1);
      // res->Draw(Form("%s_coil%d_res/%s_coil%d_sens>>hgauge",
      // 		     det_name[idet].Data(),icoil,det_name[idet].Data(),icoil),
      // 			  arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil),"");
      // double approx_width =hgauge->GetRMS();
      // hgauge->SetName("");
      // TH1D *hfrac = new TH1D("hfrac","Fractional Residual",100,-3*approx_width,3*approx_width);
      TH1D *hfrac = new TH1D("hfrac","Fractional Residual",100,-1,1);
      res->Draw(Form("%s_coil%d_res/%s_coil%d_sens>>hfrac",
		     det_name[idet].Data(),icoil,det_name[idet].Data(),icoil),
			  arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&!"+redundant_cut[icoil-1],"");
      hfrac->SetName("Dit");

      pad3rr->cd();
      // TH1D *hfrac_red = new TH1D("hfrac_red","Fractional Residual",100,-3*approx_width,3*approx_width);
      TH1D *hfrac_red = new TH1D("hfrac_red","Fractional Residual",100,-1,1);
      res->Draw(Form("%s_coil%d_res/%s_coil%d_sens>>hfrac_red",
		     det_name[idet].Data(),icoil,det_name[idet].Data(),icoil),
		arm_cut+Form("&& %s_coil%d_flag",det_name[idet].Data(),icoil)+"&&"+redundant_cut[icoil-1],"");
      hfrac_red->SetName("Redundant");
      hfrac_red->SetLineColor(kRed);
      TString my_title = det_name[idet];
      my_title.ReplaceAll("_"," ");
      fprintf(table," %s vs coil %d & ",my_title.Data(), icoil);
      fprintf(table_red," %s vs coil %d & ",my_title.Data(), icoil);
      TH1D* harray[]={hres, hsens_dit, hred, hsens_rd};
      // dithering coils
      if(hres->GetEntries()>0){
	if(fabs(hres->GetMean())<0.01)
	  fprintf(table,"%.2f & %4.1e & %.2f \\\\ \n",
		  hsens_dit->GetMean(),hres->GetMean(),hres->GetRMS());
	else
	  fprintf(table,"%.2f & %.2f & %.2f \\\\ \n",
		  hsens_dit->GetMean(),hres->GetMean(),hres->GetRMS());
      }else{
	fprintf(table," &  \\\\ \n");
      }

      // redandunt coils
      if(hred->GetEntries()>0){
	if(fabs(hred->GetMean())<0.01)
	  fprintf(table_red,"%.2f & %4.1e & %.2f \\\\ \n",
		  hsens_rd->GetMean(),hred->GetMean(),hred->GetRMS());
	else
	  fprintf(table_red,"%.2f & %.2f & %.2f \\\\ \n",
		  hsens_rd->GetMean(),hred->GetMean(),hred->GetRMS());
      }else{
	fprintf(table_red," &  \\\\ \n");
      }

      c1->Print(Form("./plots/residual_%s_all_run_%s.pdf",
		     full_tag.Data(),det_name[idet].Data()));
    }
    c1->Print(Form("./plots/residual_%s_all_run_%s.pdf]",
		   full_tag.Data(),det_name[idet].Data()));
  }
  fclose(table);
  fclose(table_red);
}

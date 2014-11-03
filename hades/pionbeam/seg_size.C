void seg_size() {

  const int ncase= 4;
  TFile *f[ncase];
  string fName[ncase] = {
    "output/report.v2/ms10_ds11_xy0.5mm_dth10_dph50_p1.30.root_graphs.root",
    "output/report.v2/ms10_ds11_seg7.0_xy0.5mm_dth10_dph50_p1.30.root_graphs.root",
    "output/report.v2/ms10_ds11_seg3.5_xy0.5mm_dth10_dph50_p1.30.root_graphs.root",
    //"output/report.v2/ms10_ds11_seg1.0_xy0.5mm_dth10_dph50_p1.30.root_graphs.root",
    "output/report.v2/ms10_ds10_xy0.5mm_dth10_dph50_p1.30.root_graphs.root"
  };

  TGraphErrors *tge_pres[ncase+1];
  TGraphErrors *tge_xres[ncase+1];
  TGraphErrors *tge_yres[ncase+1];

  cout << "Done reading in tgraphs" << endl;
  for (int ii=0; ii<ncase; ++ii) {
    cout << "File : " << fName[ii] << endl;
    f[ii] = TFile::Open(fName[ii].c_str());
    if (ii==0) {
      tge_pres[ncase] = (TGraphErrors*) f[ii]->Get("tge_pres_TDR");
      tge_xres[ncase] = (TGraphErrors*) f[ii]->Get("tge_xres_had_TDR");
      tge_yres[ncase] = (TGraphErrors*) f[ii]->Get("tge_yres_had_TDR");
    }
    tge_pres[ii] = (TGraphErrors*) f[ii]->Get("tge_pres_MINUIT");
    tge_xres[ii] = (TGraphErrors*) f[ii]->Get("tge_xres_had_MINUIT");
    tge_yres[ii] = (TGraphErrors*) f[ii]->Get("tge_yres_had_MINUIT");
  }

  cout << "Done reading in tgraphs" << endl;

  string sLeg_pres[ncase] = {
    "MINUIT, Seg(DIA)= 14mm",
    "MINUIT, Seg(DIA)= 7mm",
    "MINUIT, Seg(DIA)= 3.5mm",
    //"MINUIT, Seg(DIA)= 1mm",
    "MINUIT, Seg(DIA)= 0mm",
    "TDR"
  };

  TLegend *tl_pres = new TLegend(.4,.5,.9,.9);
  tl_pres->SetBorderSize(0);
  tl_pres->SetFillStyle(0);
  TMultiGraph *tmg_pres = new TMultiGraph();
  tmg_pres->SetTitle("Momentum Resolution; mom offset(%%); mom res (%%)");
  for (int i=0; i<ncase+1; ++i) {
    if (i<ncase) tge_pres[i]->SetMarkerStyle(20+i);
    tmg_pres->Add(tge_pres[i],"p");
    tl_pres->AddEntry(tge_pres[i],sLeg_pres[i].c_str(),"pl");
  }

  TMultiGraph *tmg_xres = new TMultiGraph();
  tmg_xres->SetTitle("X Pos Resolution (mm); mom offset(%%); Position resolution (mm)");
  for (int i=0; i<ncase+1; ++i) {
    if (i<ncase) tge_xres[i]->SetMarkerStyle(20+i);
    tmg_xres->Add(tge_xres[i],"p");
  }

  TMultiGraph *tmg_yres = new TMultiGraph();
  tmg_yres->SetTitle("Y Pos Resolution (mm); mom offset(%%); Position resolution (mm)");
  for (int i=0; i<ncase+1; ++i) {
    if (i<ncase) tge_yres[i]->SetMarkerStyle(20+i);
    tmg_yres->Add(tge_yres[i],"p");
  }


  TCanvas *tc_res = new TCanvas("tc_res_var","tc_res_var",1200,600);
  tc_res->Divide(3,1);

  tc_res->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  tmg_pres->Draw("a");
  tmg_pres->SetMinimum(0.0);
  tl_pres->Draw();

  tc_res->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  tmg_xres->Draw("a");
  tmg_xres->SetMinimum(0.0);
  //tmg_xres->SetMaximum(max_xyres*1.1);
  tmg_xres->SetMaximum(20);
  //tl_xres->Draw();

  tc_res->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  tmg_yres->Draw("a");
  tmg_yres->SetMinimum(0.0);
  //tmg_yres->SetMaximum(max_xyres*1.1);
  tmg_yres->SetMaximum(20);
  //tl_yres->Draw();

  tc_res->ls();
  double xl=0, xu=0, yl=0, yu=0;
  TPad *pad1 = (TPad*) tc_res->FindObject("tc_res_var_1");
  pad1->GetPadPar(xl,yl,xu,yu);
  cout << "pad1: xl= " << xl << "xu= " << xu << "yl= " << yl  << "yu= " << yu  << endl;
  xu= 0.44;
  pad1->SetPad(xl,yl,xu,yu);
  TPad *pad2 = (TPad*) tc_res->FindObject("tc_res_var_2");
  pad2->GetPadPar(xl,yl,xu,yu);
  cout << "pad2: xl= " << xl << "xu= " << xu << "yl= " << yl  << "yu= " << yu  << endl;
  xl= 0.45; xu=0.735;
  pad2->SetPad(xl,yl,xu,yu);
  TPad *pad3 = (TPad*) tc_res->FindObject("tc_res_var_3");
  pad3->GetPadPar(xl,yl,xu,yu);
  cout << "pad3: xl= " << xl << "xu= " << xu << "yl= " << yl  << "yu= " << yu  << endl;
  xl= 0.735; xu=0.99;
  pad3->SetPad(xl,yl,xu,yu);

  double eps = 1e-9;
  pad3->SetLeftMargin(eps);
  pad3->SetRightMargin(0.1);
  pad2->SetLeftMargin(0.15);
  pad2->SetRightMargin(eps);
  pad1->SetRightMargin(0.05);
  pad1->SetLeftMargin(0.15);

  tc_res->Print(Form("slides/tc_res_var.eps"));


}

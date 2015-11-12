void bg_xsect() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  ifstream inf;
  inf.open("figs/bg_xsect_data.txt");
  const int npt = 26;
  double ecm[npt], ecm_er[npt], _plab[npt], _plab_er[npt], xsect[npt], xsect_er[npt];
  int ctr = 0, ipt;
  while(true) {
    inf >> ipt;
    if (!inf.good() || ipt!= ctr || ipt > 24) break;
    inf >> ecm[ipt] >> _plab[ipt] >>  xsect[ipt] >> xsect_er[ipt];
    ecm_er[ipt] = _plab_er[ipt] = 0.0;
    cout << "ipt= " << ipt << " " << ecm[ipt] << " " << _plab[ipt] << " " << xsect[ipt] << " " << xsect_er[ipt] << endl;
    ++ctr;
  }

  ecm[25] = 10;
  _plab[25] = 14;
  xsect[25] = 0;
  xsect_er[25] = 0;

  double ixsect[nplab] = {0.0};
  double ixsect_er[nplab] = {0.0};
  for (int iplab=0; iplab < nplab; ++iplab) {
    bool inter = false;
    for (int jplab=1; jplab<24; ++jplab) {
      if (_plab[jplab-1] < plab[iplab] && plab[iplab] < _plab[jplab]) {
	ixsect[iplab] = xsect[jplab-1] + (xsect[jplab]-xsect[jplab-1])*(plab[iplab]-_plab[jplab-1])/(_plab[jplab]-_plab[jplab-1]);
	ixsect_er[iplab] = pow(_plab[jplab]-_plab[jplab-1], 2)/8;
	inter = true;
      }
    }
    if (!inter) {
      ixsect[iplab] = xsect[jplab-2] + (plab[iplab]-_plab[jplab-2])*(xsect[jplab-1]-xsect[jplab-2])/(_plab[jplab-1]-_plab[jplab-2]);
      ixsect_er[iplab] = pow(_plab[jplab-1]-_plab[jplab-2], 2)/8;
    }
  }

  TGraphErrors *tge_ecm= new TGraphErrors(npt, ecm,xsect,ecm_er,xsect_er);
  tge_ecm->SetTitle(";E_{cm}[GeV];#sigma(#bar{p}+p#rightarrow#pi^{0}#pi^{+}#pi^{-})[mb]");
  set_style(tge_ecm, 2);

  TGraphErrors *tge_plab= new TGraphErrors(npt, _plab,xsect,_plab_er,xsect_er);
  tge_plab->SetTitle(";p_{#bar{p}} [GeV];#sigma(#bar{p}+p#rightarrow#pi^{0}#pi^{+}#pi^{-}) [mb]");
  set_style(tge_plab, 2, 20, 1, 2);

  TCanvas *tc_ecm = new TCanvas("tc_ecm","tc_ecm");
  tc_ecm->cd();
  tge_ecm->Draw("ap");
  tge_ecm->SetMaximum(12);
  gPad->SetLogy();
  tc_ecm->Update();

  TCanvas *tc_plab = new TCanvas("tc_plab","tc_plab");
  tc_plab->cd();
  tge_plab->Draw("ap");
  gPad->SetLogy();
  gPad->Update();
  TH1F* dum = (TH1F*)tge_plab->GetHistogram();
  set_style(dum);
  dum->SetMaximum(20);
  dum->GetXaxis()->SetRangeUser(0,14);
  tc_plab->GetPad(0)->SetBottomMargin(0.15);

  //TF1 *fpl = new TF1("mypl","[0]*(x**[1])",3.0,14);
  //fpl->SetParameter(0,10);
  //fpl->SetParameter(1,-1);
  //tge_plab->Fit(fpl,"r");

  TLine *tl[3],*tl2[3];
  double bg_xsect[nplab] = {0.2, 0.05, 0.02};
  //for (int iplab=0; iplab < nplab; ++iplab) {
  //  cout << "xsect (plab=" << plab[iplab] << ") = " << ixsect[iplab] << " \pm " << ixsect_er[iplab] << endl;
  //  bg_xsect[iplab] = ixsect[iplab];
  //}
  //for (int iplab=0; iplab < nplab; ++iplab) {
  //  bg_xsect[iplab] = fpl->Eval(plab[iplab]);
  //  cout << "xsect (plab=" << plab[iplab] << ") = " << bg_xsect[iplab] << " \pm " << endl;
  //}

  for (int iplab=0; iplab<nplab; ++iplab) {
    tl[iplab] = new TLine();
    tl[iplab]->SetLineWidth(2);
    tl[iplab]->SetLineColor(1);
    tl[iplab]->SetLineStyle(9);
    cout << "plab = " << plab[iplab] << endl;
    tl[iplab]->DrawLine(plab[iplab], tge_plab->GetMinimum(), plab[iplab], bg_xsect[iplab]);

    tl2[iplab] = new TLine();
    tl2[iplab]->SetLineWidth(2);
    tl2[iplab]->SetLineStyle(9);
    tl2[iplab]->SetLineColor(1);
    cout << "plab = " << bg_xsect[iplab] << endl;
    tl2[iplab]->DrawLine(0, bg_xsect[iplab], plab[iplab], bg_xsect[iplab]);

  }

  TLegend *leg = new TLegend(0.2, 0.74, 0.9, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(tge_plab,"World data on #sigma(#bar{p}+p#rightarrow#pi^{0}#pi^{+}#pi^{-})", "pl");
  leg->AddEntry(tl[0],"Cross-sections used for background MC","l");
  leg->Draw();

  tc_plab->Print(Form("%s/figs/2015.09.15/bg_xsect_world.pdf",bdir));
}

void nbrem() {

  /*
     |----+------+--------+--------|
     | id |   pT | th_min | th_max |
     |----+------+--------+--------|
     |  0 | 0.15 |   30.0 |   45.0 |
     |  1 | 0.25 |   30.0 |   45.0 |
     |  2 |  0.5 |   30.0 |   45.0 |
     |  3 |  0.5 |   45.0 |   60.0 |
     |  4 |  0.5 |   60.0 |   75.0 |
     |  5 |  0.5 |   75.0 |   90.0 |
     |  6 |  0.5 |   90.0 |  105.0 |
     |  7 |  0.5 |  105.0 |  120.0 |
     |  8 |  1.0 |   30.0 |   45.0 |
     |  9 |  1.0 |   45.0 |   60.0 |
     | 10 |  1.0 |   60.0 |   75.0 |
     | 11 |  1.0 |   75.0 |   90.0 |
     | 12 |  1.5 |   30.0 |   45.0 |
     | 13 |  1.5 |   45.0 |   60.0 |
     | 14 |  2.0 |   30.0 |   45.0 |
     | 15 | 0.15 |   10.0 |   20.0 |
     | 16 | 0.25 |   10.0 |   20.0 |
     | 17 |  0.5 |   10.0 |   20.0 |
     | 18 |  1.0 |   10.0 |   20.0 |
     | 19 |  1.5 |   10.0 |   20.0 |
     | 20 |  2.0 |   10.0 |   20.0 |
     | 21 |  2.5 |   10.0 |   20.0 |
     |----+------+--------+--------|
   */

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  Double_t config[22][3] = {
    {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
    {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
    {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
    {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
    {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
    {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

  //int nc = 3;
  //int ic[nc] = {0,1,2,8,12,14};
  //int col[nc] = {1,2,5,6,9,7};

  const int nc = 3;
  int ic[nc] = {1,8,14};
  int col[nc] = {1,2,4};
  TH1F *hnmcb_frac[nc];
  TCanvas *tc = new TCanvas("tc","tc");
  TLegend *tl1 = new TLegend(0.2,0.2,0.55,0.5);
  tl1->SetFillStyle(0);
  tl1->SetBorderSize(0);
  tl1->SetTextSize(0.07);

  for (int ii = 0; ii < nc; ++ii) {
    TFile *f = new TFile(Form("../grid.out/psim_oct14_binsong_configs/bremcorr.%d_hists.root",ic[ii]));
    hnmcb_frac[ii] = (TH1F*) f->Get("h_nmcb_gt1mev");
    hnmcb_frac[ii]->SetLineColor(col[ii]);
    hnmcb_frac[ii]->SetLineWidth(2);
    hnmcb_frac[ii]->GetXaxis()->SetRangeUser(0,5);
    tl1->AddEntry(hnmcb_frac[ii],Form("p_{T} = %3.2f GeV/c",config[ic[ii]][0]), "pl");
    hnmcb_frac[ii]->Draw(ii==0?"":"same");
  }
  tl1->SetHeader("    30^{o}<#theta<45^{o}");
  tl1->Draw();

}

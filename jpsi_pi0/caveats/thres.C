#include "pad_setup.C"
void set_style(TH1F* h) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  //  h->GetXaxis()->SetNdivisions(610);
}

void thres() {

  gStyle->SetOptStat(0);
  //gStyle->SetPadLeftMargin(0.13);
  //gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  const int ncfg = 6;
  TFile *f[ncfg];
  TH2F* dphdth[ncfg];
  TH1F* dth[ncfg];
  int cfg[ncfg] = {2, 8, 14, 17, 18, 20};
  int col[ncfg] = {1, 2, 4, 1, 2, 4};
  double mom[ncfg] = {0.5, 1, 2, 0.5, 1, 2};
  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(3,2);
  for (int icfg=0; icfg < ncfg; ++icfg) {
    f[icfg] = new TFile(Form("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.%d_dth_hists.root",cfg[icfg]));
    dphdth[icfg] = (TH2F*) f[icfg]->Get("h_dphi_dthe_brem");
    dth[icfg] = (TH1F*) (dphdth[icfg]->ProjectionX(Form("dth%d",icfg)));
    dth[icfg]->Rebin(2);
    dth[icfg]->GetXaxis()->SetRangeUser(-5,5);
    dth[icfg]->SetLineColor(col[icfg]);
    set_style(dth[icfg]);
    dth[icfg]->SetTitle(Form("p=%4.1f GeV/c, %s;#Delta#theta",mom[icfg], (icfg<3?"Barrel":"Endcap")));
    //dth[icfg]->Scale(1./double(dth[icfg]->GetEntries()));
    tc->cd(1+icfg);
    dth[icfg]->DrawNormalized(icfg==0?"":"same");
  }

}

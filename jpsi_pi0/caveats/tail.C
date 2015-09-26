void tail() {

  const int ncfg=6;
  TFile*f[ncfg];
  TH1F* hres[ncfg];

  f[0] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.2_full_range_hists.root");
  f[1] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.8_full_range_hists.root");
  f[2] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.14_full_range_hists.root");
  f[3] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.17_full_range_hists.root");
  f[4] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.18_full_range_hists.root");
  f[5] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.20_full_range_hists.root");

  double mom[ncfg] = {0.5, 1, 2, 0.5, 1, 2};
  TF1* fgaus[ncfg];
  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(3,2);
  for (int i=0; i<ncfg; ++i) {
    hres[i] = (TH1F*)f[i]->Get("h_rec");
    tc->cd(1+i);
    fgaus[i] = new TF1(Form("gaus%d",i),"gaus",-0.03,0.02);
    hres[i]->Fit(fgaus[i],"QRIME+");
    double dmin = fgaus[i]->GetParameter(1) - 2*fgaus[i]->GetParameter(2);
    double dmax = fgaus[i]->GetParameter(1) + 2*fgaus[i]->GetParameter(2);
    int imin = hres[i]->GetXaxis()->FindBin(dmin);
    int imax = hres[i]->GetXaxis()->FindBin(dmax);
    double in = hres[i]->Integral(imin,imax);
    double full = hres[i]->Integral();
    double frac= in*100/full;
    cout << "\\hline " << mom[i] << " & " << Form("%4.1f",fgaus[i]->GetParameter(2)*1000) << "~MeV & " << Form("%4.1f\\%% \\\\",frac) << endl;
  }

}

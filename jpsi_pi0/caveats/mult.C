void mult() {

  const int ncfg = 12;
  TFile*f[ncfg];
  TH1F* mult[ncfg];

  f[0] = new TFile("../grid.out/esim_oct14_constp_g3/bremcorr.all.ibs.cfg.8_hists.root");
  f[1] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.8_hists.root");
  f[2] = new TFile("../grid.out/esim_oct14_constp_steplimonlyg4/bremcorr.all.ibs.cfg.8_hists.root");

  f[3] = new TFile("../grid.out/esim_oct14_constp_g3/bremcorr.all.ibs.cfg.14_hists.root");
  f[4] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.14_hists.root");
  f[5] = new TFile("../grid.out/esim_oct14_constp_steplimonlyg4/bremcorr.all.ibs.cfg.14_hists.root");

  f[6] = new TFile("../grid.out/esim_oct14_constp_g3/bremcorr.all.ibs.cfg.18_hists.root");
  f[7] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.18_hists.root");
  f[8] = new TFile("../grid.out/esim_oct14_constp_steplimonlyg4/bremcorr.all.ibs.cfg.18_hists.root");

  f[9] = new TFile("../grid.out/esim_oct14_constp_g3/bremcorr.all.ibs.cfg.20_hists.root");
  f[10] = new TFile("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.20_hists.root");
  f[11] = new TFile("../grid.out/esim_oct14_constp_steplimonlyg4/bremcorr.all.ibs.cfg.20_hists.root");

  for (int icfg=0; icfg < ncfg; ++icfg ) {
    mult[icfg] = (TH1F*)f[icfg]->Get("h_nmcb_gt1mev");
  }

  double tot[ncfg]={0.};
  for (int i=1; i<4; ++i) {
    cout << "\\hline " << i-1 << " & ";
    for (int icfg =0; icfg < ncfg; ++icfg) {
      tot[icfg] += mult[icfg]->GetBinContent(i);
      //cout << Form("%4.1f", mult[icfg]->GetBinContent(i)) << "$\\%$" << (icfg<ncfg-1?" & ":" \\\\");
      cout << Form("%4.1f", mult[icfg]->GetBinContent(i)) << (icfg<ncfg-1?" & ":" \\\\");
    }
    cout <<endl;
  }
  cout <<  "\\hline $\\geq$3 & ";
  for (int icfg =0; icfg < ncfg; ++icfg) {
    //cout << Form("%4.1f",100. - tot[icfg]) << "$\\%$" << (icfg<ncfg-1?" & ":" \\\\");
    cout << Form("%4.1f",100. - tot[icfg]) << (icfg<ncfg-1?" & ":" \\\\");
  }
  cout << endl;
  cout << "\\hline" << endl;
}

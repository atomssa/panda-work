void xsect_scale_down(int pass=15){

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";
  gROOT->LoadMacro(Form("%s/figs/v2/ananote.C",bdir));

  double tvalidmin[nplab] = {-0.092, -1.0, -1.0};
  double tvalidmax[nplab] = {0.59, 0.43, 0.3};

  double xsect[nplab] = {0.2, 0.05, 0.02};
  double xsect_sig[nplab] = {206.8e-9, 280.5e-9, 199.6e-9};
  double xsect_tcut[nplab] = {0.};
  double xsect_tcut_mcut[nplab] = {0.};
  double xsect_tcut_mcut_wide[nplab] = {0.};
  double br = 0.0596;

  for (int iplab=0; iplab < nplab; ++iplab) {

    TFile *f = TFile::Open(Form("../../hists/note.v2.oct.2015/anav2_pi0pipm_brem_p%d_pass%d.root",iplab,pass));
    //TFile *f = TFile::Open(Form("../../hists/note.aug.2015/anav2_pi0pipm_brem_p%d_pass%d.root",iplab,pass));

    TH1F* h_tdist;
    TH1F* h_tot;
    TH1F* h_tcut;
    TH1F* h_tcut_mcut;

    h_tdist = (TH1F*)f->Get("tu/httrumc");

    if (pass>=14) {
      h_tot = (TH1F*)f->Get("tu/htrupi0thlab");
      h_tcut = (TH1F*)f->Get("tu/htrupi0thlab_tcut");
      h_tcut_mcut = (TH1F*)f->Get("tu/htrupi0thlab_tcut_mcut");
    } else {
      h_tot = (TH1F*)f->Get("tu/htrupi0thlab");
      h_tcut = (TH1F*)f->Get("tu/htrupi0thlab_tc");
      h_tcut_mcut = (TH1F*)f->Get("tu/htrupi0thlab_tc_mc");
    }

    //double tot = h_tot->GetEntries();
    //double pass_tcut = h_tcut->GetEntries();

    int i_tmin = h_tdist->GetXaxis()->FindBin(tvalidmin[iplab]);
    int i_tmax = h_tdist->GetXaxis()->FindBin(tvalidmax[iplab]);
    double tot = h_tdist->GetEntries();
    double pass_tcut = h_tdist->Integral(i_tmin,i_tmax);

    double frac_tcut = 2*pass_tcut/tot;
    xsect_tcut[iplab] = frac_tcut*xsect[iplab];

    double tot_tcut_wide = h_tcut->Integral();
    double tot_tcut_wide_mcut = h_tcut_mcut->Integral();
    double frac_mcut_wide = tot_tcut_wide_mcut/tot_tcut_wide;
    xsect_tcut_mcut_wide[iplab] = frac_tcut*frac_mcut_wide*xsect[iplab];


    double pass_mcut_tcut = h_tcut_mcut->GetEntries();
    double frac_mcut_tcut = pass_mcut_tcut/tot;
    xsect_tcut_mcut[iplab] = frac_mcut_tcut*xsect[iplab];

    cout << "xsect= " << xsect[iplab]
      //<< " & frac_tcut = " << frac_tcut
	 << " & xsect_tcut = " << xsect_tcut[iplab] << " s/b_tcut = " << br*xsect_sig[iplab]/xsect_tcut[iplab]
      //<< " & frac_mcut = " << frac_mcut_tcut
      // << " & xsect_tmcut = " << xsect_tcut_mcut[iplab] << " s/b_tmcut = " << br*xsect_sig[iplab]/xsect_tcut_mcut[iplab]
      //<< " & frac_mcut_wide = " << frac_mcut_wide
      //<< " & frac_tcut_mcut_wide = " << frac_tcut*frac_mcut_wide
	 << " & xsect_tmcut_wide = " << xsect_tcut_mcut_wide[iplab] << " s/b_tmcut_wide = " << br*xsect_sig[iplab]/xsect_tcut_mcut_wide[iplab]
	 << endl;

  }


}

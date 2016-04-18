double count_jpsim(TH1F* s) {
  int mmin = s->GetXaxis()->FindBin(2.8);
  int mmax = s->GetXaxis()->FindBin(3.3);
  double _sig = s->Integral(mmin,mmax);
  return _sig>0? _sig: -1;
}

double stob_r(TH1F* s, TH1F*b) {
  int mmin = s->GetXaxis()->FindBin(2.8);
  int mmax = s->GetXaxis()->FindBin(3.3);
  double _sig = s->Integral(mmin,mmax);
  double _bg = b->Integral(mmin,mmax);
  return _bg>0?_sig/_bg: -1;
}

static const double Mp = 0.938;
static const double Mp_sq = Mp*Mp;
static const double Mpi0 = 0.0; //0.135;
static const double Mpi0_sq = Mpi0*Mpi0;
static const double Mjpsi = 3.097;
static const double Mjpsi2 = Mjpsi*Mjpsi;
static const double msqTot = 2*Mp_sq + Mjpsi2 + Mpi0_sq;
double mirror(double t, double s) {
  return msqTot - s - t;
}

void kin_cuts2(int iaft = 8) {

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.1);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  gROOT->LoadMacro(Form("%s/figs/pv0/ananote.C",bdir));

  const int nbg = 4;
  int col_bg[nbg] = {30, 41, 46, 38};
  TFile *fsg[nplab], *fbg[nbg][nplab];
  TH1F *hmep_fg_b[nplab], *hmep_valid_fg_b[nplab], *hmep_sg_b[nplab], *hmep_valid_sg_b[nplab];
  TH1F *hmep_bg_b[nbg][nplab], *hmep_valid_bg_b[nbg][nplab], *hmep_bg_tot_b[nplab], *hmep_valid_bg_tot_b[nplab];
  TH1F *hmep_fg_a[nplab], *hmep_valid_fg_a[nplab], *hmep_sg_a[nplab], *hmep_valid_sg_a[nplab];
  TH1F *hmep_bg_a[nbg][nplab], *hmep_valid_bg_a[nbg][nplab], *hmep_bg_tot_a[nplab], *hmep_valid_bg_tot_a[nplab];

  TH1F *hpi0jpsi_chi24c[nbg+1][nplab], *hpi0jpsi_chi24cc[nbg+1][nplab];
  TH1F *hpi0pi0jpsi_chi24c[nbg+1][nplab], *hpi0pi0jpsi_chi24cc[nbg+1][nplab];
  TH2F *hpi0vs2pi0_chi24c[nbg+1][nplab], *hpi0vs2pi0_chi24cc[nbg+1][nplab];

  const double plab[nplab] = {5.513, 8., 12.};
  const double _s[nplab] = {12.25, 16.87, ,24.35};

  const double tvalidmin[nplab] = {-0.092, -1.0, -1.0};
  const double tvalidmax[nplab] = {0.59, 0.43, 0.3};

  double purity[4][3] = { { 90.4, 89.3, 89.1}, {98.8, 98.7, 99.7}, {98.8, 98.7, 98.8}, {99.0, 99.0, 99.0} };

  double nevt_sim_sg[3] = {32780.0,50142.0,51860.0};
  double nevt_xsect_sg[3] = {32780.0,50142.0,51860.0};
  double nevt_sim_bg[4][3] = {{814794.0,888292.0,898721.0}, {214780.0,174864.0,160099.0}, {570751.0,609044.0,527506.0}, {1.98e6,2.0e6,1.97e7}};
  double nevt_xsect_bg[4][3] = {{4.0e11, 1e11, 2e10}, {1.15e12, 3.15e11, 6.84e10}, {3.19e12, 1.14e12, 2.92e11}, {94243.0, 157947.3, 177361.2}};
  // Undo event scaling applied in anav2 module, to do the scaling anew with proper number of events
  //double pi0pi0jpsi_unscale_anav2[3] = {200000.0/94243.0,200000.0/157947.3,200000.0/177361.2};
  double pi0pi0jpsi_unscale_anav2[3] = {1.0,1.0,1.0};
  double pi0pi0jpsi_re_scale[3] = {1.0};
  for (int iplab=0; iplab < nplab; ++iplab) {
    pi0pi0jpsi_re_scale[iplab] = pi0pi0jpsi_unscale_anav2[iplab]*nevt_xsect_bg[3][iplab]/nevt_sim_bg[3][iplab];
  }

  //double nevt_sim[5][3] = {{814794.0,888292.0,898721.0}, {32780.0,50142.0,51860.0}, {214780.0,174864.0,160099.0}, {570751.0,609044.0,527506.0}, {1.98e6,2.0e6,1.97e7}};
  ////double nevt_sim[5][3] = {{814794.0,888292.0,898721.0}, {32780.0,50142.0,51860.0}, {214780.0,174864.0,160099.0}, {570751.0,609044.0,527506.0}, {2.0e5,2.0e5,1.98e5}};
  //double nevt_xsect[5][3] = {{4.0e11, 1e11, 2e10}, {32780.0,50142.0,51860.0}, {1.15e12, 3.15e11, 6.84e10}, {3.19e12, 1.14e12, 2.92e11}, {94243.0, 157947.3, 177361.2}};
  ////temporary, until fixed in effing macro
  //double pi0pi0jpsi_re_scale[3] = {1.0,1.0,1.0};
  //for (int iplab=0; iplab < nplab; ++iplab) {
  //  pi0pi0jpsi_re_scale[iplab] = 200000.0/nevt_sim[4][iplab];
  //}

  double nevt_sig_sim_valid[3] = {24561.88, 33327.48, 23742.89};
  double nevt_pi0pi0jpsi_sim_valid[3] = {1.98e6,2.0e6,1.98e6};
  //double nevt_pi0pi0jpsi_sim_valid[3] = {200000.0,200000.0,200000.0};
  //double nevt_pi0pipm_sim_valid[3] = {814794.0,888292.0,898721.0};
  double nevt_pi0pipm_sim_valid[3] = {4.0e11, 1e11, 2e10};
  double pi0pipm_frac_valid[3] = {1.0,1.0,1.0};

  //double tmp_scale[3] = {0.0};
  //for (int ii=0; ii < 3; ++ii) {
  //  tmp_scale[ii] = nevt_xsect[0][ii]/nevt_xsect[2][ii];
  //  cout << "tmp_scale" << ii << " = " << tmp_scale[ii] << endl;
  //}
  //return;
  //double tmp_scale[3] = {32780.0/200000.0,50142.0/200000.0,51860.0/200000.0};
  //double tmp_scale[3]= {94243.0/200000.0, 157947.3/200000.0, 177361.2/200000.0};

  //const char *leg_bg[nbg] = {"#pi^{0}#pi^{+}#pi^{-}","#pi^{0}#pi^{0}#pi^{+}#pi^{-}", "#pi^{0}#pi^{+}#pi^{-}#pi^{+}#pi^{-}", "#pi^{0}J/#psi(#pi^{+}#pi^{-})"};
  //const char *leg_bg[nbg] = {"2(#pi^{0})J/#psi(#rightarrow e^{+}e^{-})", "#pi^{0}#pi^{+}#pi^{-}", "2(#pi^{0})#pi^{+}#pi^{-}", "#pi^{0}2(#pi^{+}#pi^{-})"};
  const char *leg_bg[nbg] = {"#pi^{0}#pi^{0}J/#psi(#rightarrow e^{+}e^{-})", "#pi^{0}#pi^{+}#pi^{-}", "#pi^{0}#pi^{0}#pi^{+}#pi^{-}", "#pi^{0}#pi^{+}#pi^{-}#pi^{+}#pi^{-}"};
  //const double scale[nplab][nbg] = {{1.,1.,1.,0.037817}, {1.,1.,1.,0.037817}, {1.,1.,1.,0.037817}}; //  (30pb for pi0pi0jpsi )
  //const double scale[nplab][nbg] = {{1.,1.,1.,0.64}, {1.,1.,1.,0.64}, {1.,1.,1.,0.64}}; // (30pb for pi0pi0jpsi->e+e-4gamma)

  TLegend *tl = new TLegend(0.17,0.54,0.45,0.89);
  tl->SetTextSize(0.04);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);

  TLegend *tl2 = new TLegend(0.17,0.5,0.65,0.85);
  tl2->SetBorderSize(0);
  tl2->SetFillStyle(0);

  TLegend *tl3 = new TLegend(0.17,0.5,0.65,0.85);
  tl3->SetBorderSize(0);
  tl3->SetFillStyle(0);

  //TLegend *tl4 = new TLegend(0.11,0.6,0.65,0.9);
  //TLegend *tl4 = new TLegend(0.15,0.65,0.7,0.93);
  TLegend *tl4 = new TLegend(0.172052,0.635701,0.635701,0.916199);
  tl4->SetTextSize(0.06);
  tl4->SetBorderSize(0);
  tl4->SetFillStyle(1001);
  tl4->SetFillColor(0);

  int ibef = 5;

  cout << "iaft= " << iaft << endl;
  for (int iplab=0; iplab < nplab; ++iplab) {

    //int pass=14;
    //const char *hist_name = "hmep";
    //fsg[iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    //fbg[0][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    //fbg[1][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    //fbg[2][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi02pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    //fbg[3][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pipm2_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    //hmep_sg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("%s_%d",hist_name,ibef)) ->Clone(Form("hmep_sg_b_p%d",iplab));
    //hmep_sg_a[iplab] = (TH1F*) fsg[iplab]->Get(Form("%s_%d",hist_name,iaft)) ->Clone(Form("hmep_sg_a_p%d",iplab));
    //hmep_fg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("%s_%d",hist_name,ibef)) ->Clone(Form("hmep_fg_b_p%d",iplab));
    //hmep_fg_a[iplab] = (TH1F*) fsg[iplab][ii]->Get(Form("%s_%d",hist_name,iaft)) ->Clone(Form("hmep_fg_a_p%d",iplab));

    int pass=18;
    fsg[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[0][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    //fbg[0][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pi0jpsi_%s4comp_p0_pass%d.root",bdir,(ibrem==0?"raw":"brem"), pass));
    fbg[1][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[2][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi02pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[3][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pipm2_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));

    hmep_sg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",ibef)) ->Clone(Form("hmep_sg_b_p%d",iplab));
    hmep_sg_a[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",iaft)) ->Clone(Form("hmep_sg_a_p%d",iplab));
    hmep_valid_sg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_valid_%d",ibef)) ->Clone(Form("hmep_valid_sg_b_p%d",iplab));
    hmep_valid_sg_a[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_valid_%d",iaft)) ->Clone(Form("hmep_valid_sg_a_p%d",iplab));
    set_style(hmep_sg_b[iplab],2,2);
    set_style(hmep_sg_a[iplab],2,2);
    set_style(hmep_valid_sg_b[iplab],2,2);
    set_style(hmep_valid_sg_a[iplab],2,2);

    hmep_fg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",ibef)) ->Clone(Form("hmep_fg_b_p%d",iplab));
    hmep_fg_a[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",iaft)) ->Clone(Form("hmep_fg_a_p%d",iplab));
    hmep_valid_fg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_valid_%d",ibef)) ->Clone(Form("hmep_valid_fg_b_p%d",iplab));
    hmep_valid_fg_a[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_valid_%d",iaft)) ->Clone(Form("hmep_valid_fg_a_p%d",iplab));
    set_style(hmep_fg_b[iplab],1,2);
    set_style(hmep_fg_a[iplab],1,2);
    set_style(hmep_valid_fg_b[iplab],1,2);
    set_style(hmep_valid_fg_a[iplab],1,2);

    //if (iplab==0) {
    //  //hmep_fg_b[iplab]->SetTitle(Form("Before cut on 4C fit #chi^{2} (p_{#bar{p}} = %5.4g GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
    //  //hmep_fg_a[iplab]->SetTitle(Form("After cuts on 4C fit #chi^{2} (p_{#bar{p}} = %5.4g GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
    //  hmep_fg_b[iplab]->SetTitle(Form("p_{#bar{p}} = %5.4g GeV/c; M_{ee}[GeV/c^{2}]", plab[iplab]));
    //  hmep_fg_a[iplab]->SetTitle(Form("p_{#bar{p}} = %5.4g GeV/c; M_{ee}[GeV/c^{2}]", plab[iplab]));
    //} else {
    //  //hmep_fg_b[iplab]->SetTitle(Form("Before cut on 4C fit #chi^{2} (p_{#bar{p}} = %4.1f GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
    //  //hmep_fg_a[iplab]->SetTitle(Form("After cuts on 4C fit #chi^{2} (p_{#bar{p}} = %4.1f GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
    //  hmep_fg_b[iplab]->SetTitle(Form("p_{#bar{p}} = %4.1f GeV/c; M_{ee}[GeV/c^{2}]", plab[iplab]));
    //  hmep_fg_a[iplab]->SetTitle(Form("p_{#bar{p}} = %4.1f GeV/c; M_{ee}[GeV/c^{2}]", plab[iplab]));
    //}

    double dm_ee = 1000*(hmep_fg_b[iplab]->GetXaxis()->GetBinCenter(2) - hmep_fg_b[iplab]->GetXaxis()->GetBinCenter(1));

    hmep_fg_b[iplab]->SetTitle(Form("; M_{ee}[GeV/c^{2}];dN/dm_{ee}[(%3.0f MeV/c^{2})^{-1}]", dm_ee, plab[iplab]));
    hmep_fg_a[iplab]->SetTitle(Form("; M_{ee}[GeV/c^{2}]", plab[iplab]));
    hmep_valid_fg_b[iplab]->SetTitle(Form("; M_{ee}[GeV/c^{2}];dN/dm_{ee}[(%3.0f MeV/c^{2})^{-1}]", dm_ee, plab[iplab]));
    hmep_valid_fg_a[iplab]->SetTitle(Form("; M_{ee}[GeV/c^{2}]", plab[iplab]));

    if (iplab==0) {
      //tl->AddEntry(hmep_fg_a[iplab],"FG+BG","l");
      tl2->AddEntry(hmep_fg_a[iplab],"FG+BG","l");
      tl4->AddEntry(hmep_fg_a[iplab],"FG+BG","l");
      tl->AddEntry(hmep_sg_a[iplab],"#pi^{0}J/#psi(#rightarrow e^{+}e^{-})","l");
      tl2->AddEntry(hmep_sg_a[iplab],"#pi^{0}J/#psi(#rightarrow e^{+}e^{-})","l");
      tl3->AddEntry(hmep_sg_a[iplab],"#pi^{0}J/#psi(#rightarrow e^{+}e^{-})","l");
      tl4->AddEntry(hmep_sg_a[iplab],"#pi^{0}J/#psi(#rightarrow e^{+}e^{-})","l");

    }

    hpi0jpsi_chi24c[0][iplab] = (TH1F*) fsg[iplab]->Get("hpi0jpsi_chi24c")->Clone(Form("hpi0jpsi_chi24c_t0_p%d",iplab));
    hpi0jpsi_chi24cc[0][iplab] = (TH1F*) fsg[iplab]->Get("hpi0jpsi_chi24c_c")->Clone(Form("hpi0jpsi_chi24cc_t0_p%d",iplab));
    set_style(hpi0jpsi_chi24c[0][iplab],2,4);
    set_style(hpi0jpsi_chi24cc[0][iplab],2);

    hpi0pi0jpsi_chi24c[0][iplab] = (TH1F*) fsg[iplab]->Get("hpi0pi0jpsi_chi24c")->Clone(Form("hpi0pi0jpsi_chi24c_t0_p%d",iplab));
    hpi0pi0jpsi_chi24cc[0][iplab] = (TH1F*) fsg[iplab]->Get("hpi0pi0jpsi_chi24c_c")->Clone(Form("hpi0pi0jpsi_chi24cc_t0_p%d",iplab));
    set_style(hpi0pi0jpsi_chi24c[0][iplab],2,4);
    set_style(hpi0pi0jpsi_chi24cc[0][iplab],2);
    hpi0vs2pi0_chi24c[0][iplab] = (TH2F*) fsg[iplab]->Get("hpi0vs2pi0_chi24c")->Clone(Form("hpi0vs2pi0_chi24c_t0_p%d",iplab));
    hpi0vs2pi0_chi24cc[0][iplab] = (TH2F*) fsg[iplab]->Get("hpi0vs2pi0_chi24c_c")->Clone(Form("hpi0vs2pi0_chi24cc_t0_p%d",iplab));
    set_style(hpi0vs2pi0_chi24c[0][iplab]);
    set_style(hpi0vs2pi0_chi24cc[0][iplab]);

    TH1F *ttrumc = (TH1F*) fbg[1][iplab]->Get("tu/httrumc");
    double i_tmin = ttrumc->GetXaxis()->FindBin(tvalidmin[iplab]);
    double i_tmax = ttrumc->GetXaxis()->FindBin(tvalidmax[iplab]);
    double i_umin = ttrumc->GetXaxis()->FindBin(mirror(tvalidmax[iplab],_s[iplab]));
    double i_umax = ttrumc->GetXaxis()->FindBin(mirror(tvalidmin[iplab],_s[iplab]));
    pi0pipm_frac_valid[iplab] = (ttrumc->Integral(i_tmin,i_tmax)+ttrumc->Integral(i_tmin,i_tmax))/ttrumc->Integral();

    for (int ibg=0; ibg < nbg; ++ibg) {
      hmep_bg_b[ibg][iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",ibef))->Clone(Form("hmep_bg_b_p%d",iplab));
      hmep_bg_a[ibg][iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",iaft))->Clone(Form("hmep_bg_a_p%d",iplab));
      hmep_valid_bg_b[ibg][iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_valid_%d",ibef))->Clone(Form("hmep_valid_bg_b_p%d",iplab));
      hmep_valid_bg_a[ibg][iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_valid_%d",iaft))->Clone(Form("hmep_valid_bg_a_p%d",iplab));

      set_style(hmep_bg_b[ibg][iplab],col_bg[ibg],2);
      set_style(hmep_bg_a[ibg][iplab],col_bg[ibg],2);
      set_style(hmep_valid_bg_b[ibg][iplab],col_bg[ibg],2);
      set_style(hmep_valid_bg_a[ibg][iplab],col_bg[ibg],2);

      if (ibg==0){
      	hmep_bg_b[ibg][iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);
      	hmep_bg_a[ibg][iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);
      	hmep_valid_bg_b[ibg][iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);
      	hmep_valid_bg_a[ibg][iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);
      }

      hmep_fg_b[iplab]->Add(hmep_bg_b[ibg][iplab]);
      hmep_fg_a[iplab]->Add(hmep_bg_a[ibg][iplab]);
      hmep_valid_fg_b[iplab]->Add(hmep_valid_bg_b[ibg][iplab]);
      hmep_valid_fg_a[iplab]->Add(hmep_valid_bg_a[ibg][iplab]);

      if (ibg==0) {
      	hmep_bg_tot_b[iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",ibef))->Clone(Form("hmep_bg_tot_b_p%d",iplab));
      	hmep_bg_tot_a[iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",iaft))->Clone(Form("hmep_bg_tot_b_p%d",iplab));
      	hmep_valid_bg_tot_b[iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_valid_%d",ibef))->Clone(Form("hmep_bg_tot_b_p%d",iplab));
      	hmep_valid_bg_tot_a[iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_valid_%d",iaft))->Clone(Form("hmep_bg_tot_a_p%d",iplab));

      	set_style(hmep_bg_tot_b[iplab],col_bg[ibg],2);
      	set_style(hmep_bg_tot_a[iplab],col_bg[ibg],2);
      	set_style(hmep_valid_bg_tot_b[iplab],col_bg[ibg],2);
      	set_style(hmep_valid_bg_tot_a[iplab],col_bg[ibg],2);

      	hmep_bg_tot_b[iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);
      	hmep_bg_tot_a[iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);
      	hmep_valid_bg_tot_b[iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);
      	hmep_valid_bg_tot_a[iplab]->Scale(pi0pi0jpsi_re_scale[iplab]);

      } else {
      	hmep_bg_tot_b[iplab]->Add(hmep_bg_b[ibg][iplab]);;
      	hmep_bg_tot_a[iplab]->Add(hmep_bg_a[ibg][iplab]);;
      	hmep_valid_bg_tot_b[iplab]->Add(hmep_valid_bg_b[ibg][iplab]);;
      	hmep_valid_bg_tot_a[iplab]->Add(hmep_valid_bg_a[ibg][iplab]);;
      }

      if (iplab==0) {
	tl->AddEntry(hmep_bg_a[ibg][iplab],leg_bg[ibg],"l");
	tl2->AddEntry(hmep_bg_a[ibg][iplab],leg_bg[ibg],"l");
	tl4->AddEntry(hmep_bg_a[ibg][iplab],leg_bg[ibg],"l");
	if (ibg==3) {
	  tl3->AddEntry(hmep_bg_a[ibg][iplab],leg_bg[ibg],"l");
	}
      }

      hpi0jpsi_chi24c[ibg+1][iplab] = (TH1F*) fbg[ibg][iplab]->Get("hpi0jpsi_chi24c")->Clone(Form("hpi0jpsi_chi24c_t%d_p%d",ibg+1,iplab));
      hpi0jpsi_chi24cc[ibg+1][iplab] = (TH1F*) fbg[ibg][iplab]->Get("hpi0jpsi_chi24c_c")->Clone(Form("hpi0jpsi_chi24cc_t%d_p%d",ibg+1,iplab));
      set_style(hpi0jpsi_chi24c[ibg+1][iplab],col_bg[ibg],4);
      set_style(hpi0jpsi_chi24cc[ibg+1][iplab],col_bg[ibg]);

      if (ibg==0) {
	hpi0pi0jpsi_chi24c[ibg+1][iplab] = (TH1F*) fbg[ibg][iplab]->Get("hpi0pi0jpsi_chi24c")->Clone(Form("hpi0pi0jpsi_chi24c_t%d_p%d",ibg+1,iplab));
	hpi0pi0jpsi_chi24cc[ibg+1][iplab] = (TH1F*) fbg[ibg][iplab]->Get("hpi0pi0jpsi_chi24c_c")->Clone(Form("hpi0pi0jpsi_chi24cc_t%d_p%d",ibg+1,iplab));
	set_style(hpi0pi0jpsi_chi24c[ibg+1][iplab],col_bg[ibg],4);
	set_style(hpi0pi0jpsi_chi24cc[ibg+1][iplab],col_bg[ibg]);
	hpi0vs2pi0_chi24c[ibg+1][iplab] = (TH2F*) fbg[ibg][iplab]->Get("hpi0vs2pi0_chi24c")->Clone(Form("hpi0vs2pi0_chi24c_t%d_p%d",ibg+1,iplab));
	hpi0vs2pi0_chi24cc[ibg+1][iplab] = (TH2F*) fbg[ibg][iplab]->Get("hpi0vs2pi0_chi24c_c")->Clone(Form("hpi0vs2pi0_chi24cc_t%d_p%d",ibg+1,iplab));
	set_style(hpi0vs2pi0_chi24c[ibg+1][iplab]);
	set_style(hpi0vs2pi0_chi24cc[ibg+1][iplab]);
      }

    }

    int mmin = hmep_sg_b[iplab]->GetXaxis()->FindBin(2.8);
    int mmax = hmep_sg_b[iplab]->GetXaxis()->FindBin(3.3);
    //cout << "==================================" << endl;
    //cout << "plab= " << plab[iplab] << " step= " << iaft << "  CountSig(3sig-jpsi) = " << hmep_sg_a[iplab]->Integral(mmin,mmax) << endl;
    //cout << "EffSig = " << Form("%4.1f\\%%", 100*hmep_sg_a[iplab]->Integral(mmin,mmax) / nevt_sig_sim_valid[iplab]) << endl;
    //cout << "EffPi0pi0jpsi " << Form("%3.1e\\%%", 100*hmep_bg_a[0][iplab]->Integral(mmin,mmax)/ nevt_pi0pi0jpsi_sim_valid[iplab]) << endl;
    //cout << "Effpi0pipm " << Form("%3.1e\\%%", 100*hmep_bg_a[1][iplab]->Integral(mmin,mmax)/ nevt_pi0pipm_sim_valid[iplab]) << endl;
    ////cout << "EffPi0pi0jpsi " << hmep_bg_a[0][iplab]->Integral(mmin,mmax)/ nevt_pi0pi0jpsi_sim_valid[iplab] << endl;
    ////cout << "Effpi0pipm " << hmep_bg_a[1][iplab]->Integral(mmin,mmax)/ nevt_pi0pipm_sim_valid[iplab] << endl;
    //cout << "ContPi0pi0jpsi " << Form("%3.1f\\%%", 100*hmep_bg_a[0][iplab]->Integral(mmin,mmax)/ hmep_sg_a[iplab]->Integral(mmin,mmax)) << endl;
    //cout << "Contpi0pipm " << Form("%3.1f\\%%", 100*hmep_bg_a[1][iplab]->Integral(mmin,mmax)/ hmep_sg_a[iplab]->Integral(mmin,mmax)) << endl;

    cout << "===========" << endl;
    cout << "& " << Form("%4.1f\\%%", 100*hmep_valid_sg_a[iplab]->Integral(mmin,mmax) / nevt_sig_sim_valid[iplab]) << endl;

    //cout << "& " << Form("%3.1e", nevt_pi0pi0jpsi_sim_valid[iplab] / hmep_bg_a[0][iplab]->Integral(mmin,mmax) ) << endl;
    //cout << "& " << Form("%3.1e",  pi0pipm_frac_valid[iplab] * nevt_pi0pipm_sim_valid[iplab] / hmep_bg_a[1][iplab]->Integral(mmin,mmax) ) << endl;

    //cout << "& " << Form("%3.1e", nevt_xsect[3][iplab] / hmep_bg_a[0][iplab]->Integral(mmin,mmax) ) << endl;
    //cout << "& " << Form("%3.1e",  nevt_xsect[0][iplab] / hmep_bg_a[1][iplab]->Integral(mmin,mmax) ) << endl;

    //cout << "& " << Form("%3.1e\\%%", 100*hmep_bg_a[0][iplab]->Integral(mmin,mmax)/ nevt_pi0pi0jpsi_sim_valid[iplab]) << endl;
    //cout << "& " << Form("%3.1e\\%%", 100*hmep_bg_a[1][iplab]->Integral(mmin,mmax)/ nevt_pi0pipm_sim_valid[iplab]) << endl;
    //cout << "EffPi0pi0jpsi " << hmep_bg_a[0][iplab]->Integral(mmin,mmax)/ nevt_pi0pi0jpsi_sim_valid[iplab] << endl;
    //cout << "Effpi0pipm " << hmep_bg_a[1][iplab]->Integral(mmin,mmax)/ nevt_pi0pipm_sim_valid[iplab] << endl;
    cout << "& " << Form("%3.1f\\%%", 100*hmep_valid_bg_a[0][iplab]->Integral(mmin,mmax)/ hmep_valid_sg_a[iplab]->Integral(mmin,mmax)) << endl;
    cout << "& " << Form("%3.1f\\%%", 100*hmep_valid_bg_a[1][iplab]->Integral(mmin,mmax)/ hmep_valid_sg_a[iplab]->Integral(mmin,mmax)) << endl;

    cout << leg_bg[2] << " & " << Form("%10.8f\\%%", 100*hmep_valid_bg_a[2][iplab]->Integral(mmin,mmax)/ hmep_valid_sg_a[iplab]->Integral(mmin,mmax)) << endl;
    cout << leg_bg[3] << " & " << Form("%10.8f\\%%", 100*hmep_valid_bg_a[3][iplab]->Integral(mmin,mmax)/ hmep_valid_sg_a[iplab]->Integral(mmin,mmax)) << endl;

    cout << "& " << Form("%4.1f\\%%", purity[iaft-5][iplab]) << endl;

    double bg_tot = hmep_valid_bg_tot_a[iplab]->Integral(mmin,mmax);
    cout << "bgtot = " << bg_tot << endl;
    cout << "bg sum = " << hmep_valid_sg_a[iplab]->Integral(mmin,mmax) + hmep_valid_sg_a[iplab]->Integral(mmin,mmax) << endl;
    cout << "& " << Form("%3.1f",  hmep_valid_sg_a[iplab]->Integral(mmin,mmax)/bg_tot) << endl;
    //double bg_tot = hmep_bg_tot_a[iplab]->Integral(mmin,mmax);
    //cout << "& " << Form("%3.1f",  hmep_sg_a[iplab]->Integral(mmin,mmax)/bg_tot) << endl;

    //cout << "COUNTS:" << endl;
    //cout << "BEF(ip="<<iplab<<"): Sig (3sig-jpsi) = " << hmep_sg_b[iplab]->Integral(mmin,mmax) << endl;
    //cout << "BEF(ip="<<iplab<<"): pi0pi0jpsi (3sig-jpsi) = "  << hmep_bg_b[3][iplab]->Integral(mmin,mmax) << endl;
    //cout << "BEF(ip="<<iplab<<"): pi0pipm (3sig-jpsi) = "  << hmep_bg_b[0][iplab]->Integral(mmin,mmax) << endl;
    //cout << "BEF(ip="<<iplab<<"): pi02pipm (3sig-jpsi) = "  << hmep_bg_b[1][iplab]->Integral(mmin,mmax) << endl;
    //cout << "BEF(ip="<<iplab<<"): pi0pipm2 (3sig-jpsi) = "  << hmep_bg_b[2][iplab]->Integral(mmin,mmax) << endl;
    //cout << "BEF(ip="<<iplab<<"): BgTot (3sig-jpsi) = "  << hmep_bg_tot_b[iplab]->Integral(mmin,mmax) << endl;
    //
    //cout << "AFT(ip="<<iplab<<"): Sig (3sig-jpsi) = " << hmep_sg_a[iplab]->Integral(mmin,mmax) << endl;
    //cout << "AFT(ip="<<iplab<<"): pi0pi0jspi (3sig-jpsi) = "  << hmep_bg_a[3][iplab]->Integral(mmin,mmax) << endl;
    //cout << "AFT(ip="<<iplab<<"): pi0jspi (3sig-jpsi) = "  << hmep_bg_a[0][iplab]->Integral(mmin,mmax) << endl;
    //cout << "AFT(ip="<<iplab<<"): pi02pipm (3sig-jpsi) = "  << hmep_bg_a[1][iplab]->Integral(mmin,mmax) << endl;
    //cout << "AFT(ip="<<iplab<<"): pi0pipm2 (3sig-jpsi) = "  << hmep_bg_a[2][iplab]->Integral(mmin,mmax) << endl;
    //cout << "AFT(ip="<<iplab<<"): BgTot (3sig-jpsi) = "  << hmep_bg_tot_a[iplab]->Integral(mmin,mmax) << endl;
    //
    //cout << "S/B:" << endl;
    //cout << "BEF(ip="<<iplab<<"): Sig/pi0pi0jpsi = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[3][iplab]) << endl;
    //cout << "BEF(ip="<<iplab<<"): Sig/pi0pipm = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[0][iplab]) << endl;
    //cout << "BEF(ip="<<iplab<<"): Sig/pi02pipm = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[1][iplab]) << endl;
    //cout << "BEF(ip="<<iplab<<"): Sig/pi0pipm2 = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[2][iplab]) << endl;
    //cout << "BEF(ip="<<iplab<<"): Sig/Bgtot = " << stob_r( hmep_sg_b[iplab], hmep_bg_tot_b[iplab]) << endl;
    //
    //cout << "AFT(ip="<<iplab<<"): Sig/pi0pi0jpsi = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[3][iplab]) << endl;
    //cout << "AFT(ip="<<iplab<<"): Sig/pi0pipm = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[0][iplab]) << endl;
    //cout << "AFT(ip="<<iplab<<"): Sig/pi02pipm = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[1][iplab]) << endl;
    //cout << "AFT(ip="<<iplab<<"): Sig/pi0pipm2 = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[2][iplab]) << endl;
    //cout << "AFT(ip="<<iplab<<"): Sig/Bgtot = " << stob_r( hmep_sg_a[iplab], hmep_bg_tot_a[iplab]) << endl;
    //
    //
    //cout << "Rejection:" << endl;
    ////
    ////double nevt_xsect[5][3] = {{4.0e11, 1e11, 2e10}, {32780.0,50142.0,51860.0}, {1.15e12, 3.15e11, 6.84e10}, {3.19e12, 1.14e12, 2.92e11}, {94243.0, 157947.3, 177361.2}};
    //cout << "Rejection:" << endl;
    //cout << "BEF(ip="<<iplab<<"): pi0pi0jpsi = " << nevt_xsect[4][iplab]/count_jpsim( hmep_bg_b[3][iplab]) << endl;
    //cout << "BEF(ip="<<iplab<<"): pi0pipm = " <<  nevt_xsect[0][iplab]/count_jpsim( hmep_bg_b[0][iplab]) << endl;
    //cout << "BEF(ip="<<iplab<<"): pi02pipm = " << nevt_xsect[2][iplab]/count_jpsim( hmep_bg_b[1][iplab]) << endl;
    //cout << "BEF(ip="<<iplab<<"): pi0pipm2 = " << nevt_xsect[3][iplab]/count_jpsim( hmep_bg_b[2][iplab]) << endl;
    //
    //cout << "AFT(ip="<<iplab<<"): pi0pi0jpsi = " << nevt_xsect[4][iplab]/count_jpsim( hmep_bg_a[3][iplab]) << endl;
    //cout << "AFT(ip="<<iplab<<"): pi0pipm = " <<  nevt_xsect[0][iplab]/count_jpsim( hmep_bg_a[0][iplab]) << endl;
    //cout << "AFT(ip="<<iplab<<"): pi02pipm = " << nevt_xsect[2][iplab]/count_jpsim( hmep_bg_a[1][iplab]) << endl;
    //cout << "AFT(ip="<<iplab<<"): pi0pipm2 = " << nevt_xsect[3][iplab]/count_jpsim( hmep_bg_a[2][iplab]) << endl;

  }

  TCanvas *tc_bkf = new TCanvas("invm_bgsrc_bef_kinfit", "tc_bkf", 1600, 700);
  tc_bkf->Divide(3,1);
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_bkf->cd(1+iplab);

    tc_bkf->GetPad(1+iplab)->SetTopMargin(0.08);
    tc_bkf->GetPad(1+iplab)->SetLeftMargin(0.15);
    tc_bkf->GetPad(1+iplab)->SetRightMargin(0.001);
    tc_bkf->GetPad(1+iplab)->SetBottomMargin(0.13);

    hmep_fg_b[iplab]->GetYaxis()->SetTitleOffset(1.2);
    hmep_fg_b[iplab]->GetYaxis()->SetLabelOffset(0.0);
    hmep_fg_b[iplab]->SetMaximum(iplab==2?4e4:2e4);
    hmep_fg_b[iplab]->GetXaxis()->SetNdivisions(1005,false);

    hmep_fg_b[iplab]->SetLineWidth(1);
    hmep_sg_b[iplab]->SetLineWidth(1);
    for (int ibg=0; ibg < nbg; ++ibg) hmep_bg_b[ibg][iplab]->SetLineWidth(1);

    hmep_fg_b[iplab]->Draw();
    if (iplab==2) {
      tl4->Draw();
      hmep_fg_b[iplab]->Draw("same");
    }
    hmep_sg_b[iplab]->Draw("same");
    for (int ibg=0; ibg < nbg; ++ibg) hmep_bg_b[ibg][iplab]->Draw("same");
    gPad->SetLogy();

    TLatex *tlat = new TLatex();
    //tlat->SetNDC(true);
    tlat->SetTextSize(0.05);
    //tlat->DrawLatexNDC(0.65,0.88,Form("p_{#bar{p}} = %4.1f GeV/c", plab[iplab]));
    tlat->DrawLatexNDC(iplab==2?0.65:0.55,0.87,Form("p_{#bar{p}} = %4.1f GeV/c", plab[iplab]));
  }
  //tc_bkf->Print("invm_bgsrc_bef_kinfit.pdf");

  TCanvas *tc_bkf_valid = new TCanvas("invm_valid_bgsrc_bef_kinfit", "tc_bkf_valid", 1600, 700);
  tc_bkf_valid->Divide(3,1);
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_bkf_valid->cd(1+iplab);

    tc_bkf_valid->GetPad(1+iplab)->SetTopMargin(0.08);
    tc_bkf_valid->GetPad(1+iplab)->SetLeftMargin(0.15);
    tc_bkf_valid->GetPad(1+iplab)->SetRightMargin(0.001);
    tc_bkf_valid->GetPad(1+iplab)->SetBottomMargin(0.13);

    hmep_valid_fg_b[iplab]->GetYaxis()->SetTitleOffset(1.2);
    hmep_valid_fg_b[iplab]->GetYaxis()->SetLabelOffset(0.0);
    hmep_valid_fg_b[iplab]->SetMaximum(2e4);
    hmep_valid_fg_b[iplab]->GetXaxis()->SetNdivisions(1005,false);

    hmep_valid_fg_b[iplab]->SetLineWidth(1);
    hmep_valid_sg_b[iplab]->SetLineWidth(1);
    for (int ibg=0; ibg < nbg; ++ibg) hmep_valid_bg_b[ibg][iplab]->SetLineWidth(1);

    hmep_valid_fg_b[iplab]->Draw();
    if (iplab==2) {
      tl4->Draw();
      hmep_valid_fg_b[iplab]->Draw("same");
    }
    hmep_valid_sg_b[iplab]->Draw("same");
    for (int ibg=0; ibg < nbg; ++ibg) hmep_valid_bg_b[ibg][iplab]->Draw("same");
    gPad->SetLogy();

    TLatex *tlat = new TLatex();
    //tlat->SetNDC(true);
    tlat->SetTextSize(0.05);
    //tlat->DrawLatexNDC(0.65,0.88,Form("p_{#bar{p}} = %4.1f GeV/c", plab[iplab]));
    tlat->DrawLatexNDC(iplab==2?0.65:0.55,0.87,Form("p_{#bar{p}} = %4.1f GeV/c", plab[iplab]));
  }
  //tc_bkf_valid->Print("invm_valid_bgsrc_bef_kinfit.pdf");


  TCanvas *tc_akf = new TCanvas("invm_bgsrc_aft_kinfit", "tc_akf", 1600, 700);
  tc_akf->Divide(3,1);
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_akf->cd(1+iplab);
    tc_akf->GetPad(1+iplab)->SetRightMargin(0.001);
    tc_akf->GetPad(1+iplab)->SetLeftMargin(0.1);
    hmep_fg_a[iplab]->Draw();
    hmep_sg_a[iplab]->Draw("same");
    for (int ibg=0; ibg < nbg; ++ibg) hmep_bg_a[ibg][iplab]->Draw("same");
    //gPad->SetLogy();
    if (iplab==2) tl4->Draw();
  }
  //tc_akf->Print(iaft==6?"invm_bgsrc_aft_chi2sig_cut.pdf":"invm_bgsrc_aft_kincuts.pdf");
  //tc_akf->Print(iaft==6?"invm_bgsrc_aft_chi2sig_cut.png":"invm_bgsrc_aft_kincuts.png");

  TCanvas *tc[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc[iplab] = new TCanvas(Form("tc%d",iplab),Form("tc%d",iplab),1300,700);
    tc[iplab]->Divide(2,1);
    tc[iplab]->cd(1);
    hmep_fg_b[iplab]->Draw();
    hmep_sg_b[iplab]->Draw("same");
    //tl->Draw();
    //gPad->SetLogy();
    for (int ibg=0; ibg < nbg; ++ibg) hmep_bg_b[ibg][iplab]->Draw("same");
    tc[iplab]->cd(2);
    hmep_fg_a[iplab]->Draw();
    hmep_sg_a[iplab]->Draw("same");
    tl2->Draw();
    for (int ibg=0; ibg < nbg; ++ibg) hmep_bg_a[ibg][iplab]->Draw("same");
    gPad->SetLogy();
    //tc[iplab]->Print(Form("%s/marc_meeting/minv_sig_and_bghyp_pessimistic_norm_p%d.eps",bdir,iplab));
  }

  TCanvas *tc_chi2[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_chi2[iplab] = new TCanvas(Form("tc_chi2_p%d",iplab),Form("tc_chi2_p%d",iplab));

    tc_chi2[iplab]->cd();
    tc_chi2[iplab]->GetPad(0)->SetTopMargin(0.08);
    tc_chi2[iplab]->GetPad(0)->SetBottomMargin(0.14);
    tc_chi2[iplab]->GetPad(0)->SetLeftMargin(0.13);
    tc_chi2[iplab]->GetPad(0)->SetRightMargin(0.1);
    TGaxis::SetMaxDigits(3);

    for (int ii=0; ii < nbg+1; ++ii) {
      //if (iplab==0)
      //	hpi0jpsi_chi24c[ii][iplab]->SetTitle(Form("#chi^{2} of 2#gamma e^{+}e^{-} hyp. 4C kinematic fit (p_{#bar{p}} = %5.4g GeV/c); #chi^{2}; arb. unit",plab[iplab]));
      //else
      //	hpi0jpsi_chi24c[ii][iplab]->SetTitle(Form("#chi^{2} of 2#gamma e^{+}e^{-} hyp. 4C kinematic fit (p_{#bar{p}} = %4.1f GeV/c); #chi^{2}; arb. unit",plab[iplab]));

      hpi0jpsi_chi24c[ii][iplab]->SetTitle("; #chi^{2}; arb. unit");
      hpi0jpsi_chi24c[ii][iplab]->GetXaxis()->SetNdivisions(1005,false);
      hpi0jpsi_chi24c[ii][iplab]->DrawNormalized(ii==0?"":"same");
    }
    gPad->SetLogy();
    tl->Draw();

    tc_chi2[iplab]->cd();
    TPad *pad = new TPad("inset","inset",0.42,0.49,0.89,0.93);
    pad->SetRightMargin(0.02);
    pad->SetLeftMargin(0.18);
    pad->SetBottomMargin(0.14);
    pad->SetFillStyle(4000);
    pad->SetFillColorAlpha(0,0);
    pad->Draw();
    pad->cd();
    for (int ii=0; ii < nbg+1; ++ii) {
      hpi0jpsi_chi24cc[ii][iplab]->SetTitle("");
      hpi0jpsi_chi24cc[ii][iplab]->GetXaxis()->SetRangeUser(0,60);
      hpi0jpsi_chi24cc[ii][iplab]->SetTitle("; #chi^{2}; arb. unit");
      hpi0jpsi_chi24cc[ii][iplab]->DrawNormalized(ii==0?"":"same");
    }

    tc_chi2[iplab]->cd();
    TLatex *tlat = new TLatex();
    tlat->SetTextSize(0.04);
    tlat->DrawLatexNDC(0.65,0.76,"#chi^{2} (2#gamma e^{+}e^{-} fit)");
    tlat->DrawLatexNDC(0.65,0.7,Form("p_{#bar{p}} = %4.1f GeV/c",plab[iplab]));

    //tc_chi2[iplab]->Print(Form("%s/marc_meeting/chi2_sighyp_p%d.eps",bdir,iplab));
    tc_chi2[iplab]->Print(Form("kinfit_4c_chi2_dists_p%d.pdf",iplab));
  }

  return;

  TCanvas *tc_chi2bg[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_chi2bg[iplab] = new TCanvas(Form("tc_chi2bg_p%d",iplab),Form("tc_chi2bg_p%d",iplab));
    tc_chi2bg[iplab]->cd();
    if (iplab==0)
      hpi0pi0jpsi_chi24c[1][iplab]->SetTitle(Form("#chi^{2} of 4#gamma e^{+}e^{-} 4C kinematic fit (p_{#bar{p}} = %5.4g GeV/c); #chi^{2}; arb. unit",plab[iplab]));
    else
      hpi0pi0jpsi_chi24c[1][iplab]->SetTitle(Form("#chi^{2} of 4#gamma e^{+}e^{-} 4C kinematic fit (p_{#bar{p}} = %4.1f GeV/c); #chi^{2}; arb. unit",plab[iplab]));
    hpi0pi0jpsi_chi24c[1][iplab]->GetXaxis()->SetNdivisions(505);
    hpi0pi0jpsi_chi24c[1][iplab]->DrawNormalized();
    hpi0pi0jpsi_chi24c[0][iplab]->DrawNormalized("same");
    tl3->Draw();
  }

  TCanvas *tc_chi2bgvssig[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_chi2bgvssig[iplab] = new TCanvas(Form("tc_chi2bgvssig_p%d",iplab),Form("tc_chi2bgvssig_p%d",iplab),1300,700);
    tc_chi2bgvssig[iplab]->Divide(2,1);

    for (int ipad=0; ipad < 2; ++ipad) {
      tc_chi2bgvssig[iplab]->GetPad(1+ipad)->SetTopMargin(0.08);
      tc_chi2bgvssig[iplab]->GetPad(1+ipad)->SetBottomMargin(0.14);
      tc_chi2bgvssig[iplab]->GetPad(1+ipad)->SetLeftMargin(0.13);
      tc_chi2bgvssig[iplab]->GetPad(1+ipad)->SetRightMargin(0.1);
    }

    hpi0vs2pi0_chi24c[1][iplab]->SetTitle(Form("; #chi^{2} (2#gamma e^{+}e^{-} fit); #chi^{2} (4#gamma e^{+}e^{-} fit)"));
    hpi0vs2pi0_chi24c[0][iplab]->SetTitle(Form("; #chi^{2} (2#gamma e^{+}e^{-} fit); #chi^{2} (4#gamma e^{+}e^{-} fit)"));

    //if (iplab==0) {
    //  hpi0vs2pi0_chi24c[1][iplab]->SetTitle(Form("#chi^{2} 2#gamma e^{+}e^{-} vs 4#gamma e^{+}e^{-} 4C fit. (2(#pi^{0})J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %5.4g GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
    //  hpi0vs2pi0_chi24c[0][iplab]->SetTitle(Form("#chi^{2} 2#gamma e^{+}e^{-} vs 4#gamma e^{+}e^{-} 4C fit. (#pi^{0}J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %5.4g GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
    //} else {
    //  hpi0vs2pi0_chi24c[1][iplab]->SetTitle(Form("#chi^{2} 2#gamma e^{+}e^{-} vs 4#gamma e^{+}e^{-} 4C fit. (2(#pi^{0})J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %4.1f GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
    //  hpi0vs2pi0_chi24c[0][iplab]->SetTitle(Form("#chi^{2} 2#gamma e^{+}e^{-} vs 4#gamma e^{+}e^{-} 4C fit. (#pi^{0}J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %4.1f GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
    //}
    TGaxis::SetMaxDigits(3);
    hpi0vs2pi0_chi24c[1][iplab]->GetXaxis()->SetNdivisions(505);
    hpi0vs2pi0_chi24c[1][iplab]->GetYaxis()->SetNdivisions(504);
    hpi0vs2pi0_chi24c[0][iplab]->GetXaxis()->SetNdivisions(505);
    hpi0vs2pi0_chi24c[0][iplab]->GetYaxis()->SetNdivisions(504);

    TLatex *tlat = new TLatex();
    tlat->SetTextSize(0.04);

    tc_chi2bgvssig[iplab]->cd(1);
    hpi0vs2pi0_chi24c[0][iplab]->Draw("box");
    tlat->DrawLatexNDC(0.35,0.86,Form("#pi^{0}J/#psi(#rightarrow e^{+}e^{-})"));
    tlat->DrawLatexNDC(0.35,0.8,Form("p_{#bar{p}} = %3.1f GeV/c", plab[iplab]));

    tc_chi2bgvssig[iplab]->cd(2);
    hpi0vs2pi0_chi24c[1][iplab]->Draw("col");
    tlat->DrawLatexNDC(0.45,0.86,Form("#pi^{0}#pi^{0}J/#psi(#rightarrow e^{+}e^{-})"));
    tlat->DrawLatexNDC(0.45,0.8,Form("p_{#bar{p}} = %3.1f GeV/c", plab[iplab]));

    //tc_chi2bgvssig[iplab]->Print(Form("%s/marc_meeting/chi2_sig_vs_bg_hyp_p%d.eps",bdir,iplab));
    //tc_chi2bgvssig[iplab]->Print(Form("chi2_sig_vs_bg_hyp_p%d.pdf",iplab));
  }

}

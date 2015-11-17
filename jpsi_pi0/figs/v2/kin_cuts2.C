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

void kin_cuts2(int pass=14) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  gROOT->LoadMacro(Form("%s/figs/v2/ananote.C",bdir));

  const int nbg = 4;
  int col_bg[nbg] = {41, 46, 38, 30};
  TFile *fsg[nplab], *fbg[nbg][nplab];
  TH1F *hmep_fg_b[nplab], *hmep_sg_b[nplab], *hmep_bg_b[nbg][nplab], *hmep_bg_tot_b[nplab];
  TH1F *hmep_fg_a[nplab], *hmep_sg_a[nplab], *hmep_bg_a[nbg][nplab], *hmep_bg_tot_a[nplab];

  TH1F *hpi0jpsi_chi24c[nbg+1][nplab], *hpi0jpsi_chi24cc[nbg+1][nplab];
  TH1F *hpi0pi0jpsi_chi24c[nbg+1][nplab], *hpi0pi0jpsi_chi24cc[nbg+1][nplab];
  TH2F *hpi0vs2pi0_chi24c[nbg+1][nplab], *hpi0vs2pi0_chi24cc[nbg+1][nplab];

  // temporary, until fixed in effing macro
  double nevt_sim[5][3] = {{814794.0,888292.0,898721.0}, {32780.0,50142.0,51860.0}, {214780.0,174864.0,160099.0}, {570751.0,609044.0,527506.0}, {200000.0,200000.0,200000.0}};
  double nevt_xsect[5][3] = {{4.0e11, 1e11, 2e10}, {32780.0,50142.0,51860.0}, {1.15e12, 3.15e11, 6.84e10}, {3.19e12, 1.14e12, 2.92e11}, {94243.0, 157947.3, 177361.2}};
  double pi0pi0jpsi_scale[3] = {1.0};
  for (int iplab=0; iplab < nplab; ++iplab) {
    nevt_xsect[4][iplab] = nevt_xsect[1][iplab]*(nevt_xsect[2][iplab]/nevt_xsect[0][iplab]);
    pi0pi0jpsi_scale[iplab] = nevt_xsect[4][iplab]/nevt_sim[4][iplab];
    cout << "nxsect ip=" << iplab << " = " << nevt_xsect[4][iplab] << endl;
    cout << "nsim ip=" << iplab << " = " << nevt_sim[4][iplab] << endl;
    cout << "scale ip=" << iplab << " = " << pi0pi0jpsi_scale[iplab] << endl;
  }

  //double tmp_scale[3] = {0.0};
  //for (int ii=0; ii < 3; ++ii) {
  //  tmp_scale[ii] = nevt_xsect[0][ii]/nevt_xsect[2][ii];
  //  cout << "tmp_scale" << ii << " = " << tmp_scale[ii] << endl;
  //}
  //return;
  //double tmp_scale[3] = {32780.0/200000.0,50142.0/200000.0,51860.0/200000.0};
  //double tmp_scale[3]= {94243.0/200000.0, 157947.3/200000.0, 177361.2/200000.0};

  //const char *leg_bg[nbg] = {"#pi^{0}#pi^{+}#pi^{-}","#pi^{0}#pi^{0}#pi^{+}#pi^{-}", "#pi^{0}#pi^{+}#pi^{-}#pi^{+}#pi^{-}", "#pi^{0}J/#psi(#pi^{+}#pi^{-})"};
  const char *leg_bg[nbg] = {"#pi^{0}#pi^{+}#pi^{-}","2(#pi^{0})#pi^{+}#pi^{-}", "#pi^{0}2(#pi^{+}#pi^{-})", "2(#pi^{0})J/#psi(#rightarrow e^{+}e^{-})"};
  //const double scale[nplab][nbg] = {{1.,1.,1.,0.037817}, {1.,1.,1.,0.037817}, {1.,1.,1.,0.037817}}; //  (30pb for pi0pi0jpsi )
  //const double scale[nplab][nbg] = {{1.,1.,1.,0.64}, {1.,1.,1.,0.64}, {1.,1.,1.,0.64}}; // (30pb for pi0pi0jpsi->e+e-4gamma)


  TLegend *tl = new TLegend(0.17,0.5,0.45,0.85);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);

  TLegend *tl2 = new TLegend(0.17,0.5,0.65,0.85);
  tl2->SetBorderSize(0);
  tl2->SetFillStyle(0);

  TLegend *tl3 = new TLegend(0.17,0.5,0.65,0.85);
  tl3->SetBorderSize(0);
  tl3->SetFillStyle(0);

  TLegend *tl4 = new TLegend(0.11,0.5,0.65,0.85);
  tl4->SetTextSize(0.06);
  tl4->SetBorderSize(0);
  tl4->SetFillStyle(0);

  int ibef = 5;
  int iaft = 8;

  for (int iplab=0; iplab < nplab; ++iplab) {

    fsg[iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[0][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[1][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi02pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[2][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pipm2_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[3][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));

    hmep_sg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",ibef)) ->Clone(Form("hmep_sg_b_p%d",iplab));
    hmep_sg_a[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",iaft)) ->Clone(Form("hmep_sg_a_p%d",iplab));
    hmep_fg_b[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",ibef)) ->Clone(Form("hmep_fg_b_p%d",iplab));
    hmep_fg_a[iplab] = (TH1F*) fsg[iplab]->Get(Form("hmep_%d",iaft)) ->Clone(Form("hmep_fg_a_p%d",iplab));

    set_style(hmep_fg_b[iplab],1);
    set_style(hmep_fg_a[iplab],1);
    set_style(hmep_sg_b[iplab],2);
    set_style(hmep_sg_a[iplab],2);

    if (iplab==0) {
      hmep_fg_b[iplab]->SetTitle(Form("Before cut on 4C fit #chi^{2} (p_{#bar{p}} = %5.4g GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
      hmep_fg_a[iplab]->SetTitle(Form("After cuts on 4C fit #chi^{2} (p_{#bar{p}} = %5.4g GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
    } else {
      hmep_fg_b[iplab]->SetTitle(Form("Before cut on 4C fit #chi^{2} (p_{#bar{p}} = %4.1f GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
      hmep_fg_a[iplab]->SetTitle(Form("After cuts on 4C fit #chi^{2} (p_{#bar{p}} = %4.1f GeV/c); M_{ee}[GeV/c^{2}]", plab[iplab]));
    }

    if (iplab==0) {
      //tl->AddEntry(hmep_fg_a[iplab],"FG+BG","l");
      tl2->AddEntry(hmep_fg_a[iplab],"FG+BG","l");
      tl4->AddEntry(hmep_fg_a[iplab],"FG+BG","l");
      tl->AddEntry(hmep_sg_a[iplab],"#pi^{0}J/#psi(#rightarrow e^{+}e^{-})","l");
      tl2->AddEntry(hmep_sg_a[iplab],"#pi^{0}J/#psi(#rightarrow e^{+}e^{-})","l");
      tl3->AddEntry(hmep_sg_a[iplab],"#pi^{0}J/#psi(#rightarrow e^{+}e^{-})","l");
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

    for (int ibg=0; ibg < nbg; ++ibg) {
      hmep_bg_b[ibg][iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",ibef))->Clone(Form("hmep_bg_b_p%d",iplab));
      hmep_bg_a[ibg][iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",iaft))->Clone(Form("hmep_bg_a_p%d",iplab));
      set_style(hmep_bg_b[ibg][iplab],col_bg[ibg]);
      set_style(hmep_bg_a[ibg][iplab],col_bg[ibg]);

      //if (ibg==nbg-1){
      //	hmep_bg_b[ibg][iplab]->Scale(pi0pi0jpsi_scale[iplab]);
      //	hmep_bg_a[ibg][iplab]->Scale(pi0pi0jpsi_scale[iplab]);
      //}

      hmep_fg_b[iplab]->Add(hmep_bg_b[ibg][iplab]);
      hmep_fg_a[iplab]->Add(hmep_bg_a[ibg][iplab]);

      if (ibg==0) {
	hmep_bg_tot_b[iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",ibef))->Clone(Form("hmep_bg_tot_b_p%d",iplab));
	hmep_bg_tot_a[iplab] = (TH1F*) fbg[ibg][iplab]->Get(Form("hmep_%d",iaft))->Clone(Form("hmep_bg_tot_b_p%d",iplab));
      } else {
	hmep_bg_tot_b[iplab]->Add(hmep_bg_b[ibg][iplab]);;
	hmep_bg_tot_a[iplab]->Add(hmep_bg_a[ibg][iplab]);;
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

      if (ibg==3) {
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
    cout << "plab= " << plab[iplab] << endl;

    cout << "COUNTS:" << endl;
    cout << "BEF(ip="<<iplab<<"): Sig (3sig-jpsi) = " << hmep_sg_b[iplab]->Integral(mmin,mmax) << endl;
    cout << "BEF(ip="<<iplab<<"): pi0pi0jpsi (3sig-jpsi) = "  << hmep_bg_b[3][iplab]->Integral(mmin,mmax) << endl;
    cout << "BEF(ip="<<iplab<<"): pi0pipm (3sig-jpsi) = "  << hmep_bg_b[0][iplab]->Integral(mmin,mmax) << endl;
    cout << "BEF(ip="<<iplab<<"): pi02pipm (3sig-jpsi) = "  << hmep_bg_b[1][iplab]->Integral(mmin,mmax) << endl;
    cout << "BEF(ip="<<iplab<<"): pi0pipm2 (3sig-jpsi) = "  << hmep_bg_b[2][iplab]->Integral(mmin,mmax) << endl;
    cout << "BEF(ip="<<iplab<<"): BgTot (3sig-jpsi) = "  << hmep_bg_tot_b[iplab]->Integral(mmin,mmax) << endl;

    cout << "AFT(ip="<<iplab<<"): Sig (3sig-jpsi) = " << hmep_sg_a[iplab]->Integral(mmin,mmax) << endl;
    cout << "AFT(ip="<<iplab<<"): pi0pi0jspi (3sig-jpsi) = "  << hmep_bg_a[3][iplab]->Integral(mmin,mmax) << endl;
    cout << "AFT(ip="<<iplab<<"): pi0jspi (3sig-jpsi) = "  << hmep_bg_a[0][iplab]->Integral(mmin,mmax) << endl;
    cout << "AFT(ip="<<iplab<<"): pi02pipm (3sig-jpsi) = "  << hmep_bg_a[1][iplab]->Integral(mmin,mmax) << endl;
    cout << "AFT(ip="<<iplab<<"): pi0pipm2 (3sig-jpsi) = "  << hmep_bg_a[2][iplab]->Integral(mmin,mmax) << endl;
    cout << "AFT(ip="<<iplab<<"): BgTot (3sig-jpsi) = "  << hmep_bg_tot_a[iplab]->Integral(mmin,mmax) << endl;

    cout << "S/B:" << endl;
    cout << "BEF(ip="<<iplab<<"): Sig/pi0pi0jpsi = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[3][iplab]) << endl;
    cout << "BEF(ip="<<iplab<<"): Sig/pi0pipm = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[0][iplab]) << endl;
    cout << "BEF(ip="<<iplab<<"): Sig/pi02pipm = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[1][iplab]) << endl;
    cout << "BEF(ip="<<iplab<<"): Sig/pi0pipm2 = " << stob_r( hmep_sg_b[iplab], hmep_bg_b[2][iplab]) << endl;
    cout << "BEF(ip="<<iplab<<"): Sig/Bgtot = " << stob_r( hmep_sg_b[iplab], hmep_bg_tot_b[iplab]) << endl;

    cout << "AFT(ip="<<iplab<<"): Sig/pi0pi0jpsi = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[3][iplab]) << endl;
    cout << "AFT(ip="<<iplab<<"): Sig/pi0pipm = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[0][iplab]) << endl;
    cout << "AFT(ip="<<iplab<<"): Sig/pi02pipm = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[1][iplab]) << endl;
    cout << "AFT(ip="<<iplab<<"): Sig/pi0pipm2 = " << stob_r( hmep_sg_a[iplab], hmep_bg_a[2][iplab]) << endl;
    cout << "AFT(ip="<<iplab<<"): Sig/Bgtot = " << stob_r( hmep_sg_a[iplab], hmep_bg_tot_a[iplab]) << endl;


    cout << "Rejection:" << endl;

    cout << "Rejection:" << endl;
    cout << "BEF(ip="<<iplab<<"): pi0pi0jpsi = " << nevt_sim[3][iplab]/count_jpsim( hmep_bg_b[3][iplab]) << endl;
    cout << "BEF(ip="<<iplab<<"): pi0pipm = " <<  nevt_sim[0][iplab]/count_jpsim( hmep_bg_b[0][iplab]) << endl;
    cout << "BEF(ip="<<iplab<<"): pi02pipm = " << nevt_sim[1][iplab]/count_jpsim( hmep_bg_b[1][iplab]) << endl;
    cout << "BEF(ip="<<iplab<<"): pi0pipm2 = " << nevt_sim[2][iplab]/count_jpsim( hmep_bg_b[2][iplab]) << endl;

    cout << "AFT(ip="<<iplab<<"): pi0pi0jpsi = " << nevt_sim[3][iplab]/count_jpsim( hmep_bg_a[3][iplab]) << endl;
    cout << "AFT(ip="<<iplab<<"): pi0pipm = " <<  nevt_sim[0][iplab]/count_jpsim( hmep_bg_a[0][iplab]) << endl;
    cout << "AFT(ip="<<iplab<<"): pi02pipm = " << nevt_sim[1][iplab]/count_jpsim( hmep_bg_a[1][iplab]) << endl;
    cout << "AFT(ip="<<iplab<<"): pi0pipm2 = " << nevt_sim[2][iplab]/count_jpsim( hmep_bg_a[2][iplab]) << endl;

  }

  TCanvas *tc_bkf = new TCanvas("invm_bgsrc_bef_kinfit", "tc_bkf", 1400, 700);
  tc_bkf->Divide(3,1);
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_bkf->cd(1+iplab);
    tc_bkf->GetPad(1+iplab)->SetRightMargin(0.001);
    tc_bkf->GetPad(1+iplab)->SetLeftMargin(0.1);
    hmep_fg_b[iplab]->Draw();
    hmep_sg_b[iplab]->Draw("same");
    for (int ibg=0; ibg < nbg; ++ibg) hmep_bg_b[ibg][iplab]->Draw("same");
    //gPad->SetLogy();
    if (iplab==2) tl4->Draw();
  }
  tc_bkf->Print("invm_bgsrc_bef_kinfit.pdf");

  TCanvas *tc_akf = new TCanvas("invm_bgsrc_aft_kinfit", "tc_akf", 1400, 700);
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
  tc_akf->Print(iaft==6?"invm_bgsrc_aft_chi2sig_cut.pdf":"invm_bgsrc_aft_kincuts.pdf");

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
    //gPad->SetLogy();
    //tc[iplab]->Print(Form("%s/marc_meeting/minv_sig_and_bghyp_pessimistic_norm_p%d.eps",bdir,iplab));
  }

  TCanvas *tc_chi2[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_chi2[iplab] = new TCanvas(Form("tc_chi2_p%d",iplab),Form("tc_chi2_p%d",iplab));

    tc_chi2[iplab]->cd();
    for (int ii=0; ii < nbg+1; ++ii) {
      if (iplab==0)
	hpi0jpsi_chi24c[ii][iplab]->SetTitle(Form("#chi^{2} of SIG hyp. 4C kinematic fit (p_{#bar{p}} = %5.4g GeV/c); #chi^{2}; arb. unit",plab[iplab]));
      else
	hpi0jpsi_chi24c[ii][iplab]->SetTitle(Form("#chi^{2} of SIG hyp. 4C kinematic fit (p_{#bar{p}} = %4.1f GeV/c); #chi^{2}; arb. unit",plab[iplab]));
      hpi0jpsi_chi24c[ii][iplab]->GetXaxis()->SetNdivisions(505);
      hpi0jpsi_chi24c[ii][iplab]->DrawNormalized(ii==0?"":"same");
    }
    gPad->SetLogy();
    tl->Draw();

    tc_chi2[iplab]->cd();
    TPad *pad = new TPad("inset","inset",0.4,0.4,0.9,0.9);
    pad->SetRightMargin(0.02);
    pad->SetLeftMargin(0.18);
    pad->SetFillStyle(4000);
    pad->Draw();
    pad->cd();
    for (int ii=0; ii < nbg+1; ++ii) {
      hpi0jpsi_chi24cc[ii][iplab]->SetTitle("");
      hpi0jpsi_chi24cc[ii][iplab]->GetXaxis()->SetRangeUser(0,60);
      hpi0jpsi_chi24cc[ii][iplab]->DrawNormalized(ii==0?"":"same");
    }
    //tc_chi2[iplab]->Print(Form("%s/marc_meeting/chi2_sighyp_p%d.eps",bdir,iplab));
    tc_chi2[iplab]->Print(Form("kinfit_4c_chi2_dists_p%d.pdf",iplab));
  }

  return;


  TCanvas *tc_chi2bg[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_chi2bg[iplab] = new TCanvas(Form("tc_chi2bg_p%d",iplab),Form("tc_chi2bg_p%d",iplab));
    tc_chi2bg[iplab]->cd();
    if (iplab==0)
      hpi0pi0jpsi_chi24c[4][iplab]->SetTitle(Form("#chi^{2} of BG hyp. 4C kinematic fit (p_{#bar{p}} = %5.4g GeV/c); #chi^{2}; arb. unit",plab[iplab]));
    else
      hpi0pi0jpsi_chi24c[4][iplab]->SetTitle(Form("#chi^{2} of BG hyp. 4C kinematic fit (p_{#bar{p}} = %4.1f GeV/c); #chi^{2}; arb. unit",plab[iplab]));
    hpi0pi0jpsi_chi24c[4][iplab]->GetXaxis()->SetNdivisions(505);
    hpi0pi0jpsi_chi24c[4][iplab]->DrawNormalized();
    hpi0pi0jpsi_chi24c[0][iplab]->DrawNormalized("same");
    tl3->Draw();
  }


  TCanvas *tc_chi2bgvssig[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    tc_chi2bgvssig[iplab] = new TCanvas(Form("tc_chi2bgvssig_p%d",iplab),Form("tc_chi2bgvssig_p%d",iplab),1300,700);
    tc_chi2bgvssig[iplab]->Divide(2,1);
    if (iplab==0) {
      hpi0vs2pi0_chi24c[4][iplab]->SetTitle(Form("#chi^{2} Bg. vs Sig. hyp. (2(#pi^{0})J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %5.4g GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
      hpi0vs2pi0_chi24c[0][iplab]->SetTitle(Form("#chi^{2} Bg. vs Sig. hyp. (#pi^{0}J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %5.4g GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
    } else {
      hpi0vs2pi0_chi24c[4][iplab]->SetTitle(Form("#chi^{2} Bg. vs Sig. hyp. (2(#pi^{0})J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %4.1f GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
      hpi0vs2pi0_chi24c[0][iplab]->SetTitle(Form("#chi^{2} Bg. vs Sig. hyp. (#pi^{0}J/#psi(#rightarrow e^{+}e^{-}), p_{#bar{p}} = %4.1f GeV/c); #chi^{2} (Sig. Hyp.); #chi^{2} (Bg. Hyp.)    ",plab[iplab]));
    }
    TGaxis::SetMaxDigits(3);
    hpi0vs2pi0_chi24c[4][iplab]->GetXaxis()->SetNdivisions(505);
    hpi0vs2pi0_chi24c[4][iplab]->GetYaxis()->SetNdivisions(504);
    hpi0vs2pi0_chi24c[0][iplab]->GetXaxis()->SetNdivisions(505);
    hpi0vs2pi0_chi24c[0][iplab]->GetYaxis()->SetNdivisions(504);

    tc_chi2bgvssig[iplab]->cd(1);
    hpi0vs2pi0_chi24c[0][iplab]->Draw("box");
    tc_chi2bgvssig[iplab]->cd(2);
    hpi0vs2pi0_chi24c[4][iplab]->Draw("box");
    tc_chi2bgvssig[iplab]->Print(Form("%s/marc_meeting/chi2_sig_vs_bg_hyp_p%d.eps",bdir,iplab));
  }

}

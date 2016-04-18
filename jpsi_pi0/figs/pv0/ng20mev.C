#include "ananote.C"

TLegend *get_legend(TH1F*hsg, TH1F*hbg, TString hdr) {
  TLegend *tl = new TLegend(0.5, 0.55, 0.9, 0.85 );
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);
  tl->SetTextSize(0.05);
  tl->SetHeader(hdr);
  tl->AddEntry(hsg,"#pi^{0}J/#psi (#rightarrow e^{+}e^{-})","l");
  tl->AddEntry(hbg,"#pi^{0}#pi^{0}J/#psi (#rightarrow e^{+}e^{-})","l");
  return tl;
}

void ng20mev() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  bool verbose = false;
  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  int iplab = 0;
  int pass = 24;
  TFile *fsg = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_%s4ng20mev_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
  TFile *fbg = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pi0jpsi_%s4ng20mev_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));

  TH1F *hng20_nocut_sg = (TH1F*) fsg->Get("hng20mev_nocut");
  TH1F *hng20_nocut_bg = (TH1F*) fbg->Get("hng20mev_nocut");

  TH1F *hng20_excl_sg = (TH1F*) fsg->Get("hng20mev_excl");
  TH1F *hng20_excl_bg = (TH1F*) fbg->Get("hng20mev_excl");

  TH1F *hng20_chi2sg_sg = (TH1F*) fsg->Get("hng20mev_chi2sig");
  TH1F *hng20_chi2sg_bg = (TH1F*) fbg->Get("hng20mev_chi2sig");

  TH1F *hng20_chi2bg_sg = (TH1F*) fsg->Get("hng20mev_chi2bg");
  TH1F *hng20_chi2bg_bg = (TH1F*) fbg->Get("hng20mev_chi2bg");

  TH1F *hng20_kinall_sg = (TH1F*) fsg->Get("hng20mev_kinall");
  TH1F *hng20_kinall_bg = (TH1F*) fbg->Get("hng20mev_kinall");

  TCanvas *tc_nocut = new TCanvas("tc_nocut", "tc_nocut");
  tc_nocut->cd();
  set_style(hng20_nocut_sg,1,0,false);
  //hng20_nocut_sg->SetTitle(Form("N_{#gamma>20MeV}/evt (p_{#bar{p}}=%3.1f GeV/c);N_{#gamma>20MeV}/evt;arb. unit",plab[iplab]));
  hng20_nocut_sg->SetTitle(";N_{#gamma>20MeV};arb. unit");
  hng20_nocut_sg->DrawNormalized();
  set_style(hng20_nocut_bg,2,0,false);
  hng20_nocut_bg->DrawNormalized("same");
  get_legend(hng20_nocut_sg,hng20_nocut_bg, Form(" All events, p_{#bar{p}} = %3.1f GeV/c", plab[iplab]))->Draw();
  tc_nocut->Print("ng20mev_allevents_p0.pdf");
  return;

  TCanvas *tc_excl = new TCanvas("tc_excl", "tc_excl");
  tc_excl->cd();
  set_style(hng20_excl_sg,1,0,false);
  hng20_excl_sg->SetTitle(Form("N_{#gamma>20MeV}/evt (p_{#bar{p}}=%3.1f GeV/c);N_{#gamma>20MeV}/evt;arb. unit",plab[iplab]));
  hng20_excl_sg->DrawNormalized();
  set_style(hng20_excl_bg,2,0,false);
  hng20_excl_bg->DrawNormalized("same");
  get_legend(hng20_excl_sg,hng20_excl_bg, "Most back-to-back")->Draw();

  TCanvas *tc_chi2sg = new TCanvas("tc_chi2sg", "tc_chi2sg");
  tc_chi2sg->cd();
  set_style(hng20_chi2sg_sg,1,0,false);
  hng20_chi2sg_sg->SetTitle(Form("N_{#gamma>20MeV}/evt (p_{#bar{p}}=%3.1f GeV/c);N_{#gamma>20MeV}/evt;arb. unit",plab[iplab]));
  hng20_chi2sg_sg->DrawNormalized();
  set_style(hng20_chi2sg_bg,2,0,false);
  hng20_chi2sg_bg->DrawNormalized("same");
  get_legend(hng20_chi2sg_sg,hng20_chi2sg_bg, "#chi^{2}_{2#gamma e^{+}e^{-}} cut")->Draw();

  TCanvas *tc_chi2bg = new TCanvas("tc_chi2bg", "tc_chi2bg");
  tc_chi2bg->cd();
  set_style(hng20_chi2bg_sg,1,0,false);
  hng20_chi2bg_sg->SetTitle(Form("N_{#gamma>20MeV}/evt (p_{#bar{p}}=%3.1f GeV/c);N_{#gamma>20MeV}/evt;arb. unit",plab[iplab]));
  hng20_chi2bg_sg->DrawNormalized();
  set_style(hng20_chi2bg_bg,2,0,false);
  hng20_chi2bg_bg->DrawNormalized("same");
  get_legend(hng20_chi2bg_sg,hng20_chi2bg_bg, "After #chi^{2}_{2#gamma e^{+}e^{-}} and #chi^{2}_{4#gamma e^{+}e^{-}} cuts")->Draw();
  //tc_chi2bg->Print("ng20mev_aft_chi2bg.pdf");

  TCanvas *tc_kinall = new TCanvas("tc_kinall", "tc_kinall");
  tc_kinall->cd();
  set_style(hng20_kinall_bg,2,0,false);
  hng20_kinall_bg->SetTitle(Form("N_{#gamma>20MeV}/evt (p_{#bar{p}}=%3.1f GeV/c);N_{#gamma>20MeV}/evt;arb. unit",plab[iplab]));
  hng20_kinall_bg->DrawNormalized();
  set_style(hng20_kinall_sg,1,0,false);
  hng20_kinall_sg->DrawNormalized("same");
  get_legend(hng20_chi2bg_sg,hng20_chi2bg_bg, "#chi^{2} and N_{#gamma>20MeV} cut")->Draw();

}

void set_style(TAxis *axis) {
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.05);
}

void mom_vs_the_epi() {
  gStyle->SetOptStat(0);

  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  bool eff = false;
  TFile *f_bg;
  if (eff)
    f_bg = TFile::Open("hists/bg_eff.root");
  else
    f_bg = TFile::Open("hists/bg_noeff.root");
  bg(f_bg,eff);

  TFile *f_sig;
  if (eff)
    f_sig = TFile::Open("hists/sig_eff.root");
  else
    f_sig = TFile::Open("hists/sig_noeff.root");
  sig(f_sig,eff);

}

void bg(TFile *f, bool eff) {

  TH1F *pipm = (TH1F*) f->Get("lab_mom_pip_the_pip")->Clone("lab_mom_pipm_the_pipm");
  pipm->Add((TH1F*) f->Get("lab_mom_pim_the_pim"));
  set_style(pipm->GetXaxis());
  set_style(pipm->GetYaxis());
  pipm->SetTitle("#pi^{#pm} #theta vs. momentum;p_{#pi^{#pm}} [GeV/c]; #theta_{#pi^{#pm}}[#circ]");
  TCanvas *tc_bg = new TCanvas("tc_bg","tc_bg",700,700);

  pipm->Draw("colz");

  TText *tt = new TText();
  tt->SetTextSize(0.08);
  //tt->SetTextColor(2);
  tt->DrawTextNDC(0.3,0.6,"Background");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.38,0.53,"(#bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-})");

  if (eff)
    tc_bg->Print("figs/bg_mom_vs_the_dists_pipm.png");
  else
    tc_bg->Print("figs/bg_mom_vs_the_dists_pipm_noeff.png");

}


void sig(TFile *f, bool eff) {

  TH1F *epm = (TH1F*) f->Get("lab_mom_ep_the_ep")->Clone("lab_mom_epm_the_epm");
  epm->Add((TH1F*) f->Get("lab_mom_em_the_em"));
  set_style(epm->GetXaxis());
  set_style(epm->GetYaxis());
  epm->SetTitle("e^{#pm} #theta vs. momentum;p_{e^{#pm}} [GeV/c]; #theta_{e^{#pm}}[#circ]");

  TCanvas *tc_sig = new TCanvas("tc_sig","tc_sig",700,700);

  epm->Draw("colz");

  TText *tt = new TText();
  tt->SetTextSize(0.08);
  //tt->SetTextColor(2);
  tt->DrawTextNDC(0.3,0.6,"Signal");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.38,0.53,"(#bar{p}p#rightarrow#pi^{0}J/#psi)");

  if (eff)
    tc_sig->Print("figs/sig_mom_vs_the_dists_epm.png");
  else
    tc_sig->Print("figs/sig_mom_vs_the_dists_epm_noeff.png");

}

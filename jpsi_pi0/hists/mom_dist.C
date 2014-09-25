void set_style(TAxis *axis) {
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.05);
}

void mom_dist() {

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

  TH1F *pi0 = (TH1F*) f->Get("lab_mom_pi0");
  TH1F *pipm_fwd = (TH1F*) f->Get("lab_p_mom_fwd_pip_pim");
  TH1F *pipm_bwd = (TH1F*) f->Get("lab_p_mom_bwd_pip_pim");
  TH1F *pipm_p = (TH1F* ) f->Get("lab_p_mom_pip_pim");
  pipm_bwd->SetTitle("Momentum distributions;p[GeV/c]");
  set_style(pipm_bwd->GetXaxis());
  set_style(pipm_bwd->GetYaxis());

  pi0->SetLineColor(1);
  pipm_p->SetLineColor(9);
  pipm_fwd->SetLineColor(3);
  pipm_bwd->SetLineColor(7);

  pi0->SetLineWidth(2);
  pipm_p->SetLineWidth(2);
  pipm_fwd->SetLineWidth(2);
  pipm_bwd->SetLineWidth(2);

  TLegend *tl_bg = new TLegend(0.5,0.72,0.9,0.9);
  tl_bg->AddEntry(pi0,"#pi^{0}","pl");
  tl_bg->AddEntry(pipm_p,"#pi^{+}-#pi^{+} pair","pl");
  tl_bg->AddEntry(pipm_fwd,"Fwd going #pi^{#pm}","pl");
  tl_bg->AddEntry(pipm_bwd,"Bwd going #pi^{#pm}","pl");
  tl_bg->SetFillStyle(0);
  tl_bg->SetBorderSize(0);

  TCanvas *tc_bg = new TCanvas("tc_bg","tc_bg",700,700);

  tc_bg->SetLogy();
  if (eff)
    pipm_bwd->SetMinimum(20);
  else
    pipm_bwd->SetMinimum(2e2);
  pipm_bwd->Draw();
  pi0->Draw("same");
  pipm_p->Draw("same");
  pipm_fwd->Draw("same");

  tc_bg->Draw();
  tl_bg->Draw();

  TText *tt = new TText();
  tt->SetTextSize(0.08);
  //tt->SetTextColor(2);
  tt->DrawTextNDC(0.3,0.67,"Background");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.38,0.6,"(#bar{p}p#rightarrow#pi^{0}#pi^{#pm})");

  if (eff)
    tc_bg->Print("figs/bg_mom_dists_eff.pdf");
  else
    tc_bg->Print("figs/bg_mom_dists_noeff.pdf");

}

void sig(TFile *f, bool eff) {

  TH1F *pi0 = (TH1F*) f->Get("lab_mom_pi0");
  TH1F *epm_fwd = (TH1F*) f->Get("lab_p_mom_fwd_ep_em");
  TH1F *epm_bwd = (TH1F*) f->Get("lab_p_mom_bwd_ep_em");
  TH1F *jpsi = (TH1F* ) f->Get("lab_mom_jpsi");
  jpsi->SetTitle("Momentum distributions;p[GeV/c]");
  set_style(jpsi->GetXaxis());
  set_style(jpsi->GetYaxis());

  pi0->SetLineColor(1);
  jpsi->SetLineColor(9);
  epm_fwd->SetLineColor(3);
  epm_bwd->SetLineColor(7);

  pi0->SetLineWidth(2);
  jpsi->SetLineWidth(2);
  epm_fwd->SetLineWidth(2);
  epm_bwd->SetLineWidth(2);

  TLegend *tl_sig = new TLegend(0.3,0.7,0.7,0.9);
  tl_sig->AddEntry(pi0,"#pi^{0}","pl");
  tl_sig->AddEntry(jpsi,"J/#psi","pl");
  tl_sig->AddEntry(epm_fwd,"Fwd going e^{#pm}","pl");
  tl_sig->AddEntry(epm_bwd,"Bwd going e^{#pm}","pl");
  tl_sig->SetFillStyle(0);
  tl_sig->SetBorderSize(0);

  TCanvas *tc_sig = new TCanvas("tc_sig","tc_sig",720,10,700,700);
  tc_sig->SetLogy();
  jpsi->SetMinimum(1e2);
  jpsi->Draw();
  pi0->Draw("same");
  epm_fwd->Draw("same");
  epm_bwd->Draw("same");

  tl_sig->Draw();

  TText *tt = new TText();
  tt->SetTextSize(0.08);
  //tt->SetTextColor(4);
  tt->DrawTextNDC(0.7,0.8,"Signal");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.7,0.73,"(#bar{p}p#rightarrow#pi^{0}J/#psi)");

  if (eff)
    tc_sig->Print("figs/sig_mom_dists_eff.pdf");
  else
    tc_sig->Print("figs/sig_mom_dists_noeff.pdf");

}

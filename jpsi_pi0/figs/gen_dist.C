void set_style(TAxis *axis) {
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.05);
}

void gen_dist() {

  gStyle->SetOptStat(0);

  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  bool eff = false;
  TFile *f_bg;
  if (eff)
    f_bg = TFile::Open("hists/bg_eff.root");
  else
    f_bg = TFile::Open("hists/bg_noeff.root");
  //bg(f_bg,eff);
  bg_cm(f_bg,eff);

  TFile *f_sig;
  if (eff)
    f_sig = TFile::Open("hists/sig_eff.root");
  else
    f_sig = TFile::Open("hists/sig_noeff.root");
  //sig(f_sig,eff);
  sig_cm(f_sig,eff);

}

void bg(TFile *f, bool eff) {

  TH1F *pi0 = (TH1F*) f->Get("lab_the_pi0");
  TH1F *pipm_fwd = (TH1F*) f->Get("lab_p_the_fwd_pip_pim");
  TH1F *pipm_bwd = (TH1F*) f->Get("lab_p_the_bwd_pip_pim");
  TH1F *pipm_p = (TH1F* ) f->Get("lab_p_the_pip_pim");
  pipm_p->SetTitle("#theta distributions;#theta[#circ]");
  set_style(pipm_p->GetXaxis());
  set_style(pipm_p->GetYaxis());

  pi0->SetLineColor(1);
  pipm_p->SetLineColor(9);
  pipm_fwd->SetLineColor(3);
  pipm_bwd->SetLineColor(7);

  pi0->SetLineWidth(2);
  pipm_p->SetLineWidth(2);
  pipm_fwd->SetLineWidth(2);
  pipm_bwd->SetLineWidth(2);

  TLegend *tl_bg = new TLegend(0.5,0.68,0.9,0.9);
  tl_bg->AddEntry(pi0,"#pi^{0}","pl");
  tl_bg->AddEntry(pipm_p,"#pi^{+}-#pi^{+} pair","pl");
  tl_bg->AddEntry(pipm_fwd,"Fwd going #pi^{+}#pi^{-}","pl");
  tl_bg->AddEntry(pipm_bwd,"Bwd going #pi^{+}#pi^{-}","pl");
  tl_bg->SetFillStyle(0);
  tl_bg->SetBorderSize(0);

  TCanvas *tc_bg = new TCanvas("tc_bg","tc_bg",700,700);

  tc_bg->SetLogy();
  if (eff)
    pipm_p->SetMinimum(20);
  else
    pipm_p->SetMinimum(2e2);
  pipm_p->Draw();
  pi0->Draw("same");
  pipm_fwd->Draw("same");
  pipm_bwd->Draw("same");

  tc_bg->Draw();
  tl_bg->Draw();

  TText *tt = new TText();
  tt->SetTextSize(0.08);
  //tt->SetTextColor(2);
  tt->DrawTextNDC(0.3,0.6,"Background");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.38,0.53,"(#bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-})");

  if (eff)
    tc_bg->Print("figs/bg_ang_dists_eff.pdf");
  else
    tc_bg->Print("figs/bg_ang_dists_noeff.pdf");

}

void bg_cm(TFile *f, bool eff) {

  TH1F *pi0 = (TH1F*) f->Get("cm_the_pi0");
  TH1F *pipm_p = (TH1F* ) f->Get("cm_p_the_pip_pim");
  pi0->SetTitle("#pi^{0} #theta* distribution (backward kinematics);#pi^{0}#theta*[#circ]");
  set_style(pipm_p->GetXaxis());
  set_style(pipm_p->GetYaxis());

  pi0->GetXaxis()->SetRangeUser(50,180);
  pipm_p->GetXaxis()->SetRangeUser(50,180);

  pi0->SetLineColor(1);
  pipm_p->SetLineColor(9);

  pi0->SetLineWidth(2);
  pipm_p->SetLineWidth(2);

  TLegend *tl_bg = new TLegend(0.5,0.68,0.9,0.9);
  tl_bg->AddEntry(pi0,"#pi^{0}","pl");
  tl_bg->AddEntry(pipm_p,"#pi^{+}-#pi^{+} pair","pl");
  tl_bg->SetFillStyle(0);
  tl_bg->SetBorderSize(0);

  TCanvas *tc_bg_cm = new TCanvas("tc_bg_cm","tc_bg_cm",700,700);

  //tc_bg_cm->SetLogy();
  if (eff)
    pipm_p->SetMinimum(20);
  else
    pipm_p->SetMinimum(2e2);
  //pipm_p->Draw();
  pi0->Draw();

  tc_bg_cm->Draw();
  //tl_bg->Draw();

  TText *tt = new TText();
  tt->SetTextSize(0.08);
  //tt->SetTextColor(2);
  tt->DrawTextNDC(0.37,0.4,"Background");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.48,0.3,"(#bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-})");

  if (eff)
    tc_bg_cm->Print("figs/bg_cm_ang_dists_eff.pdf");
  else
    tc_bg_cm->Print("figs/bg_cm_ang_dists_noeff.pdf");

}


void sig(TFile *f, bool eff) {

  TH1F *pi0 = (TH1F*) f->Get("lab_the_pi0");
  TH1F *epm_fwd = (TH1F*) f->Get("lab_p_the_fwd_ep_em");
  TH1F *epm_bwd = (TH1F*) f->Get("lab_p_the_bwd_ep_em");
  TH1F *jpsi = (TH1F* ) f->Get("lab_the_jpsi");
  jpsi->SetTitle("#theta distributions;#theta[#circ]");
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

  TLegend *tl_sig = new TLegend(0.5,0.65,0.9,0.9);
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
  tt->DrawTextNDC(0.54,0.54,"Signal");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.51,0.47,"(#bar{p}p#rightarrow#pi^{0}J/#psi)");

  if (eff)
    tc_sig->Print("figs/sig_ang_dists_eff.pdf");
  else
    tc_sig->Print("figs/sig_ang_dists_noeff.pdf");

}


void sig_cm(TFile *f, bool eff) {

  TH1F *pi0 = (TH1F*) f->Get("cm_the_pi0");
  TH1F *jpsi = (TH1F* ) f->Get("cm_the_jpsi");
  pi0->SetTitle("#pi^{0} #theta* distribution (backward kinematics);#pi^{0}#theta*[#circ]");
  set_style(jpsi->GetXaxis());
  set_style(jpsi->GetYaxis());

  pi0->SetLineColor(1);
  jpsi->SetLineColor(2);

  pi0->GetXaxis()->SetRangeUser(60,180);
  jpsi->GetXaxis()->SetRangeUser(60,180);

  pi0->SetLineWidth(2);
  jpsi->SetLineWidth(2);

  TLegend *tl_sig = new TLegend(0.5,0.65,0.9,0.9);
  tl_sig->AddEntry(pi0,"#pi^{0}","pl");
  tl_sig->AddEntry(jpsi,"J/#psi","pl");
  tl_sig->SetFillStyle(0);
  tl_sig->SetBorderSize(0);

  TCanvas *tc_sig_cm = new TCanvas("tc_sig_cm","tc_sig_cm",720,10,700,700);
  //tc_sig_cm->SetLogy();
  jpsi->SetMinimum(1e2);
  //jpsi->Draw();
  pi0->Draw();

  TText *tt = new TText();
  tt->SetTextSize(0.08);
  //tt->SetTextColor(4);
  tt->DrawTextNDC(0.54,0.54,"Signal");
  TLatex *tt2 = new TLatex();
  tt2->SetNDC();
  tt2->SetTextSize(0.06);
  tt2->DrawLatex(0.51,0.47,"(#bar{p}p#rightarrow#pi^{0}J/#psi)");

  //tl_sig->Draw();

  if (eff)
    tc_sig_cm->Print("figs/sig_cm_ang_dists_eff.pdf");
  else
    tc_sig_cm->Print("figs/sig_cm_ang_dists_noeff.pdf");

}

void sig_bg_cm(TFile *f_sig, TFile *f_bg, bool eff) {


}

void mom_dist() {

  gStyle->SetOptStat(0);

  TFile *f_bg = TFile::Open("hists/bg_noeff.root");
  bg(f_bg);

  TFile *f_sig = TFile::Open("hists/sig_noeff.root");
  sig(f_sig);

}

void bg(TFile *f) {

  TH1F *pi0 = (TH1F*) f->Get("lab_mom_pi0");
  TH1F *pipm_fwd = (TH1F*) f->Get("lab_p_mom_fwd_pip_pim");
  TH1F *pipm_bwd = (TH1F*) f->Get("lab_p_mom_bwd_pip_pim");
  TH1F *pipm_p = (TH1F* ) f->Get("lab_p_mom_pip_pim");
  pipm_bwd->SetTitle("Momentum distributions, Lab frame;p[GeV/c]");

  pi0->SetLineColor(1);
  pipm_p->SetLineColor(2);
  pipm_fwd->SetLineColor(3);
  pipm_bwd->SetLineColor(7);

  pi0->SetLineWidth(2);
  pipm_p->SetLineWidth(2);
  pipm_fwd->SetLineWidth(2);
  pipm_bwd->SetLineWidth(2);

  TLegend *tl_bg = new TLegend(0.2,0.15,0.6,0.4);
  tl_bg->AddEntry(pi0,"#pi^{0}","pl");
  tl_bg->AddEntry(pipm_p,"#pi^{+}-#pi^{+} pair","pl");
  tl_bg->AddEntry(pipm_fwd,"Fwd going #pi^{#pm}","pl");
  tl_bg->AddEntry(pipm_bwd,"Bwd going #pi^{#pm}","pl");
  tl_bg->SetFillStyle(0);
  tl_bg->SetBorderSize(0);

  TCanvas *tc_bg = new TCanvas("tc_bg","tc_bg");
  //tc_bg->SetLogy();
  pipm_bwd->Draw();
  pi0->Draw("same");
  pipm_p->Draw("same");
  pipm_fwd->Draw("same");

  tc_bg->Draw();
  tl_bg->Draw();

}

void sig(TFile *f) {

  TH1F *pi0 = (TH1F*) f->Get("lab_mom_pi0");
  TH1F *epm_fwd = (TH1F*) f->Get("lab_p_mom_fwd_ep_em");
  TH1F *epm_bwd = (TH1F*) f->Get("lab_p_mom_bwd_ep_em");
  TH1F *jpsi = (TH1F* ) f->Get("lab_mom_jpsi");
  jpsi->SetTitle("Momentum distributions, Lab frame;p[GeV/c]");

  pi0->SetLineColor(1);
  jpsi->SetLineColor(2);
  epm_fwd->SetLineColor(3);
  epm_fwd->SetLineColor(7);

  pi0->SetLineWidth(2);
  jpsi->SetLineWidth(2);
  epm_fwd->SetLineWidth(2);
  epm_fwd->SetLineWidth(2);

  TLegend *tl_sig = new TLegend(0.3,0.55,0.7,0.8);
  tl_sig->AddEntry(pi0,"#pi^{0}","pl");
  tl_sig->AddEntry(jpsi,"J/#psi","pl");
  tl_sig->AddEntry(epm_fwd,"Fwd going e^{#pm}","pl");
  tl_sig->AddEntry(epm_bwd,"Bwd going e^{#pm}","pl");
  tl_sig->SetFillStyle(0);
  tl_sig->SetBorderSize(0);

  TCanvas *tc_sig = new TCanvas("tc_sig","tc_sig",720,10,700,500);
  tc_sig->SetLogy();
  jpsi->Draw();
  pi0->Draw("same");
  epm_fwd->Draw("same");
  epm_bwd->Draw("same");

  tc_sig->Draw();
  tl_sig->Draw();

}

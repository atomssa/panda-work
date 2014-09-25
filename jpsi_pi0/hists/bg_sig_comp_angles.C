void set_style(TAxis *axis) {
  axis->SetTitleSize(0.05);
  //axis->SetTitleSize(0.05,"Y");
  //axis->SetTitleSize(0.05,"Z");
  //axis->SetTitleSize(0.05,"T");
  axis->SetLabelSize(0.05);
  //axis->SetLabelSize(0.05,"Y");
  //axis->SetLabelSize(0.05,"Z");
  //axis->SetLabelSize(0.05,"T");
}

void bg_sig_comp_angles(){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  bool eff= true;

  TFile *f_bg;
  if (eff) {
    f_bg = TFile::Open("hists/bg_eff.root");
  } else {
    f_bg = TFile::Open("hists/bg_noeff.root");
  }

  TH1F *bg_goa = (TH1F*)f_bg->Get("lab_p_oa_g1_g2");
  TH1F *bg_pi0a = (TH1F*)f_bg->Get("lab_the_pi0");
  bg_goa->SetLineWidth(2);
  bg_pi0a->SetLineWidth(2);
  bg_goa->SetLineColor(2);
  bg_pi0a->SetLineColor(2);

  TFile *f_sig;
  double scale = 1;
  if (eff) {
    f_sig = TFile::Open("hists/sig_eff.root");
    //scale = 0.07;
  } else {
    f_sig = TFile::Open("hists/sig_noeff.root");
    //scale = 0.7;
  }
  TH1F *sig_goa = (TH1F*)f_sig->Get("lab_p_oa_g1_g2");
  TH1F *sig_pi0a = (TH1F*)f_sig->Get("lab_the_pi0");
  sig_goa->SetLineWidth(2);
  sig_pi0a->SetLineWidth(2);
  sig_goa->Scale(scale);
  sig_pi0a->Scale(scale);

  TLegend *tl = new TLegend(0.4,0.6,0.9,0.9);
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);
  tl->AddEntry(bg_goa, "Background","pl");
  tl->AddEntry(sig_goa, "Signal","pl");

  TCanvas *oa = new TCanvas("oa","oa",700,700);
  oa->cd();
  gPad->SetLogy();
  bg_goa->SetMinimum(10);
  if (eff) {
    set_style(sig_goa->GetXaxis());
    set_style(sig_goa->GetYaxis());
    sig_goa->Draw();
    bg_goa->Draw("same");
  } else {
    set_style(bg_goa->GetXaxis());
    set_style(bg_goa->GetYaxis());
    bg_goa->Draw();
    sig_goa->Draw("same");
  }
  tl->Draw();

  TCanvas *th = new TCanvas("th","th",720,10,700,700);
  th->cd();
  gPad->SetLogy();
  set_style(sig_pi0a->GetXaxis());
  set_style(sig_pi0a->GetYaxis());
  sig_pi0a->SetMinimum(10);
  sig_pi0a->Draw();
  bg_pi0a->Draw("same");
  tl->Draw();

  if (eff) {
    oa->Print("figs/bg_sig_comp_gg_oa_eff.pdf");
    th->Print("figs/bg_sig_comp_pi0_the_eff.pdf");
  } else {
    oa->Print("figs/bg_sig_comp_gg_oa_noeff.pdf");
    th->Print("figs/bg_sig_comp_pi0_the_noeff.pdf");
  }

}

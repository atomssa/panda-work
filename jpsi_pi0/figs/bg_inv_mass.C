void set_style(TAxis *axis) {
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.05);
}

void bg_inv_mass(){

  gStyle->SetOptStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  TFile *f_bg = TFile::Open("hists/bg_nofilt_noeff.root");

  TH1F *pip_pim = (TH1F*) f_bg->Get("inv_p_mass_pip_pim");
  pip_pim->SetLineWidth(2);
  pip_pim->SetTitle("#pi^{+}-#pi^{+} invariant mass");
  set_style(pip_pim->GetXaxis());
  set_style(pip_pim->GetYaxis());
  TH1F *pi0_pip = (TH1F*) f_bg->Get("inv_p_mass_pi0_pip");
  pi0_pip->SetTitle("#pi^{0}-#pi^{+} and #pi^{0}-#pi^{-} invariant mass");
  pi0_pip->SetLineWidth(2);
  pi0_pip->SetLineColor(2);
  set_style(pi0_pip->GetXaxis());
  set_style(pi0_pip->GetYaxis());
  TH1F *pi0_pim = (TH1F*) f_bg->Get("inv_p_mass_pi0_pim");
  //pi0_pim->SetTitle("Invariant Mass;M_{inv}[GeV/c^{2}]");
  pi0_pim->SetTitle("#pi^{0}-#pi^{+} and #pi^{0}-#pi^{-} invariant mass");
  pi0_pim->SetLineWidth(2);
  pi0_pim->SetLineColor(3);
  set_style(pi0_pim->GetXaxis());
  set_style(pi0_pim->GetYaxis());

  TLegend *tl = new TLegend(0.5,0.6,0.9,0.9);
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);
  tl->AddEntry(pip_pim, "#pi^{+}-#pi^{-}","pl");

  TLegend *tl2 = new TLegend(0.5,0.6,0.9,0.9);
  tl2->SetFillStyle(0);
  tl2->SetBorderSize(0);
  tl2->AddEntry(pi0_pim, "#pi^{0}-#pi^{-}","pl");
  tl2->AddEntry(pi0_pip, "#pi^{0}-#pi^{+}","pl");

  TCanvas *tc = new TCanvas("tc","tc");
  tc->cd();

  pip_pim->Draw();
  TLine *tll = new TLine(2.96,0,2.96,30000);
  TLine *tlu = new TLine(3.22,0,3.22,30000);

  tll->SetLineWidth(2);
  tlu->SetLineWidth(2);

  tll->Draw();
  tlu->Draw();

  tl->Draw();
  tc->Print("figs/bg_inv_mass_pip_pim.pdf");

  TCanvas *tc2 = new TCanvas("tc2","tc2");
  tc2->cd();
  pi0_pip->Draw();
  pi0_pim->Draw("same");
  tl2->Draw();
  tc2->Print("figs/bg_inv_mass_pi0_pipm.pdf");

  TCanvas *td = new TCanvas("td","td");
  TH2F *dal = (TH2F*) f_bg->Get("inv_p_mass_pi0_pip_p_mass_pip_pim");
  dal->SetTitle("#pi^{0}#pi^{+}#pi^{-} Dalitz");
  set_style(dal->GetXaxis());
  set_style(dal->GetYaxis());
  dal->Draw("colz");
  td->Print("figs/pipm_pi0_dalitz.png");

  TCanvas *td2 = new TCanvas("td2","td2");
  TH2F *dal2 = (TH2F*) f_bg->Get("inv_p_mass_sq_pi0_pip_p_mass_sq_pip_pim");
  dal2->SetTitle("#pi^{0}#pi^{+}#pi^{-} Dalitz");
  set_style(dal2->GetXaxis());
  set_style(dal2->GetYaxis());
  dal2->Draw("colz");
  td2->Print("figs/pipm_pi0_dalitz_m2.png");



}

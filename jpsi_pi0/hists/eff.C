void set_style(TAxis *axis) {
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.05);
}

void eff() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  bool sig = false;

  string rxn = sig?"#pi^{0}J/#psi":"#pi^{0}#pi^{+}#pi^{-}";

  TFile *f_bef;
  if (sig)
    f_bef = TFile::Open("hists/sig_noeff.root");
  else
    f_bef = TFile::Open("hists/bg_noeff.root");
  TH1F *pi0b = (TH1F*) f_bef->Get("lab_the_pi0");
  pi0b->SetTitle(Form("#theta_{#pi^{0}} in #bar{p}p#rightarrow%s", rxn.c_str()));
  if (!sig)
    pi0b->SetMinimum(1e-6);

  TFile* f_aft;
  if (sig)
    f_aft = TFile::Open("hists/sig_eff.root");
  else
    f_aft = TFile::Open("hists/bg_eff.root");
  TH1F *pi0a = (TH1F*) f_aft->Get("lab_the_pi0");
  pi0a->SetTitle(Form("#theta_{#pi^{0}} in #bar{p}p#rightarrow%s", rxn.c_str()));

  pi0a->SetLineWidth(2);
  pi0a->SetLineColor(2);
  set_style(pi0a->GetXaxis());
  set_style(pi0a->GetYaxis());

  pi0b->SetLineWidth(2);
  pi0b->SetLineColor(4);
  set_style(pi0b->GetXaxis());
  set_style(pi0b->GetYaxis());

  TCanvas *tc0 = new TCanvas("tc0","tc0",700,700);
  if (!sig) gPad->SetLogy();
  pi0b->Draw();
  pi0a->Draw("same");

  TLegend *tl;
  if (sig)
    tl = new TLegend(0.5,0.65,0.9,0.8);
  else
    tl = new TLegend(0.2,0.55,0.6,0.7);
  tl->AddEntry(pi0b,"Before Eff. Weight","pl");
  tl->AddEntry(pi0a,"After Eff. Weight","pl");
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);

  TLatex *tt = new TLatex();
  tt->SetNDC();
  tt->SetTextSize(0.06);
  double eff = double(pi0a->Integral())/double(pi0b->Integral());
  double err = sqrt(eff*(1-eff))/sqrt( double(pi0b->Integral()));
  cout << "eff(1-eff) = " <<  eff*(1-eff)<< endl;
  cout << "sqrt(eff(1-eff)) = " << sqrt(eff*(1-eff))<< endl;
  cout << "sqrt(N) = " << sqrt( double(pi0b->Integral())) << endl;
  cout << "err = " << err << endl;

  // if (sig)
  //   tt->DrawLatex(0.4,0.83,Form("#varepsilon^{tot}_{%s} #approx (%4.2f#pm%3.2f)%%",rxn.c_str(), 100.*eff, 100.*err));
  // else
  //   tt->DrawLatex(0.15,0.73,Form("#varepsilon^{mis-id}_{%s} #approx (%3.1f#pm%3.1f)#times 10^{-8}",rxn.c_str(), 1e8*eff, 1e8*err));
  // //tt->DrawLatex(0.4,0.83,Form("#varepsilon_{%s} #approx (%4.2f#pm%3.2f)%%",rxn.c_str(),100.*eff,100.*err));

  if (sig)
    tt->DrawLatex(0.4,0.83,Form("#varepsilon^{PID}_{%s} #approx %4.2f%%",rxn.c_str(), 100.*eff));
  else
    tt->DrawLatex(0.15,0.73,Form("#varepsilon^{MIS-ID}_{%s} #approx %3.1f#times 10^{-8}",rxn.c_str(), 1e8*eff));


  TText *tt2 = new TText();
  tt2->SetTextSize(0.08);
  //tt->SetTextColor(4);
  if (sig)
    tt2->DrawTextNDC(0.7,0.57,"Signal");
  else
    tt2->DrawTextNDC(0.2,0.47,"Background");

  tl->Draw();

  if (sig)
    tc0->Print("figs/eff_sig_pi0_the.pdf");
  else
    tc0->Print("figs/eff_bg_pi0_the.pdf");

}

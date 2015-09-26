void set_style(TH1* h, int col) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);

  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);

  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void draw_eloss() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.08);
  tt->SetNDC(kTRUE);

  TFile *fpandarootg4 = TFile::Open("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.8_hists.root");
  TH1F* eloss_pandarootg4 = (TH1F*) fpandarootg4->Get("h_eloss_all");
  set_style(eloss_pandarootg4, 1);
  eloss_pandarootg4->SetTitle(";#Delta E (GeV);counts");

  TFile *fg3 = TFile::Open("../grid.out/esim_oct14_constp_g3/bremcorr.all.ibs.cfg.8_hists.root");
  TH1F* eloss_g3 = (TH1F*) fg3->Get("h_eloss_all");
  set_style(eloss_g3, 3);
  eloss_g3->SetTitle(";#Delta E (GeV);counts");

  TFile *fdefaultg4 = TFile::Open("../grid.out/esim_oct14_constp_defaultg4/bremcorr.all.ibs.cfg.8_hists.root");
  TH1F* eloss_defaultg4 = (TH1F*) fdefaultg4->Get("h_eloss_all");
  eloss_defaultg4->Scale(1.8/1.5); // some missing file
  set_style(eloss_defaultg4, 4);
  eloss_defaultg4->SetTitle(";#Delta E (GeV);counts");

  TFile *fsteplimg4 = TFile::Open("../grid.out/esim_oct14_constp_steplimonlyg4/bremcorr.all.ibs.cfg.8_hists.root");
  TH1F* eloss_steplimg4 = (TH1F*) fsteplimg4->Get("h_eloss_all");
  set_style(eloss_steplimg4, 2);
  eloss_steplimg4->SetTitle(";#Delta E (GeV);counts");

  TCanvas *tc = new TCanvas("tc","tc");
  eloss_pandarootg4->SetMinimum(10);
  eloss_pandarootg4->Draw();
  eloss_steplimg4->Draw("same");
  eloss_defaultg4->Draw("same");
  eloss_g3->Draw("same");
  tt->DrawLatex(0.2,0.8,"30#circ < #theta < 45#circ");
  tt->DrawLatex(0.2,0.7,"p = 1 GeV/c");
  gPad->SetLogy();

  TLegend *tl = new TLegend(0.5,0.55,0.9,0.85);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  tl->AddEntry(eloss_g3,"Geant3");
  tl->AddEntry(eloss_pandarootg4,"Geant4, PandaROOT");
  tl->AddEntry(eloss_steplimg4,"Geant4, stepLimiter only");
  tl->AddEntry(eloss_defaultg4,"Geant4, No Opt.");
  tl->Draw();

  tc->Print("eloss.pdf");

  double tot = 0.0;
  for (int i=1; i<eloss_pandarootg4->GetXaxis()->GetNbins(); ++i) {
    double e = eloss_pandarootg4->GetXaxis()->GetBinCenter(i);
    tot += eloss_pandarootg4->GetBinContent(i);
    double frac = tot/double(eloss_pandarootg4->Integral());
    //cout << "Frac = " << frac << " @ Eloss = " << e << endl;
    if (frac>0.5) {
      cout << "50%E = " << e << endl;
      break;
    }
  }

  tot = 0.0;
  for (int i=1; i<eloss_pandarootg4->GetXaxis()->GetNbins(); ++i) {
    double e = eloss_pandarootg4->GetXaxis()->GetBinCenter(i);
    tot += eloss_pandarootg4->GetBinContent(i);
    double frac = tot/double(eloss_pandarootg4->Integral());
    if (e >= 0.1) {
      cout << "Eloss= 100MeV -> remaining frac= " << 1.0-frac << endl;
      break;
    }
  }


}

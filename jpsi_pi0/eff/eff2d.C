void eff2d() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetTitleOffset(1.5,"Y");
  gStyle->SetTitleFontSize(0.08);

  //gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  TFile *fe = TFile::Open("hadd_out/hadd.e.root");
  TH2F *histe = (TH2F*) fe->Get("prob_cut_5/eff2d_e_id");
  histe->SetTitle("#varepsilon(e^{+},e^{-}) for combined eID prob > 50%;p_{MC}[GeV/c];#theta_{MC}[deg];");
  //histe->GetXaxis()->SetTitleSize(0.06);
  //histe->GetXaxis()->SetLabelSize(0.06);
  //histe->GetYaxis()->SetTitleSize(0.06);
  //histe->GetYaxis()->SetLabelSize(0.06);
  TCanvas *tce = new TCanvas("tce","tce",1000,1000);
  tce->cd();
  histe->Draw("colz");
  tce->Print("eff2de.png");


  TFile *fpi = TFile::Open("hadd_out/hadd.pi.root");
  TH2F *histpi = (TH2F*) fpi->Get("prob_cut_5/eff2d_e_id");
  histpi->SetTitle("#varepsilon(#pi^{+},#pi^{-}) for combined eID prob > 50%;p_{MC}[GeV/c];#theta_{MC}[deg];");
  //histpi->GetXaxis()->SetTitleSize(0.06);
  //histpi->GetXaxis()->SetLabelSize(0.06);
  //histpi->GetYaxis()->SetTitleSize(0.06);
  //histpi->GetYaxis()->SetLabelSize(0.06);
  TCanvas *tcpi = new TCanvas("tcpi","tcpi",1000,1000);
  tcpi->cd();
  histpi->Draw("colz");
  tcpi->Print("eff2dpi.png");


}

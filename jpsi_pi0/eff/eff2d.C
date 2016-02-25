void eff2d() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetTitleOffset(1.5,"Y");
  gStyle->SetTitleFontSize(0.08);

  //gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  TFile *f = TFile::Open("hadd_out/hadd.e.root");
  TH2F *hist = (TH2F*) f->Get("prob_cut_9/eff2d_e_id");

  hist->SetTitle("#varepsilon(e^{#pm}) for combined eID prob > 90%;p_{MC}[GeV/c];#theta_{MC}[deg];");

  //hist->GetXaxis()->SetTitleSize(0.06);
  //hist->GetXaxis()->SetLabelSize(0.06);
  //hist->GetYaxis()->SetTitleSize(0.06);
  //hist->GetYaxis()->SetLabelSize(0.06);

  TCanvas *tc = new TCanvas("tc","tc",1000,1000);
  tc->cd();
  hist->Draw("colz");

  tc->Print("eff2d.pdf");

}

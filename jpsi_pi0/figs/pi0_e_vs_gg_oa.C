void pi0_e_vs_gg_oa() {

  gStyle->SetOptStat(0);

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);

  //gStyle->SetPadTopMargin(0.12);
  //gStyle->SetTitleSize(0.08,"A");

  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Z");

  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelSize(0.05,"Z");


  TFile*f = TFile::Open("hists/bg_noeff.root");
  TH2F* h = (TH2F*) f->Get("lab_e_pi0_p_oa_g1_g2");
  TCanvas *tc = new TCanvas("tc","tc",700,700);
  tc->cd();
  h->Draw("col");

  //tc->Update();
  tc->UseCurrentStyle();
  tc->Print("figs/pi0_e_vs_gg_oa.png");
  gSystem->Exec("convert figs/pi0_e_vs_gg_oa.png figs/pi0_e_vs_gg_oa.pdf");

}

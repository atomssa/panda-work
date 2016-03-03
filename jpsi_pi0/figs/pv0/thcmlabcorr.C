#include "ananote.C"

void thcmlabcorr() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetLabelOffset(-2.5,"y");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  bool verbose = false;
  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  int ibrem = 1;
  int iplab = 2;
  int pass = 18;
  TFile *fsg = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_%s4eff_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));

  const char* kin = "itu3";

  TH2F* epcthcm_epthlab = (TH2F*) fsg->Get(Form("epeff/hepcosth_jpsi_vs_epthlab_rec_%s", kin))->Clone("epcthcm_epthlab");
  TH2F* epcthcm_emthlab = (TH2F*) fsg->Get(Form("epeff/hepcosth_jpsi_vs_emthlab_rec_%s", kin))->Clone("epcthcm_emthlab");
  TH2F* epcthcm_epmthlab = (TH2F*) fsg->Get(Form("epeff/hepcosth_jpsi_vs_epthlab_rec_%s", kin))->Clone("epcthcm_epmthlab");
  epcthcm_epmthlab->Add((TH2F*) fsg->Get(Form("epeff/hepcosth_jpsi_vs_emthlab_rec_%s", kin)));

  set_style(epcthcm_epthlab);
  set_style(epcthcm_emthlab);
  epcthcm_epthlab->GetYaxis()->SetTitleOffset(1.2);
  epcthcm_emthlab->GetYaxis()->SetTitleOffset(1.2);

  //TCanvas *tc = new TCanvas("epcthcm_epthlab","epcthcm_epthlab",800,400);
  TCanvas *tc = new TCanvas("epcthcm_epthlab","epcthcm_epthlab",1600,800);
  tc->Divide(2,1);

  tc->cd(1);
  epcthcm_epthlab->SetTitle(";cos(#theta^{e^{+}}_{J#psi});#theta^{e^{+}}_{lab}");
  epcthcm_epthlab->Draw("colz");
  tc->cd(2);
  epcthcm_emthlab->SetTitle(";cos(#theta^{e^{+}}_{J#psi});#theta^{e^{-}}_{lab}");
  epcthcm_emthlab->Draw("colz");

  tc->Print(Form("epcthcm_epthlab_%s_p%d.png",kin, iplab));

}

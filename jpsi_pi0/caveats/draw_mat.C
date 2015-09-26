void draw_mat() {

  gStyle->SetOptStat(0);

  TText *tt = new TText();
  tt->SetTextSize(0.08);

  // p=1GeV/c
  TFile *fbar = TFile::Open("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.8_hists.root");
  TH1F* h_rad = (TH1F*) fbar->Get("h_rad_true");
  h_rad->SetLineWidth(2);
  h_rad->SetTitle("Brem #gamma radial emission point (Barrel); R [cm]");
  TCanvas *tc_rad = new TCanvas("tc_rad","tc_rad", 1000, 1000);
  tc_rad->cd();
  h_rad->Draw();
  tt->DrawTextNDC(0.3,0.8,"Barrel");
  tt->DrawTextNDC(0.3,0.7,"p = 1 GeV/c");
  //tc_rad->Print("mat_budget_barrel_rad_dep.pdf");
  tc_rad->Print("rgamma.pdf");

  int ent1 = h_rad->GetXaxis()->FindBin(15.9);
  int ent2 = h_rad->GetXaxis()->FindBin(16);
  int exit1 = h_rad->GetXaxis()->FindBin(41);
  int exit2 = h_rad->GetXaxis()->FindBin(41.1);
  double entrance = h_rad->GetBinContent(ent1) + h_rad->GetBinContent(ent2);
  double exit = h_rad->GetBinContent(exit1) + h_rad->GetBinContent(exit2);
  double tot = h_rad->Integral();

  cout << "#Entrance= " << entrance << " = " << entrance*100/tot << "%" << endl;
  cout << "#Exit = " << exit << " = " << exit*100/tot << "%" << endl;

  TFile *ffwd = TFile::Open("../grid.out/esim_oct14_constp/bremcorr.all.ibs.constp.cfg.18_hists.root");
  TH1F* h_zed = (TH1F*) ffwd->Get("h_zed_true");
  h_zed->SetLineWidth(2);
  h_zed->SetTitle("Brem #gamma longitudial emission point (Fwd); Z [cm]");
  TCanvas *tc_zed = new TCanvas("tc_zed","tc_zed", 1000, 1000);
  tc_zed->cd();
  h_zed->Draw();
  tt->DrawTextNDC(0.3,0.8,"Fwd endcap");
  tt->DrawTextNDC(0.3,0.7,"p = 1 GeV/c");
  //tc_zed->Print("mat_budget_fwdendcap_zed_dep.pdf");
  tc_zed->Print("zgamma.pdf");

}

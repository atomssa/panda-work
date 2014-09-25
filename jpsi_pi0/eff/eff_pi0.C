void eff_pi0() {

  TF1* eff_E = new TF1("eff_E","-0.0553 * x + 0.918", 0., 10.);
  eff_E->SetLineWidth(2);
  TF1* eff_T = new TF1("eff_T","0.812 * (1-exp (-0.398 * x) )", 0., 180.);
  eff_T->SetLineWidth(2);

  TCanvas* tcE = new TCanvas("tcE","tcE");
  tcE->cd();
  eff_E->Draw();

  TCanvas* tcT = new TCanvas("tcT","tcT");
  tcT->cd();
  eff_T->Draw();

}

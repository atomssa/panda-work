{

  gStyle->SetOptStat(0);

  TFile *_file0 = TFile::Open("anav2_pi0jpsi_brem4eff_p0_pass17.root");
  TFile *_file1 = TFile::Open("anav2_pi0pi0jpsi_brem_p0_pass17.root");

  TLegend *tl0 = new TLegend(0.5,0.5,0.9,0.9);
  tl0->SetFillStyle(0);
  tl0->SetBorderSize(0);

  TCanvas *tc1 = new TCanvas("tc1","tc1");
  tc1->cd();

  _file1->cd();
  hpi0jpsi_prob4c->SetTitle("probability of #chi^{2};prob");
  hpi0jpsi_prob4c->DrawNormalized();
  hpi0jpsi_prob4c->SetLineWidth(2);
  tl0->AddEntry(hpi0jpsi_prob4c,"prob of #chi^{2} (#pi^{0}J/#psi)","l");

  _file0->cd();
  hpi0jpsi_prob4c->SetLineColor(2);
  hpi0jpsi_prob4c->SetLineWidth(2);
  hpi0jpsi_prob4c->DrawNormalized("same");
  tl0->AddEntry(hpi0jpsi_prob4c,"prob of #chi^{2} (#pi^{0}#pi^{0}J/#psi)","l");

  tl0->Draw();

  gPad->SetLogy();

  //tc1->Print("prob.png");

  TCanvas *tc2 = new TCanvas("tc2","tc2");
  tc2->cd();

  TLegend *tl = new TLegend(0.4,0.5,0.9,0.9);
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);

  _file0->cd();
  hpi0jpsi_chi24c_c->GetXaxis()->SetRangeUser(0,50);
  hpi0jpsi_chi24c_c->SetTitle("#chi^{2} dist (#pi^{0}J/#psi events); #chi^{2}");
  hpi0jpsi_chi24c_c->Draw();
  tl->AddEntry(hpi0jpsi_chi24c_c,"#chi^{2} from kin. fit","l");

  TF1 *f1 = new TF1("chi2","[1]*TMath::Power(x,[0]/2-1)*TMath::Power(TMath::E(),-x/2)/TMath::Gamma([0]/2)/TMath::Power(2,[0]/2)",0,50);
  f1->FixParameter(0,4);
  hpi0jpsi_chi24c_c->Fit("chi2");
  tl->AddEntry(f1,"ideal #chi^{2} function (ndf=4)","l");

  tl->Draw();
  //tc2->Print("chi2.png");

  TCanvas *tc3 = new TCanvas("tc3","tc3");
  tc3->cd();

  _file1->cd();
  hpi0jpsi_chi24c_c->GetXaxis()->SetRangeUser(0,500);
  hpi0jpsi_chi24c_c->SetTitle("#chi^{2} dist (#pi^{0}#pi^{0}J/#psi events); #chi^{2}");
  hpi0jpsi_chi24c_c->Rebin(5);
  hpi0jpsi_chi24c_c->Draw();

  tc3->Print("chi2_bg.png");


}

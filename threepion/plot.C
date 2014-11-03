void set_style(TH1* h, int col=0) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.2);

  if (col>0) {
    h->SetLineWidth(3);
    h->SetLineColor(col);
  }
}

plot() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);


  TFile *dpm = TFile::Open("output_DPM_PipPimPi0.dat.root");

  TFile *ftf = TFile::Open("output_FTF_PipPimPi0.dat.root");

  TH1F* the_pi0_cm_dpm = (TH1F*) dpm->Get("the_pi0_cm_DPM");
  TH1F* the_pi0_cm_ftf = (TH1F*) ftf->Get("the_pi0_cm_FTF");
  TH1F* minv_pippim_dpm = (TH1F*) dpm->Get("minv_pippim_lab_DPM");
  TH1F* minv_pippim_ftf = (TH1F*) ftf->Get("minv_pippim_lab_FTF");

  int iRebin=2;
  the_pi0_cm_dpm->Rebin(iRebin);
  the_pi0_cm_ftf->Rebin(iRebin);
  minv_pippim_dpm->Rebin(iRebin);
  minv_pippim_ftf->Rebin(iRebin);

  set_style(minv_pippim_dpm);
  set_style(minv_pippim_ftf);

  the_pi0_cm_dpm->SetLineColor(4);
  minv_pippim_dpm->SetLineColor(4);

  TLegend *tl_the = new TLegend(0.35,0.6,0.7,0.9);
  tl_the->AddEntry(the_pi0_cm_dpm,"#theta^{*}_{#pi^{0}} in DPM","l");
  tl_the->AddEntry(the_pi0_cm_ftf,"#theta^{*}_{#pi^{0}} in FTF","l");
  tl_the->SetBorderSize(0);
  tl_the->SetFillStyle(0);


  TCanvas *minv = new TCanvas("minv","minv",700,700);
  minv->cd();
  minv_pippim_ftf->SetTitle(";M_{inv}^{#pi^{+}#pi^{-}} [GeV/c^{2}]; dN/dM_{inv}");
  minv_pippim_ftf->Draw();
  minv_pippim_dpm->Draw("same");
  minv_pippim_dpm->SetTitle("");

  TLatex *tt = new TLatex();
  tt->DrawLatexNDC(0.25,0.3,"#bar{p}p#rightarrow#pi^{+}#pi^{-}#pi^{0}");
  tt->SetTextSize(0.035);
  tt->DrawLatexNDC(0.25,0.25,"#sqrt{s} = 3.5 GeV");
  //2.96 < mom_pich.M() && mom_pich.M() < 3.22
  TLine *tlb = new TLine(2.96,0,2.96,0.65*minv_pippim_dpm->GetMaximum());
  tlb->SetLineWidth(2);
  tlb->Draw();
  TLine *tlu = new TLine(3.22,0,3.22,0.65*minv_pippim_dpm->GetMaximum());
  tlu->SetLineWidth(2);
  tlu->Draw();
  TLegend *tl_minv = new TLegend(0.5,0.65,0.9,0.9);
  tl_minv->AddEntry(minv_pippim_dpm,"M_{inv}^{#pi^{+}#pi^{-}} in DPM","l");
  tl_minv->AddEntry(minv_pippim_ftf,"M_{inv}^{#pi^{+}#pi^{-}} in FTF","l");
  tl_minv->AddEntry(tlu,"2#sigma M cut (J/#psi)","l");
  tl_minv->SetBorderSize(0);
  tl_minv->SetFillStyle(0);
  tl_minv->Draw();
  minv->Print("minv.pdf");

  return;

  TCanvas *the = new TCanvas("the","the",700,700);
  the->cd();
  the_pi0_cm_ftf->SetTitle("#bar{p}p#rightarrow#pi^{+}#pi^{-}#pi^{0} @ mom(#bar{p})=5.513 GeV, s=12.25 GeV^{2};#theta^{*}_{#pi^{0}};dN/d#theta^{*}_{#pi^{0}}");
  the_pi0_cm_ftf->Draw();
  the_pi0_cm_dpm->SetTitle("");
  the_pi0_cm_dpm->Draw("same");
  tl_the->Draw();
  minv->Print("the.eps");


}

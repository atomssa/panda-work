void set_style(TAxis *axis) {
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.05);
}

void t(){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  TFile *f_bef = TFile::Open("hists/sig_nofilt_noeff.root");
  TH2F *tb = (TH2F*) f_bef->Get("inv_p_t_pbar_pi0_lab_the_pi0");
  tb->SetTitle("#theta_{#pi^{0}} vs t;t[(GeV)^{2}];#theta_{#pi^{0}}[#circ]");
  set_style(tb->GetXaxis());
  set_style(tb->GetYaxis());


  TFile* f_aft = TFile::Open("hists/sig_noeff.root");
  TH1F *ta = (TH1F*) f_aft->Get("inv_p_t_pbar_pi0");
  //  tb->SetTitle(Form("#theta_{J/#psi} in #bar{p}p#rightarrow%s", rxn.c_str()));
  ta->SetTitle("t;t[(GeV)^{2}]");
  ta->GetXaxis()->SetRangeUser(-0.8,0.8);
  set_style(ta->GetXaxis());
  set_style(ta->GetYaxis());



  TCanvas *tc0 = new TCanvas("tc0","tc0",1100,500);
  tc0->Divide(2,1);
  tc0->cd(1);
  tb->Draw("colz");

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.06);
  tt->SetNDC();
  tt->DrawLatex(0.5,0.5,"s=12.25 GeV^{2}");

  tc0->cd(2);
  ta->Draw();
  tt->DrawLatex(0.19,0.7,"s=12.25 GeV^{2}");
  tt->DrawLatex(0.19,0.6,"After MC rejection");

  tc0->Print("figs/t_vs_the_pi0.pdf");

}

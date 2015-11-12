void purity() {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  //gROOT->LoadMacro(Form("%s/tda/tda.C",bdir));
  gROOT->LoadMacro(Form("%s/figs/v2/ananote.C",bdir));

  int ibrem = 1;
  int pass = 14;
  int iplab = 0;

  TFile *fsig = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_%s4eff_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));

  TH1F* hmep = (TH1F* ) fsig->Get("hmep_4");
  TH1F* hmep_mct = (TH1F* ) fsig->Get("mct/hmep_mct_4");
  hmep->SetTitle(";M_{e^{+}e^{-}}[GeV/c^{2}]");
  set_style(hmep,4);
  set_style(hmep_mct,2);

  TH1F* hmgg = (TH1F* ) fsig->Get("mct/hmgg_mct_4")->Clone("hmgg");
  hmgg->Add((TH1F* ) fsig->Get("mct/hmgg_non_mct_4"));
  TH1F* hmgg_mct = (TH1F* ) fsig->Get("mct/hmgg_mct_4");
  hmgg->SetTitle(";M_{#gamma#gamma}[GeV/c^{2}]");
  hmgg->GetXaxis()->SetRangeUser(0.1,0.164);
  set_style(hmgg,4);
  set_style(hmgg_mct,2);

  TLatex *tl = new TLatex();
  tl->SetTextSize(0.06);
  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(2);
  tc->cd(1);
  hmep->Draw();
  hmep_mct->Draw("same");
  gPad->SetLogy();
  tl->DrawLatexNDC(.15,.85, Form("purity=%4.1f%%", (100.0*hmep_mct->GetEntries())/hmep->GetEntries()));
  tc->cd(2);
  hmgg->Draw();
  hmgg_mct->Draw("same");
  //gPad->SetLogy();
  tl->DrawLatexNDC(.15,.85, Form("purity=%4.1f%%", (100.0*hmgg_mct->GetEntries())/hmgg->GetEntries()));

  tc->Print("purity_most_btob.pdf");

}

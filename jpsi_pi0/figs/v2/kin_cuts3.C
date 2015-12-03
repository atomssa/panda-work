void kin_cuts3() {

  //gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  gROOT->LoadMacro(Form("%s/figs/v2/ananote.C",bdir));

  TLegend *tl = new TLegend(0.17,0.65,0.45,0.85);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);

  TFile *fbrem = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_brem4eff_p0_pass10.root",bdir));
  TH1F *hmep_brem_bef = (TH1F*) fbrem->Get("hmep_5")->Clone("hmep_brem_bef");
  TH1F *hmep_brem_aft = (TH1F*) fbrem->Get("hmep_6")->Clone("hmep_brem_aft");
  set_style(hmep_brem_bef,2);
  set_style(hmep_brem_aft,4);
  tl->AddEntry(hmep_brem_bef,"Only EID cut","l");
  tl->AddEntry(hmep_brem_aft,"EID+#chi^{2}_{SIG} cut","l");

  TFile *fraw = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_raw4eff_p0_pass10.root",bdir));
  TH1F *hmep_raw_bef = (TH1F*) fraw->Get("hmep_5")->Clone("hmep_raw_bef");
  TH1F *hmep_raw_aft = (TH1F*) fraw->Get("hmep_6")->Clone("hmep_raw_aft");
  set_style(hmep_raw_bef,2);
  set_style(hmep_raw_aft,4);

  TGaxis::SetMaxDigits(3);
  TCanvas *tc =new TCanvas("tc","tc",1300,700);
  tc->Divide(2);
  tc->cd(2);
  hmep_brem_bef->SetTitle("Invariant mass, With Brem Corr;M_{ee}[GeV/c]");
  hmep_brem_bef->SetMaximum(30000);
  hmep_brem_bef->Draw();
  hmep_brem_aft->Draw("sames");
  tl->Draw();

  // Retrieve the stat box
  //TPaveStats *ps_brem_bef = (TPaveStats*) hmep_brem_bef->GetPrimitive("stats");
  //ps_brem_bef->SetName("TEST TEST");
  //ps_brem_bef->SetTextColor(2);
  //ps_brem_bef->SetLineColor(2);

  //gPad->SetLogy();
  tc->cd(1);
  hmep_raw_bef->SetTitle("Invariant mass, No Brem Corr;M_{ee}[GeV/c]");
  hmep_raw_bef->SetMaximum(30000);
  hmep_raw_bef->Draw();
  hmep_raw_aft->Draw("sames");
  tl->Draw();
  //gPad->SetLogy();

}

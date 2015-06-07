void set_style(TH1F* h) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  h->GetXaxis()->SetNdivisions(610);
}

void all(){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  //TFile *fsim= TFile::Open("new_test/anav2_jpsi_brem_plab5.5.root");
  //TFile *fbg= TFile::Open("new_test/anav2_pip_pim_brem_plab5.5.root");

  TFile *fsim= TFile::Open("new_test/anav2_jpsi_brem_plab8.0.root");
  TFile *fbg= TFile::Open("new_test/anav2_pip_pim_brem_plab12.0.root");

  TLegend *tl = new TLegend(0.35,0.5,0.85,0.9);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  tl->SetTextSize(0.06);

  TH1F* hmep[2][5];
  for (int i=0; i<5; ++i) {
    hmep[0][i] = (TH1F*)fsim->Get(Form("hmep_%d",i));
    set_style(hmep[0][i]);
    hmep[1][i] = (TH1F*)fbg->Get(Form("hmep_%d",i));
    set_style(hmep[1][i]);
  }

  TCanvas *tc_vert = new TCanvas("tc_vert","tc_vert",800,1200);
  tc_vert->Divide(1,2);
  for (int ii=0; ii < 2; ++ii) {
    tc_vert->cd(1+ii);
    hmep[ii][0]->Draw();
    hmep[ii][0]->SetTitle("");
    hmep[ii][0]->GetYaxis()->SetTitle("counts");
    hmep[ii][0]->SetLineColor(1);
    if (ii==0) tl->AddEntry(hmep[ii][0],"All pairs (5x10^{-7} for bg)");
    hmep[ii][1]->SetLineColor(2);
    hmep[ii][1]->Draw("same");
    if (ii==0) tl->AddEntry(hmep[ii][1],"Pairs after eid (wt for bg)");
    hmep[ii][3]->SetLineColor(4);
    hmep[ii][3]->Draw("same");
    if (ii==0) tl->AddEntry(hmep[ii][3],"N_{#pi^{0}}>0 after #pi^{0} cuts");
    hmep[ii][4]->SetLineColor(3);
    hmep[ii][4]->Draw("same");
    if (ii==0) tl->AddEntry(hmep[ii][4],"Most b-to-b #pi^{0}-(e^{+}e^{-}) pair");
    if (ii==1) tl->Draw();
    if (ii==1) {
      //gPad->SetLogy();
      //hmep[ii][0]->SetMinimum(10);
      hmep[ii][0]->Scale(0.5e-6);
    }
  }

  tc_vert->Print("figs/2015.06.08/ana_steps.pdf");

}

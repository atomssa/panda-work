void dphimax() {

  static const int nconfig = 22;
  Double_t config[nconfig][3] = {
    {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
    {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
    {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
    {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
    {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
    {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(6,4);

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);

  for (int iconfig = 0; iconfig < nconfig; ++iconfig) {
    TString fname = Form("../grid.out/psim_oct14_binsong_configs/bremcorr_merged.%d_hists.root",iconfig);
    cout << "Opening file " << fname << endl;
    TFile *f = TFile::Open(fname);
    TH1F* phi_max_rec = (TH1F*) f->Get("h_dphi_max_rec");
    TH1F* phi_max_mc = (TH1F*) f->Get("h_dphi_max_mc");
    tc->cd(1+iconfig);
    phi_max_mc->Draw();
    phi_max_rec->Draw("same");
    phi_max_rec->SetLineColor(2);
    phi_max_mc->SetMinimum(1);
    gPad->SetLogy();
    tt->SetTextColor(1);
    tt->SetTextSize(0.07);
    tt->DrawLatex(0.5,0.82,Form("p_{T}= %3.1f GeV/c",config[iconfig][0]));
    tt->DrawLatex(0.5,0.72,Form("%3.0f < #theta < %3.0f",config[iconfig][1],config[iconfig][2]));
  }

}

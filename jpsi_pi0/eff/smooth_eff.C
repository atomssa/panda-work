void smooth_eff()  {

  TFile *f = TFile::Open("eff/epem.root");

  TH2F *effh = (TH2F*) f->Get("heffelec");

  TCanvas *tc0 = new TCanvas("tc0","tc0");
  tc0->cd();
  effh->Draw("colz");

  const int nbins= 1000;

  TH2F *effhs_rad = new TH2F("eff_ep_em_rad", "Efficiency of e^{+} and e^{-}; p[GeV/c]; #theta[rad]", nbins, 0, 5, nbins, 0, TMath::Pi());
  TH2F *effhs_deg = new TH2F("eff_ep_em_deg", "Efficiency of e^{+} and e^{-}; p[GeV/c]; #theta[deg]", nbins, 0, 5, nbins, 0, 180);
  TH1F *effhsx = (TH1F*) effhs_deg->ProjectionX();
  TH1F *effhsy = (TH1F*) effhs_deg->ProjectionY();

  for (int i = 1; i <= nbins; i++) {
    for (int j = 0; j <= nbins; j++) {
      double xx = effhsx->GetBinCenter(i);
      double yy = effhsy->GetBinCenter(j);
      //cout << "xx= " << xx << " yy= " << yy << endl;
      double eff= effh->Interpolate(xx,yy);
      //double _xx = effh->GetXaxis()->FindBin(xx);
      //double _yy = effh->GetYaxis()->FindBin(yy);
      //if (eff!=0)
      //cout << "_xx= " << _xx << " _yy= " << _yy << " eff = " << eff << endl;
      // TODO -- how to get the error?
      effhs_deg->SetBinContent(i,j,eff);
      effhs_rad->SetBinContent(i,j,eff);
    }
  }

  TCanvas *tc1 = new TCanvas("tc1","tc1",720,10,700,700);
  tc1->cd();
  gStyle->SetOptStat(0);
  //effhs_deg->Draw("colz");
  effhs_rad->Draw("colz");
  tc1->Print("figs/epem_eff_mom_the.png");

  TFile *fout = TFile::Open("effic_smooth.root","RECREATE");
  fout->cd();
  effhs_rad->Write();
  effhs_deg->Write();
  fout->Write();
  fout->Close();

}

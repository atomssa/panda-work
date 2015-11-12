void num_evt() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  double sig[nplab] = {28040.,23619,15460};
  double bg[nplab] = {36000, 9000, 4500};

  TLatex *tl[6][nplab];
  TLegend *legend = new TLegend(0.3,0.6,0.99,0.8);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    for (int i=0; i<6; ++i) {
      tl[i][iplab] =new TLatex();
      //tl[i][iplab]->SetNDC(true);
      tl[i][iplab]->SetLineColor(i==0?2:1);
      //tl[i][iplab]->SetTextSize((i==1?1.1:1.4)*tl[i][iplab]->GetTextSize());
      //tl[i][iplab]->SetTextSize(((i==2)?1.4:1.5)*tl[i][iplab]->GetTextSize());
    }
  }

  TH1F* hdummy = new TH1F("my_dummy","my_dummy",10,0,15);
  set_style(hdummy,1);

  TGraph *tgsig = new TGraph(nplab, plab, sig);
  tgsig->SetTitle("Number of events in full signal MC;p_{lab}[GeV/c];num. evt");
  tgsig->SetMarkerStyle(20);
  tgsig->SetMarkerSize(2);
  tgsig->SetMarkerColor(2);
  tgsig->SetMarkerColor(2);
  TCanvas *tcsig = new TCanvas("num_evt_sig","num_evt_sig");
  tcsig->cd();
  hdummy->Draw();
  tgsig->Draw("ap");
  tgsig->SetMinimum(0);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tl[0][iplab]->DrawLatex(plab[iplab]+(iplab==2?-2.6:0.3),sig[iplab]*0.97,Form("%4.2f<|t/u|<%4.2f",tmin[iplab],tmax[iplab]));
  }
  tcsig->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir, tcsig->GetName()));

  TGraph *tgbg = new TGraph(nplab, plab, bg);
  tgbg->SetTitle("Number of events in full background MC;p_{lab}[GeV/c];num. evt");
  tgbg->SetMarkerStyle(24);
  tgbg->SetMarkerSize(2);
  tgbg->SetMarkerColor(2);
  TCanvas *tcbg = new TCanvas("num_evt_bg","num_evt_bg");
  tcbg->cd();
  hdummy->Draw();
  tgbg->Draw("ap");
  tgbg->SetMinimum(2);
  //gPad->SetLogy();
  //TGraph *tmg =  new TMultiGraph();
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tl[1][iplab]->DrawLatex(plab[iplab]+(iplab==2||iplab==1?-2.6:0.3),bg[iplab]*0.97,Form("%4.2f<|t/u|<%4.2f",tmin[iplab],tmax[iplab]));
  }
  tcbg->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir, tcbg->GetName()));

}

void draw_pad3(TCanvas *tc, int ic, const char* fname_mu, const char* fname_e, double pt, double thmin, double thmax) {

  int ipad  = ic+1;
  TFile *fmu = TFile::Open(fname_mu);
  //TH1F* rec_mu = ((TH1F*) fmu->Get("h_rec"))->Clone("h_rec_mu");
  //rec_mu->Sumw2();
  //rec_mu->SetLineWidth(1.5);
  //rec_mu->SetLineColor(1);
  //rec_mu->SetLineStyle(6);

  cout << "Opening file " << fname_e << endl;
  TFile *f = TFile::Open(fname_e);

  //TH1F* rec_e = ((TH1F*) f->Get("h_rec"))->Clone("h_rec_e");
  //rec_e->SetLineWidth(1);
  //rec_e->SetLineColor(2);
  //rec_e->SetLineStyle(6);
  //rec_e->Sumw2();

  TH1F* wbfcor_mw_bf_pc_e[nbs];
  for (int ibs=0; ibs < nbs; ++ibs) {
    wbfcor_mw_bf_pc_e[ibs] = (TH1F*)((TH1F*)f->Get(Form("b%d_h_wbfcor_mw_bf_pc",ibs)))->Clone(Form("b%d_h_wbfcor_mw_bf_pc_e",ibs));
  }

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);
  tt->SetTextSize(0.04);
  tc->cd(1+ic);
  //rec_mu->Scale(double(rec_e->Integral())/double(rec_mu->Integral()));
  //rec_mu->Scale(double(rec_e->GetEntries())/double(rec_mu->GetEntries()));
  //rec_e->SetTitle(Form("p_{T}= %4.2f GeV/c, %3.0f < #theta < %3.0f; (p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  //rec_mu->Draw("hist");
  //rec_e->Draw("hist");
  //int bs[5] = {0,4,5,6,7};

  int bs[5] = {1, 3, 0, 6, 7};
  wbfcor_mw_bf_pc_e[bs[2]]->SetTitle(Form("p_{T}= %4.2f GeV/c, %3.0f < #theta < %3.0f; (p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  wbfcor_mw_bf_pc_e[bs[2]]->Draw("hist");

  for (int ibs=0; ibs < 5; ++ibs) {
    //for (int ibs=0; ibs < 8; ++ibs) {
    int col = (ibs!=4)?(1+ibs):7;

    wbfcor_mw_bf_pc_e[bs[ibs]]->Draw("hist,same");
    wbfcor_mw_bf_pc_e[bs[ibs]]->SetLineColor(col);
    //int col = 1+ibs;
    //wbfcor_mw_bf_pc_e[ibs]->Draw("hist,same");
    //wbfcor_mw_bf_pc_e[ibs]->SetLineColor(col);
    tt->SetTextColor(col);
    tt->DrawLatex(0.159,0.82-(ibs*0.08),Form("Bin size: %4.3g#circ", epbs_binSize[bs[ibs]]));
  }
  tt->SetTextSize(0.08);
  tt->SetTextColor(1);
  bool fwd = thmin < 15.0 && 15.0<thmax;
  tt->DrawLatex(fwd?0.62:0.65,0.75,fwd?"Endcap":"Barrel");

  //gPad->SetLogy();
}

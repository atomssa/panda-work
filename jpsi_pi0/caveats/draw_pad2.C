void draw_pad2(int ibs, TCanvas *tc, int ic, const char* fname_mu, const char* fname_e, double pt, double thmin, double thmax) {
  int ipad  = ic+1;
  TFile *fmu = TFile::Open(fname_mu);
  TH1F* rec_mu = ((TH1F*) fmu->Get("h_rec"))->Clone("h_rec_mu");
  rec_mu->Sumw2();
  rec_mu->SetLineWidth(1.5);
  rec_mu->SetLineColor(1);
  rec_mu->SetLineStyle(6);

  cout << "Opening file " << fname_e << endl;
  TFile *f = TFile::Open(fname_e);

  TH1F* rec_e = ((TH1F*) f->Get("h_rec"))->Clone("h_rec_e");
  rec_e->SetLineWidth(1);
  rec_e->SetLineColor(2);
  rec_e->SetLineStyle(6);
  rec_e->Sumw2();

  TH1F* cor_e = ((TH1F*) f->Get(Form("b%d_h_cor",ibs)))->Clone(Form("b%d_h_cor_e",ibs));
  cor_e->SetLineWidth(2.0);
  cor_e->SetLineColor(3);
  cor_e->Sumw2();

  TH1F* wbfcor_e = ((TH1F*) f->Get(Form("b%d_h_wbfcor",ibs)))->Clone(Form("b%d_h_wbfcor_e",ibs));
  wbfcor_e->SetLineWidth(2.0);
  wbfcor_e->SetLineColor(2);
  wbfcor_e->Sumw2();

  TH1F* wbfcor_mw_bf_e = ((TH1F*) f->Get(Form("b%d_h_wbfcor_mw_bf",ibs)))->Clone(Form("b%d_h_wbfcor_mw_bf_e",ibs));
  wbfcor_mw_bf_e->SetLineWidth(2.0);
  wbfcor_mw_bf_e->SetLineColor(1);
  wbfcor_mw_bf_e->Sumw2();

  TH1F* wbfcor_mw_bf_pc_e = ((TH1F*) f->Get(Form("b%d_h_wbfcor_mw_bf_pc",ibs)))->Clone(Form("b%d_h_wbfcor_mw_bf_pc_e",ibs));
  wbfcor_mw_bf_pc_e->SetLineWidth(2.0);
  wbfcor_mw_bf_pc_e->SetLineColor(4);
  wbfcor_mw_bf_pc_e->Sumw2();

  tc->cd(1+ic);
  //rec_mu->Scale(double(rec_e->Integral())/double(rec_mu->Integral()));
  rec_mu->Scale(double(rec_e->GetEntries())/double(rec_mu->GetEntries()));
  rec_mu->SetTitle(Form("p_{T}= %4.2f GeV/c, %3.0f < #theta < %3.0f; (p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  rec_mu->Draw("hist");
  rec_e->Draw("hist,same");
  cor_e->Draw("hist,same");
  wbfcor_e->Draw("hist,same");
  wbfcor_mw_bf_e->Draw("hist,same");
  wbfcor_mw_bf_pc_e->Draw("hist,same");
  //gPad->SetLogy();

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);
  tt->SetTextSize(0.04);

  tt->SetTextColor(3);
  double chi2t_cor = rec_mu->Chi2Test(cor_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.82,Form("S+M, no wt, %5.2f%%", chi2t_cor));

  tt->SetTextColor(2);
  double chi2t_wbfcor = rec_mu->Chi2Test(wbfcor_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.72,Form("S_{wt}+M_{no wt}, %5.2f%%", chi2t_wbfcor));

  tt->SetTextColor(1);
  double chi2t_wbfcor_mw_bf = rec_mu->Chi2Test(wbfcor_mw_bf_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.62,Form("S_{wt}+M_{wt}, %5.2f%%", chi2t_wbfcor_mw_bf));

  tt->SetTextColor(4);
  double chi2t_wbfcor_mw_bf_pc = rec_mu->Chi2Test(wbfcor_mw_bf_pc_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.52,Form("S_{wt}+M_{wt,#Delta#phi cut}, %5.2f%%", chi2t_wbfcor_mw_bf_pc));
  tt->SetTextColor(4);

}

void draw_pad1(TCanvas *tc, int ic, const char* fname_mu, const char* fname_e, double pt, double thmin, double thmax) {
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
  rec_e->SetLineWidth(1.5);
  rec_e->SetLineColor(1);
  rec_e->Sumw2();

  TH1F* mrg_e = ((TH1F*) f->Get("h_mrg_w_bf_pc"))->Clone("h_mrg_e");
  mrg_e->SetLineWidth(1.5);
  mrg_e->SetLineColor(2);
  mrg_e->Sumw2();

  TH1F* sep_e = ((TH1F*) f->Get("h_sep_w_bf"))->Clone("h_sep_e");
  sep_e->SetLineWidth(1.5);
  sep_e->SetLineColor(4);
  sep_e->Sumw2();

  TH1F* cor_e = ((TH1F*) f->Get("h_wbfcor_mw_bf_pc"))->Clone("h_cor_e");
  cor_e->SetLineWidth(1.5);
  cor_e->SetLineColor(3);
  cor_e->Sumw2();

  tc->cd(1+ic);
  rec_mu->Scale(double(rec_e->GetEntries())/double(rec_mu->GetEntries()));
  rec_mu->SetTitle(Form("p_{T}= %4.2f GeV/c, %3.0f < #theta < %3.0f; (p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  rec_mu->Draw("hist");
  //cor_e->Draw();
  rec_e->Draw("hist,same");
  mrg_e->Draw("hist,same");
  sep_e->Draw("hist,same");
  cor_e->Draw("hist,same");
  //gPad->SetLogy();

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);
  tt->SetTextColor(1);
  tt->SetTextSize(0.05);
  //tt->DrawLatex(0.5,0.82,Form("p_{T}= %4.2f GeV/c",pt));
  //tt->DrawLatex(0.5,0.72,Form("%3.0f < #theta < %3.0f", thmin, thmax));

  tt->SetTextColor(1);
  double chi2t_raw = rec_mu->Chi2Test(rec_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.82,Form("Raw, %5.2f%%", chi2t_raw));

  tt->SetTextColor(2);
  double chi2t_mrg = rec_mu->Chi2Test(mrg_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.72,Form("Mrg, %5.2f%%", chi2t_mrg));

  tt->SetTextColor(4);
  double chi2t_sep = rec_mu->Chi2Test(sep_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.62,Form("Sep, %5.2f%%", chi2t_sep));

  tt->SetTextColor(3);
  double chi2t_cor = rec_mu->Chi2Test(cor_e,"UU,NORM,P");
  tt->DrawLatex(0.15,0.52,Form("M+S, %5.2f%%", chi2t_cor));

}

void draw_pad6(TCanvas *tc, int ic, const char* fname_mu, const char* fname_e, double pt, double thmin, double thmax) {
  int ipad  = ic+1;

  cout << "Opening file " << fname_e << endl;
  TFile *f = TFile::Open(fname_e);

  TH1F* rec_e = ((TH1F*) f->Get("h_rec"))->Clone("h_rec_e");
  rec_e->SetLineWidth(1.5);
  rec_e->SetLineColor(1);
  rec_e->Sumw2();

  TH1F* cor_e = ((TH1F*) f->Get("b0_h_cor"))->Clone("h_cor_e");
  cor_e->SetLineWidth(2);
  cor_e->SetLineColor(4);
  cor_e->Sumw2();

  TH1F* cor_w_e = ((TH1F*) f->Get("b0_h_wbfcor_mw_bf_pc"))->Clone("h_cor_w_e");
  cor_w_e->SetTitle(Form("; (p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  cor_w_e->SetLineWidth(2);
  cor_w_e->SetLineColor(2);
  cor_w_e->Sumw2();
  cor_w_e->GetXaxis()->SetNdivisions(605);

  tc->cd(1+ic);
  TGaxis::SetMaxDigits(3);
  cor_w_e->Draw("hist");
  cor_e->Draw("hist,same");
  rec_e->Draw("hist,same");

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);
  tt->SetTextColor(1);
  tt->SetTextSize(0.05);

  tt->SetTextColor(4);
  tt->DrawLatex(0.15,0.82,Form("KF+Unwtd. corr"));

  tt->SetTextColor(2);
  tt->DrawLatex(0.15,0.72,Form("KF+Wtd. corr"));

  tt->SetTextColor(1);
  tt->DrawLatex(0.15,0.62,Form("KF"));

  tt->SetTextColor(1);
  tt->DrawLatex(0.57,0.82,Form("p_{T}= %4.2g GeV/c",pt));
  tt->DrawLatex(0.57,0.72,Form("%3.0f < #theta < %3.0f",thmin,thmax));
  bool fwd = thmin < 15.0 && 15.0<thmax;
  tt->DrawLatex(fwd?0.60:0.63,0.62,fwd?"(Fwd Endcap)":"(Barrel)");

  gPad->SetGridx();

}

// for plots that go into analysis note
TPad* draw_pad5(int ic, TPad *tc, const char* dir, const char* ext,
	       double pt, double thmin, double thmax, TF1* rec, TF1* cor, double *tmp) {

  TCanvas *tctmp = new TCanvas("tctmp","tctmp");

  int ix= (ic%2);
  int iy= (ic/2);
  int ibs=0;

  TString fname_mu = Form("%s/mum%s",dir, ext);
  cout << "Opening file " << fname_mu  << endl;
  TFile *fmu = TFile::Open(fname_mu);
  TF1* ftmp = new TF1(Form("ftmp%d",ic),"gaus",-0.03,0.03);
  TH1F* rec_mu = ((TH1F*) fmu->Get("h_rec"))->Clone("h_rec_mu");
  rec_mu->Sumw2();
  rec_mu->SetLineWidth(2.0);
  rec_mu->SetLineColor(4);
  rec_mu->SetLineStyle(6);
  ftmp->SetLineColor(4);
  ftmp->SetLineWidth(0.5);

  TString fname_e = Form("%s/esim%s",dir, ext);
  cout << "Opening file " << fname_e  << endl;
  TFile *f = TFile::Open(fname_e);
  TH1F* rec_e = ((TH1F*) f->Get("h_rec"))->Clone("h_rec_e");
  rec_e->SetLineWidth(2);
  rec_e->SetLineColor(1);
  //  rec_e->SetLineStyle(6);
  rec_e->Sumw2();

  rec_mu->Scale(double(rec_e->GetEntries())/double(rec_mu->GetEntries()));
  rec_mu->Fit(ftmp,"+RIME");
  int min_mu = rec_mu->GetXaxis()->FindBin(ftmp->GetParameter(1)-2*ftmp->GetParameter(2));
  int max_mu = rec_mu->GetXaxis()->FindBin(ftmp->GetParameter(1)+2*ftmp->GetParameter(2));

  rec->SetLineColor(1);
  tctmp->cd();
  rec_e->Fit(rec,"+RIME");
  rec->SetParameter(1, rec->GetParameter(1)- ftmp->GetParameter(1));
  int min0 = rec_e->GetXaxis()->FindBin(rec->GetParameter(1)-2*rec->GetParameter(2));
  int max0 = rec_e->GetXaxis()->FindBin(rec->GetParameter(1)+2*rec->GetParameter(2));
  //double num0 = rec_e->Integral(min0, max0);
  double num0 = rec_e->Integral(min_mu, max_mu);
  double den0 = rec_e->GetEntries();
  tmp[0] = num0/den0;
  tmp[1] = TMath::Sqrt(tmp[0]*(1.0-tmp[0]))/TMath::Sqrt(den0);

  TH1F* cor_e = ((TH1F*) f->Get(Form("b%d_h_cor",ibs)))->Clone(Form("b%d_h_cor_e",ibs));
  cor_e->SetLineWidth(2.0);
  cor_e->SetLineColor(3);
  cor_e->Sumw2();

  TH1F* sep_w_bf = ((TH1F*) f->Get(Form("h_sep_w_bf")))->Clone(Form("h_sep_w_bf_e"));
  //TH1F* sep_w_bf = ((TH1F*) f->Get(Form("h_sep")))->Clone(Form("h_sep_w_bf_e"));
  sep_w_bf->SetLineWidth(2.0);
  sep_w_bf->SetLineColor(3);
  sep_w_bf->Sumw2();

  TH1F* wbfcor_mw_bf_pc_e = ((TH1F*) f->Get(Form("b%d_h_wbfcor_mw_bf_pc",ibs)))->Clone(Form("h_cor_e"));
  //TH1F* wbfcor_mw_bf_pc_e = ((TH1F*) f->Get(Form("b%d_h_cor",ibs)))->Clone(Form("b%d_h_wbfcor_mw_bf_pc_e",ibs));
  wbfcor_mw_bf_pc_e->SetLineWidth(2.0);
  wbfcor_mw_bf_pc_e->SetLineColor(2);
  wbfcor_mw_bf_pc_e->Sumw2();
  tctmp->cd();
  cor->SetLineColor(2);
  wbfcor_mw_bf_pc_e->Fit(cor,"+RIME");
  cor->SetParameter(1, cor->GetParameter(1)- ftmp->GetParameter(1));
  int min1 = wbfcor_mw_bf_pc_e->GetXaxis()->FindBin(cor->GetParameter(1)-2*cor->GetParameter(2));
  int max1 = wbfcor_mw_bf_pc_e->GetXaxis()->FindBin(cor->GetParameter(1)+2*cor->GetParameter(2));
  //double num1 = wbfcor_mw_bf_pc_e->Integral(min1, max1);
  double num1 = wbfcor_mw_bf_pc_e->Integral(min_mu, max_mu);
  double den1 = wbfcor_mw_bf_pc_e->GetEntries();
  tmp[2] = num1/den1;
  tmp[3] = TMath::Sqrt(tmp[2]*(1.0-tmp[2]))/TMath::Sqrt(den1);

  tctmp->Close();
  tc->cd();

  gPad->SetGridx();
  //gPad->SetLogy();

  rec_mu->SetTitle(Form(";(p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  rec_e->SetTitle(Form(";(p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  wbfcor_mw_bf_pc_e->SetTitle(Form(";(p_{MC} - p)/p_{MC}",pt,thmin,thmax));
  rec_mu->GetXaxis()->SetNdivisions(605);
  rec_e->GetXaxis()->SetNdivisions(605);
  wbfcor_mw_bf_pc_e->GetXaxis()->SetNdivisions(605);

  rec_mu->GetXaxis()->SetRangeUser(-0.1,0.2);

  rec_mu->Draw("hist");
  //rec_mu->Draw("func,same");

  rec_e->Draw("hist,same");
  //rec_e->Draw("func,same");
  sep_w_bf->Draw("hist,same");

  wbfcor_mw_bf_pc_e->Draw("hist,same");
  //wbfcor_mw_bf_pc_e->Draw("func,same");

  //return;

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);

  //tt->SetTextSize((ix==0||ix==2)?0.068:0.078);
  //double offx = ix==0?0.1:0.0;
  //double offy = iy==1?0.1:0.05;

  tt->SetTextSize(0.068);
  double offx = 0.15;
  double offy = 0.0;

  //TPaveText *tpt = new TPaveText(offx+0.43,offy+0.22,0.88,0.93,"NDC");
  TPaveText *tpt;

  if (ix==0&&iy==0)
    tpt = new TPaveText(offx+0.43,offy+0.22,0.88,0.85,"NDC");
  else
    tpt = new TPaveText(offx+0.43,offy+0.55,0.88,0.9,"NDC");
  tpt->SetFillStyle(1001);
  tpt->SetFillColor(0);
  tpt->SetBorderSize(0);
  tpt->Draw();

  if (ix==0&&iy==0) {
    tt->SetTextColor(4);
    tt->DrawLatex(offx+0.55,offy+0.45,Form("#mu^{-}, KF"));
    tt->SetTextColor(1);
    tt->DrawLatex(offx+0.55,offy+0.35,Form("e^{-}, KF"));
    tt->SetTextColor(3);
    tt->DrawLatex(offx+0.4,offy+0.25,Form("e^{-}, KF+corr. (sep. only)"));
    tt->SetTextColor(2);
    tt->DrawLatex(offx+0.5,offy+0.15,Form("e^{-}, KF+corr."));
  }

  tt->SetTextColor(1);
  tt->DrawLatex(offx+0.45,offy+0.82,Form("p_{T}= %4.2g GeV/c",pt));
  tt->DrawLatex(offx+0.45,offy+0.72,Form("%3.0f#circ < #theta < %3.0f#circ",thmin,thmax));
  bool fwd = thmin < 15.0 && 15.0<thmax;
  tt->DrawLatex(fwd?offx+0.45:offx+0.5,offy+0.62,fwd?"(Fwd Endcap)":"(Barrel)");

  return tc;
}

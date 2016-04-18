void set_style(TH1* h, int col, int rebin, bool sumw2) {
  if (rebin>0)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void jpsi(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  TFile *fb = TFile::Open(Form("%s/hists/note.aug.2015/eid90pct/pass3/anav2_jpsi_brem_plab5.5.root",bdir));
  TFile *fr = TFile::Open(Form("%s/hists/note.aug.2015/eid90pct/pass3/anav2_jpsi_raw_plab5.5.root",bdir));
  //TFile *fb = TFile::Open(Form("%s/pi.brem.root",bdir));
  //TFile *fr = TFile::Open(Form("%s/pi.raw.root",bdir));

  TH1F* jpsib = (TH1F*) fb->Get("hmep_3")->Clone("jpsib");
  //TH1F* jpsib = (TH1F*) fb->Get("hmep_3")->Clone("jpsib");
  set_style(jpsib, 4, 0, false);
  jpsib->GetXaxis()->SetRangeUser(0,4.5);
  jpsib->SetTitle(";M_{e^{+}e^{-}}[GeV/c^{2}]; counts");
  TH1F* jpsir = (TH1F*) fr->Get("hmep_3");
  //TH1F* jpsir = (TH1F*) fr->Get("hmep_3");
  set_style(jpsir, 2, 0, false);
  //jpsir->SetLineStyle(8);

  TCanvas *tc = new TCanvas("tc","tc",700,700);
  //TCanvas *tc = new TCanvas("tc","tc");
  tc->cd();
  jpsib->Draw();
  jpsir->Draw("sames");

  TLegend *tl = new TLegend(0.16,0.6,0.44,0.87);
  //TLegend *tl = new TLegend(0.42,0.6,0.7,0.87);
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);
  tl->SetTextSize(0.06);
  tl->AddEntry(jpsib, "With Brem Corr.", "pl");
  tl->AddEntry(jpsir, "No Brem Corr.", "pl");
  tl->Draw();

  //tc->Print(Form("%s/figs/2015.09.15/jpsi_mass_brem_vs_raw.pdf",bdir));
  //tc->Print(Form("%s/figs/2015.09.15/pipm_mass_brem_vs_raw.pdf",bdir));

  TCanvas *tcfit = new TCanvas("tcfit","tcfit",700,0,700,700);
  tcfit->cd();
  TH1F* jpsibf = (TH1F*) fb->Get("hmep_3")->Clone("jpsibf");
  set_style(jpsibf, 4, 0, false);
  jpsibf->GetXaxis()->SetRangeUser(2.5,4.0);
  jpsibf->GetXaxis()->SetNdivisions(505);
  jpsibf->Sumw2();
  jpsibf->Sumw2();
  jpsibf->SetTitle(";M_{e^{+}e^{-}}[GeV/c^{2}]; counts");
  TF1 *fgaus=new TF1("fgaus","gaus",3.0,3.2);
  //TF1 *fgaus=new TF1("fgaus","gaus",2.5,3.5);
  fgaus->SetParName(2,"width");
  fgaus->SetParName(1,"center");
  fgaus->SetParameter(0, 100);
  fgaus->SetParameter(1, 3.1);
  fgaus->SetParameter(1, 0.01);
  jpsibf->Fit(fgaus,"RE");

  double lb3s = fgaus->GetParameter(1)-3*fgaus->GetParameter(2);
  double ub3s = fgaus->GetParameter(1)+3*fgaus->GetParameter(2);

  TLine *ll3s = new TLine();
  ll3s->SetLineWidth(2);
  //ll3s->SetLineColor(3);
  ll3s->SetLineStyle(7);
  ll3s->DrawLine(lb3s, 0, lb3s, 200);
  TLine *lu3s = new TLine();
  lu3s->SetLineWidth(2);
  //lu3s->SetLineColor(3);
  lu3s->SetLineStyle(7);
  lu3s->DrawLine(ub3s, 0, ub3s, 200);

  TLine *ll = new TLine();
  ll->SetLineWidth(2);
  ll->DrawLine(2.8, 0, 2.8, 200);
  TLine *lu = new TLine();
  lu->SetLineWidth(2);
  lu->DrawLine(3.3, 0, 3.3, 200);

  TLegend *tlf = new TLegend(0.55,0.3,0.9,0.5);
  tlf->SetFillStyle(0);
  tlf->SetBorderSize(0);
  tlf->SetTextSize(0.05);
  //tlf->AddEntry(jpsibf, "J/#psi mass", "pl");
  tlf->AddEntry(fgaus, "Gaussian fit", "pl");
  tlf->Draw();
  tlf->AddEntry(lu3s, "3#sigma window ", "l");
  tlf->AddEntry(lu, "Analysis cut", "l");
  tlf->Draw();

  mean= fgaus->GetParameter(1);
  width= fgaus->GetParameter(2);

  cout <<"mean= " << mean << endl;
  cout <<"width= " << width << endl;
  cout << "min= " << mean-3*width << endl;
  cout << "max= " << mean+3*width << endl;
  //tcfit->Print(Form("%s/figs/2015.09.15/jpsi_mass_fit.pdf",bdir));
  //tcfit->Print(Form("%s/figs/2015.09.15/pipm_mass_fit.pdf",bdir));

}

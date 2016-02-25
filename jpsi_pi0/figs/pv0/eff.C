static const double Mp = 0.938;
static const double Mp_sq = Mp*Mp;
static const double Mpi0 = 0.0; //0.135;
static const double Mpi0_sq = Mpi0*Mpi0;
static const double Mjpsi = 3.097;
static const double Mjpsi2 = Mjpsi*Mjpsi;
static const double msqTot = 2*Mp_sq + Mjpsi2 + Mpi0_sq;

double mirror(double t, double s) {
  return msqTot - s - t;
}

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

void get_vector( vector<double> &v, TVectorD* vv) {
  cout << "FUCKING ROOTERIE. GO FIGURE! Getting vector<double> from TVectorD" << endl;
  if (vv!=0) {
    vv->Print();
    int nelt = vv->GetNoElements();
    cout << "NoElements= " << nelt << endl;
    for (int ielt=0; ielt < nelt; ++ielt) {
      TVectorD vvv = vv[0];
      double elt = vvv[ielt];
      cout << "elt." << ielt << " = " << elt << endl;
      v.push_back(elt);
    }
  } else {
    cout << "vector object is null" <<endl;
  }
}

double integrate_content(TH1F* h, double min, double max) {
  int i0 = h->GetXaxis()->FindBin(min);
  int i1 = h->GetXaxis()->FindBin(max);
  return h->Integral(i0, i1);
}

void eff(){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  const int nplab = 3;
  const double plab[nplab] = {5.513, 8., 12.};
  const double _s[nplab] = {12.25, 16.87, ,24.35};

  const double tvalidmin[nplab] = {-0.092, -1.0, -1.0};
  const double tvalidmax[nplab] = {0.59, 0.43, 0.3};

  TFile *f[nplab];
  TH1F* total[nplab];
  TH1F* passed[nplab];
  TEfficiency *eff[nplab];
  TGraphAsymmErrors* tge[nplab];
  int col[nplab] = {1,2,4};

  TLegend *legend3 = new TLegend(0.15, 0.6, 0.4, 0.9);
  legend3->SetFillStyle(0);
  legend3->SetBorderSize(0);
  legend3->SetTextSize(0.06);

  TMultiGraph *tmg = new TMultiGraph("tmg",";t[GeV^{2}];Signal Reco #varepsilon_{tot}");
  TCanvas *tceff = new TCanvas("tceff","tceff");
  tceff->cd();
  int nrebin = 10;
  for (int iplab = nplab-1; iplab >=0; --iplab) {
    //f[iplab] = TFile::Open(Form("%s/hists/note.aug.2015/eid90pct/anav2_pi0jpsi_tda_eff_ip%d_brem.root",bdir,iplab));
    int ibrem = 1;

    //int pass = 14;
    //f[iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_%s4eff_p%d_pass%d.root", bdir, (ibrem==0?"raw":"brem"), iplab, pass));

    int pass = 18;
    f[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_%s4eff_p%d_pass%d.root", bdir, (ibrem==0?"raw":"brem"), iplab, pass));

    total[iplab] = (TH1F*) f[iplab]->Get("tu/httrumc");
    passed[iplab] = (TH1F*) f[iplab]->Get("tu/htrecgg");
    total[iplab]->Rebin(nrebin);
    passed[iplab]->Rebin(nrebin);

    //total[iplab] = (TH1F*) f[iplab]->Get("tu/htrupi0thlab");
    //passed[iplab] = (TH1F*) f[iplab]->Get("pi0th_bins/hpi0th");
    //int r = total[iplab]->GetXaxis()->GetNbins()/passed[iplab]->GetXaxis()->GetNbins();
    //total[iplab]->Rebin(r*nrebin);
    //passed[iplab]->Rebin(nrebin);

    //total[iplab] = (TH1F*) f[iplab]->Get("tu/htrupi0costhcm");
    //passed[iplab] = (TH1F*) f[iplab]->Get("pi0cost_cm_bins/hpi0cost_cm");
    //int r = total[iplab]->GetXaxis()->GetNbins()/passed[iplab]->GetXaxis()->GetNbins();
    //total[iplab]->Rebin(r*nrebin);
    //passed[iplab]->Rebin(nrebin);

    int nbins = total[iplab]->GetXaxis()->GetNbins();
    double eff_avg = 0.0;
    double eff_cnt = 0.0;
    for (int i=0; i<=nbins+1; ++i) {
      double _t = passed[iplab]->GetXaxis()->GetBinCenter(i);
      double _u = mirror(_t,_s[iplab]);
      bool t_ok = (tvalidmin[iplab] < _t && _t < tvalidmax[iplab]);
      bool u_ok = (tvalidmin[iplab] < _u && _u < tvalidmax[iplab]);
      bool _valid = t_ok || u_ok;

      if (passed[iplab]->GetBinContent(i)>0.25*total[iplab]->GetBinContent(i)) {
      	passed[iplab]->SetBinContent(i,0.0);
      }
      if (_valid && total[iplab]->GetBinContent(i)>0) {
	eff_avg += passed[iplab]->GetBinContent(i)/total[iplab]->GetBinContent(i);
	eff_cnt+= 1.0;
      }
    }

    double _tmid = (tvalidmax[iplab] + mirror(tvalidmax[iplab],_s[iplab]))/2.0;
    cout << "p= " << plab[iplab] << " tmid= " <<  _tmid << " Average efficiency = " << eff_avg/eff_cnt << endl;

    eff[iplab] = new TEfficiency(*passed[iplab], *total[iplab]);
    eff[iplab]->SetLineColor(col[iplab]);
    eff[iplab]->Draw(iplab==nplab-1?"alp":"alp");
    gPad->Update();
    tge[iplab] = (TGraphAsymmErrors*) eff[iplab]->GetPaintedGraph();

    tmg->Add(tge[iplab], "lp");

  }

  legend3->AddEntry(tge[0], Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[0]), "pl");
  legend3->AddEntry(tge[1], Form("p^{LAB}_{#bar{p}} = %5.3g GeV/c",plab[1]), "pl");
  legend3->AddEntry(tge[2], Form("p^{LAB}_{#bar{p}} = %5.3g GeV/c",plab[2]), "pl");

  tmg->Draw("a");
  legend3->Draw();

  TLine *tlmin[2][nplab], *tlmax[2][nplab];
  TArrow *range[2][nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {

    int _iplab = 2 - iplab;
    double vst = (0.03*_iplab+0.01)/2;
    double vend = (0.03*(_iplab+1))/2;
    double vmid = (0.03*(_iplab+1)-0.01)/2;

    tlmin[0][iplab] = new TLine(tvalidmin[iplab],vst,tvalidmin[iplab],vend);
    tlmin[0][iplab]->SetLineColor(col[iplab]);
    tlmin[0][iplab]->SetLineWidth(2);
    tlmin[0][iplab]->Draw();
    tlmax[0][iplab] = new TLine(tvalidmax[iplab],vst,tvalidmax[iplab],vend);
    tlmax[0][iplab]->SetLineColor(col[iplab]);
    tlmax[0][iplab]->SetLineWidth(2);
    tlmax[0][iplab]->Draw();
    range[0][iplab] = new TArrow(tvalidmin[iplab], vmid, tvalidmax[iplab], vmid, 0.015, "<>");
    range[0][iplab]->SetLineColor(col[iplab]);
    range[0][iplab]->SetLineStyle(7);
    range[0][iplab]->SetLineWidth(2);
    range[0][iplab]->Draw();

    tlmin[1][iplab] = new TLine(mirror(tvalidmax[iplab],_s[iplab]),vst,mirror(tvalidmax[iplab],_s[iplab]),vend);
    tlmin[1][iplab]->SetLineColor(col[iplab]);
    tlmin[1][iplab]->SetLineWidth(2);
    tlmin[1][iplab]->Draw();
    tlmax[1][iplab] = new TLine(mirror(tvalidmin[iplab],_s[iplab]),vst,mirror(tvalidmin[iplab],_s[iplab]),vend);
    tlmax[1][iplab]->SetLineColor(col[iplab]);
    tlmax[1][iplab]->SetLineWidth(2);
    tlmax[1][iplab]->Draw();
    range[1][iplab] = new TArrow(mirror(tvalidmax[iplab],_s[iplab]), vmid, mirror(tvalidmin[iplab],_s[iplab]), vmid, 0.015, "<>");
    range[1][iplab]->SetLineColor(col[iplab]);
    //range[1][iplab]->SetLineStyle(8);
    range[1][iplab]->SetLineWidth(2);
    range[1][iplab]->Draw();

  }

  TH1F* dum = (TH1F*) tmg->GetHistogram();
  dum->GetYaxis()->SetNdivisions(508);
  set_style(dum,4,0,false);
  gPad->Update();

  //tceff->Print(Form("%s/figs/v2/efficiency_vs_t.pdf",bdir));
  tceff->Print(Form("efficiency_vs_t.pdf",bdir));

  /*
  vector<double> tu_bins;
  TVectorD *vv =  (TVectorT<double>*)f->Get("tu_binning");
  get_vector(tu_bins,vv);
  cout << "tu_binningsize= " << tu_bins.size() << endl;
  const int ntbin = tu_bins.size()>0?tu_bins.size()-1:12;
  const int ntbin_max = 20;
  double eff_cor[2][ntbin_max]={{0.0}};
  double eff_cor_er[2][ntbin_max]={{0.0}};
  double t[ntbin_max] = {0.0}, t_er[ntbin_max] = {0.};
  TH1F* h_teff_den = (TH1F*) f->Get("tu/httrumc");
  TH1F* h_ueff_den = (TH1F*) f->Get("tu/hutrumc");
  TH1F* h_teff_num = (TH1F*) f->Get("tu/htrecgg");
  TH1F* h_ueff_num = (TH1F*) f->Get("tu/hurecgg");

  for (int itbin=0; itbin < ntbin; ++itbin) {
    t[itbin] = (tu_bins[itbin+1]+tu_bins[itbin])/2.0;
    eff_cor[0][itbin] = integrate_content(h_teff_num, tu_bins[itbin], tu_bins[itbin+1]);
    double d0 = integrate_content(h_teff_den, tu_bins[itbin], tu_bins[itbin+1]);
    eff_cor[0][itbin] /= d0;
    double tmp0 = eff_cor[0][itbin];
    eff_cor_er[0][itbin] = TMath::Sqrt(tmp0*(1-tmp0))/TMath::Sqrt(d0);

    eff_cor[1][itbin] = integrate_content(h_ueff_num, tu_bins[itbin], tu_bins[itbin+1]);
    double d1 = integrate_content(h_ueff_den, tu_bins[itbin], tu_bins[itbin+1]);
    eff_cor[1][itbin] /= d1;
    double tmp1 = eff_cor[1][itbin];
    eff_cor_er[1][itbin] = TMath::Sqrt(tmp1*(1-tmp1))/TMath::Sqrt(d1);

    cout << "t= " << t[itbin] << " eff= " << eff_cor[0][itbin] << endl;
  }

  TGraphErrors *tget = new TGraphErrors(ntbin-1, t, eff_cor[0], t_er, eff_cor_er[0]);
  tget->SetMarkerStyle(20);
  tget->SetMarkerColor(2);
  TGraphErrors *tgeu = new TGraphErrors(ntbin-1, t, eff_cor[1], t_er, eff_cor_er[1]);

  TCanvas *tceff2 = new TCanvas("tceff2","tceff2");
  tceff2->cd();
  tget->Draw("ap");
  tget->SetMinimum(0);
  eff->Draw("same");
  */

}

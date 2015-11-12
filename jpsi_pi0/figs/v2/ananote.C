#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1F.h"

static int ibrem = 1;
static const int nplab = 3;
static const int nbg = 4;
static const double plab[nplab] = {5.513, 8., 12.};
// these are numbers derived from analytical formula for
// values of t at cos_theta_cm = {-1, 0, 1}
static double limits[nplab][3] =
  { {-1.50406, -0.443789, 0.616486},
    {-5.96462, -2.75368, 0.457248},
    {-13.2926, -6.48861, 0.31538} };
static double tmin[nplab]={-0.443789, -0.5, -0.5};
static double tmax[nplab]={0.616486, 0.457248, 0.31538};

void set_style(TH2* h) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);

  //h->GetXaxis()->SetTitleSize(0.08);
  //h->GetXaxis()->SetLabelSize(0.08);
  //h->GetXaxis()->SetLabelOffset(0.);
  //h->GetXaxis()->SetTitleOffset(0.8);
  //h->GetYaxis()->SetTitleSize(0.08);
  //h->GetYaxis()->SetLabelSize(0.08);
  //h->GetYaxis()->SetLabelOffset(0.02);
  //h->GetYaxis()->SetTitleOffset(0.8);

}

void set_style(TH1* h, int col, int rebin=0, bool sumw2=false) {
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

void set_style_ana(TH1* h, int col, int rebin, bool sumw2) {
  if (rebin>1)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.08);
  h->GetXaxis()->SetLabelSize(0.08);
  h->GetXaxis()->SetLabelOffset(0.);
  h->GetXaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.08);
  h->GetYaxis()->SetLabelSize(0.08);
  h->GetYaxis()->SetLabelOffset(0.02);
  h->GetYaxis()->SetTitleOffset(0.8);

  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void set_style(TGraph *g, int col, int sty, int s, int w) {
  g->SetMarkerStyle(sty);
  g->SetMarkerSize(s);
  g->SetMarkerColor(col);
  g->SetLineColor(col);
  g->SetLineWidth(w);
}

void draw_fg_sig_bg(TCanvas *canv, int ipad, TH1F*fg, TH1F*sig, TH1F*bg, bool legend=false) {
  canv->cd(ipad);
  fg->Draw("hist");
  sig->Draw("same");
  bg->Draw("same");
  if (legend) {
    TLegend *tl = new TLegend(0.2,0.6, 0.45,0.88);
    tl->AddEntry(fg,"Foreground","l");
    tl->AddEntry(sig,"Signal","pl");
    tl->AddEntry(bg,"Background","pl");
    tl->Draw();
  }
}

void fetch_fg_bg_sig(TFile*fsig, TFile *fbg, const char*dir, const char *stag, const char*dtag, TH1F* fg, TH1F* sig, TH1F*bg){
  fg = (TH1F*) fsig->Get(Form("%shmep%s",dir,stag))->Clone(Form("hmep_fg_%s",dtag));
  sig = (TH1F*) fsig->Get(Form("%shmep%s",dir,stag))->Clone(Form("hmep_sig_%s",dtag));
  bg = (TH1F*) fbg->Get(Form("%shmep%s",dir,stag))->Clone(Form("hmep_bg_%s",dtag));
  fg->Add(bg);
  set_style(fg, 1);
  set_style(sig, 2);
  set_style(bg, 4);
}

// ana2

Double_t background2(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t background3(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

const double sqrt_2pi = TMath::Sqrt(2*TMath::Pi());
double binw = 0.0;
Double_t gaussianPeak(Double_t *x, Double_t *par) {
  double exp = (x[0]-par[1])/par[2];
  return binw*par[0] * TMath::Exp(-exp*exp/2.)/par[2]/sqrt_2pi;
}

Double_t fitFunctionPol2(Double_t *x, Double_t *par) {
  //return background(x,par) + langaufun(x,&par[3]);
  return background2(x,par) + gaussianPeak(x,&par[3]);
}

Double_t fitFunctionPol3(Double_t *x, Double_t *par) {
  //return background(x,par) + langaufun(x,&par[3]);
  return background3(x,par) + gaussianPeak(x,&par[4]);
}

void get_vector( vector<double> &v, TVectorD* vv) {
  if (vv!=0) {
    //vv->Print();
    int nelt = vv->GetNoElements();
    //cout << "NoElements= " << nelt << endl;
    for (int ielt=0; ielt < nelt; ++ielt) {
      TVectorD vvv = vv[0];
      double elt = vvv[ielt];
      //cout << "elt." << ielt << " = " << elt << endl;
      v.push_back(elt);
    }
  } else {
    //cout << "vector object is null" <<endl;
  }
}

/// ana3

double integrate_content(TH1F* h, double min, double max) {
  int i0 = h->GetXaxis()->FindBin(min);
  int i1 = h->GetXaxis()->FindBin(max);
  return h->Integral(i0, i1);
}

double get_mean(TH1F* h, double min, double max) {
  int i0 = h->GetXaxis()->FindBin(min);
  int i1 = h->GetXaxis()->FindBin(max);
  double num = 0;
  double den = 0;
  for (int i=i0; i<=i1; ++i) {
    num += h->GetXaxis()->GetBinCenter(i)*h->GetBinContent(i);
    den += h->GetBinContent(i);
  }
  return num/den;
}

double calc_err_r(double n, double d, double ne, double de) {
  return n*TMath::Hypot(ne/n, de/d)/d;
}

///// tvsthcm
static const double mpi = 0.135;
static const double mj = 3.096;
static const double mp = 0.938;
double s(double _plab) {
  return 2.*mp*mp + 2.*mp*TMath::Hypot(mp,_plab);
}
double pcm(double _plab) {
  return _plab*mp/TMath::Sqrt(s(_plab));
}
double q(double _plab) {
  return TMath::Sqrt( pow(s(_plab),2) + pow(mj*mj-mpi*mpi,2) - (2*s(_plab)*(mj*mj+mpi*mpi)))/2./TMath::Sqrt(s(_plab));
}
double costh(double t, double _plab) {
  return ( t + pow(pcm(_plab),2) + pow(q(_plab),2) - pow(TMath::Hypot(q(_plab),mpi) - TMath::Hypot(pcm(_plab),mp), 2) )/2./q(_plab)/pcm(_plab);
}

void figs(int _page) {
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  switch(_page) {
  case 0: fig_tvsthcm_sig(); break;
  case 1: fig_bg_xsect(); break;
  case 2: fig_ana(); break;
  case 3: fig_pi0cut(); break;
  case 4: fig_npair(); break;
  case 5: fig_gen_dists(); break;
  case 6: fig_kin_cuts(); break;
  case 7: fig_num_evt(); break;
  case 8: fig_ana2(); break;
  case 9: fig_ana3(); break;
  default: return;
  }
}

#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TText.h"
#include "TLatex.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"

#include "Riostream.h"
#include <string>
#include <iostream>

TH1* def_hist(const char *name, const int &nbin, const double &bmin, const double &bmax, const int &col) {
  cout << "Hist def: " << name << " bmin= " << bmin << " bmax= "<< bmax << endl;
  TH1F* h = new TH1F(name, name, nbin, bmin, bmax);
  h->SetLineColor(col);
  h->SetLineWidth(2);
  return h;
  cout << "Done Hist def: " << name << endl;
}

void style_tgraph(TGraph *tg, const char *name, const char *title, int col, int mar_sty, double mar_size) {
  tg->SetTitle(title);
  tg->SetName(name);
  tg->SetMarkerStyle(mar_sty);
  tg->SetMarkerSize(mar_size);
  tg->SetMarkerColor(col);
  tg->SetLineColor(col);
}

TF1* fit_xy_res(double &x, double &xe, double &max, bool do_fit, int idel, int ipad, TCanvas *tc, TH1F* h, const char *meth, const char *var_tag, const char *det_tag, const char *var, const char *det_desc) {

  tc->cd(1+ipad);

  cout << "var_tag = " << var_tag << endl;
  TLatex *tt = new TLatex();
  tt->SetTextSize(0.08);
  tt->SetTextColor(4);
  tt->SetNDC(true);

  double fmx = 2.2;
  double fit_min= h->GetMean()-fmx*h->GetRMS();
  double fit_max= h->GetMean()+fmx*h->GetRMS();
  TF1* f = new TF1(Form("fxres_%s_%d_%s",det_tag,idel,meth), "gaus",fit_min,fit_max);

  if (do_fit) {

    h->Fit(f,"R+Q");
    x = f->GetParameter(2);
    xe = f->GetParError(2);

    tt->DrawLatex(0.65,0.8, meth);
    tt->DrawLatex(0.65,0.7, var);
    tt->DrawLatex(0.65,0.6, det_desc);
    tt->DrawLatex(0.65,0.5, Form("#delta = %d%%",idel-5 ));

  } else {

    max = TMath::Max(max, x);
    h->Draw();
    tt->DrawLatex(0.12,0.15, h->GetTitle());

  }

  max = TMath::Max(max, x);

  return f;
}

static const int ndel =  11;
static const int nmeth = 2; // 0->TDR, 1->Minuit
static const int nmeth_max = 2;

// pos resolution distribution maximia
static Double_t xmax_dia[nmeth][ndel] = {{0.0}};
static Double_t ymax_dia[nmeth][ndel] = {{0.0}};
static Double_t xmax_had[nmeth][ndel] = {{0.0}};
static Double_t ymax_had[nmeth][ndel] = {{0.0}};

static TH1F *hpres[nmeth][ndel];
static TH1F *hxres_dia[nmeth][ndel];
static TH1F *hyres_dia[nmeth][ndel];
static TH1F *hxres_had[nmeth][ndel];
static TH1F *hyres_had[nmeth][ndel];

static TF1 *fxres_dia[nmeth][ndel];
static TF1 *fyres_dia[nmeth][ndel];
static TF1 *fxres_had[nmeth][ndel];
static TF1 *fyres_had[nmeth][ndel];

void get_hist_limits(const string &fName) {
  cout << "Fname = " << fName <<endl;
  ifstream inf;
  inf.open( Form("%s.bins", fName.c_str()) );
  int dummy = 0;
  inf >> dummy;
  if (inf.good())  {
    for (int idel=0; idel<ndel; ++idel)
      for (int imeth=0; imeth<nmeth_max; ++imeth)
	inf >> xmax_dia[imeth][idel] >> ymax_dia[imeth][idel] >> xmax_had[imeth][idel] >> ymax_had[imeth][idel];
  } else {
    for (int idel=0; idel<ndel; ++idel) {
      for (int imeth=0; imeth<nmeth_max; ++imeth) {
	xmax_dia[imeth][idel] = 25;  ymax_dia[imeth][idel] = 25;  xmax_had[imeth][idel] = 25;  ymax_had[imeth][idel] = 25;
      }
    }
  }
  inf.close();
}

void save_hist_limits(const string &fName) {
  ofstream outf;
  outf.open( Form("%s.bins", fName.c_str() ) );
  cout << Form("%s.bins", fName.c_str() ) << endl;
  int dummy = 0;
  outf << dummy << "  " << endl;
  for (int idel=0; idel<ndel; ++idel) {
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      cout << fabs(fxres_dia[imeth][idel]->GetParameter(1) - 3*fxres_dia[imeth][idel]->GetParameter(2)) << "  "
  	   << fabs(fyres_dia[imeth][idel]->GetParameter(1) - 3*fyres_dia[imeth][idel]->GetParameter(2)) << "  "
  	   << fabs(fxres_had[imeth][idel]->GetParameter(1) - 3*fxres_had[imeth][idel]->GetParameter(2)) << "  "
  	   << fabs(fyres_had[imeth][idel]->GetParameter(1) - 3*fyres_had[imeth][idel]->GetParameter(2)) << endl;
      outf << fabs(fxres_dia[imeth][idel]->GetParameter(1) - 3*fxres_dia[imeth][idel]->GetParameter(2)) << "  "
  	   << fabs(fyres_dia[imeth][idel]->GetParameter(1) - 3*fyres_dia[imeth][idel]->GetParameter(2)) << "  "
  	   << fabs(fxres_had[imeth][idel]->GetParameter(1) - 3*fxres_had[imeth][idel]->GetParameter(2)) << "  "
  	   << fabs(fyres_had[imeth][idel]->GetParameter(1) - 3*fyres_had[imeth][idel]->GetParameter(2)) << endl;
    }
  }
  cout << endl;
  outf << endl;
  outf.close();
}

double get_max_val_and_error(TGraphErrors *tge, double &max, double &max_e) {
  double *val = tge->GetY();
  double *err = tge->GetEY();
  int npt = tge->GetN();
  for (int ipt=0; ipt<npt; ++ipt) {
    if ( val[ipt] > max ) {
      max = val[ipt];
      max_e = err[ipt];
    }
  }
  return max_e;
}

double get_max_val_error(TGraphErrors *tge) {
  double *val = tge->GetY();
  double *err = tge->GetEY();
  double max_e = 0.0;
  double max = 0.0;
  int npt = tge->GetN();
  for (int ipt=0; ipt<npt; ++ipt) {
    if ( val[ipt] > max ) {
      max = val[ipt];
      max_e = err[ipt];
    }
  }
  return max_e;
}

double get_max(const int& n, TGraphErrors *tge[]) {
  double max = 0.0;
  double max_e = 0.0;
  for (int igr=0; igr<n; ++igr) {
    double _max=0, _max_e=0;
    get_max_val_and_error(tge[igr],_max,_max_e);
    cout << "max = " << max << " tge[igr]->GetMaximum()= " << _max;
    if ( _max > max ) {
      max = _max;
      max_e = _max_e;
      //max_e = get_max_val_error(tge[igr]);
      //max_e = get_max_val_and_error(tge[igr]);
      cout << " max_val_error= " << max_e << endl;
    }
    cout << endl;
  }
  return max + max_e;
}

void plot_pres(int CASE=0, int SUB_CASE=0) {

  gStyle->SetOptStat(0);
  //gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(2.0,"y");
  //gStyle->SetTitleSize(0.1,"y");
  gStyle->SetLabelSize(0.05,"y");
  string fName;

  bool do_fit = true;

  const char *path = "output/report.v4/";

  int idx_dia = 4;
  int idx_had = 5;
  string s_tag = "";
  double manual_max_xyres = 6.0;
  if(CASE == 0 ) {
    fName = Form("%sms00_ds00_xy0.0mm_dth0_dph0_p1.30.root",path);
    s_tag = "0";
    do_fit = false;
  } else if (CASE == 1 && SUB_CASE == 0) {
    fName = Form("%sms00_ds00_xy0.5mm_dth10_dph50_p1.30.root",path);
    s_tag = "1";
  } else if (CASE == 1 && SUB_CASE == 1) { // 1a
    fName = Form("%sms00_ds00_xy0.5mm_dth0_dph0_p1.30.root",path);
    s_tag = "1a";
  } else if (CASE == 1 && SUB_CASE == 2) { // 1b
    fName = Form("%sms00_ds00_xy0.0mm_dth10_dph50_p1.30.root",path);
    s_tag = "1b";
  } else if (CASE == 2 && SUB_CASE == 0) {
    fName = Form("%sms00_ds10_xy0.5mm_dth10_dph50_p1.30.root",path);
    s_tag = "2";
    manual_max_xyres = 6;
  } else if (CASE == 2 && SUB_CASE == 1) { // 2a
    fName = Form("%sms10_ds10_xy0.5mm_dth10_dph50_p1.30.root",path);
    s_tag = "2a";
  } else if (CASE == 2 && SUB_CASE == 2) { // 2a
    fName = Form("%sms11_ds10_xy0.5mm_dth10_dph50_p1.30.root",path);
    s_tag = "2b";
  } else if (CASE == 3 && SUB_CASE == 0) {
    fName = Form("%sms00_ds11_xy0.5mm_dth10_dph50_p1.30.root",path);
    s_tag = "3";
    manual_max_xyres = 6;
  } else if (CASE == 3 && SUB_CASE == 1) { // 3a
    fName = Form("%sms10_ds11_xy0.5mm_dth10_dph50_p1.30.root",path);
    idx_dia = 6; idx_had = 7;
    s_tag = "3a";
  } else if (CASE == 3 && SUB_CASE == 2) { // 3b
    fName = Form("%sms11_ds11_xy0.5mm_dth10_dph50_p1.30.root",path);
    idx_dia = 6; idx_had = 8;
    s_tag = "3b";
  } else {
    return;
  }

  get_hist_limits(fName);
  double d_del[ndel] = {0.0}, d_del_e[ndel] = {0.0};
  const char  *cmeth[2] = {"TDR","MINUIT"};
  cout << "cmeth[0]= " << cmeth[0] << endl;
  cout << "cmeth[1]= " << cmeth[1] << endl;

  for (int idel=0; idel<ndel; ++idel) {
    d_del[idel] = idel - 5.0;
    d_del_e[idel] = 0.0;
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      hpres[imeth][idel] = (TH1F*) def_hist(Form("h_pres_dp%d_%s",idel,cmeth[imeth]), 1000, 1.2, 1.4, idel!=9?idel+1:42);
      //hpres[imeth][idel] = (TH1F*) def_hist(Form("h_pres_dp%d_%s",idel,cmeth[imeth]), 1000, 1.0, 1.6, idel!=9?idel+1:42);
      hxres_dia[imeth][idel] = (TH1F*) def_hist(Form("h_xres_dia_dp%d_%s",idel,cmeth[imeth]),100,-xmax_dia[imeth][idel],xmax_dia[imeth][idel], 1);
      hyres_dia[imeth][idel] = (TH1F*) def_hist(Form("h_yres_dia_dp%d_%s",idel,cmeth[imeth]),100,-ymax_dia[imeth][idel],ymax_dia[imeth][idel], 1);
      hxres_had[imeth][idel] = (TH1F*) def_hist(Form("h_xres_had_dp%d_%s",idel,cmeth[imeth]),100,-xmax_had[imeth][idel],xmax_had[imeth][idel], 1);
      hyres_had[imeth][idel] = (TH1F*) def_hist(Form("h_yres_had_dp%d_%s",idel,cmeth[imeth]),100,-ymax_had[imeth][idel],ymax_had[imeth][idel], 1);
    }
  }

  TFile *f = TFile::Open(fName.c_str());
  TTree *t = (TTree*) f->Get("t");

  int ndet;      t->SetBranchAddress("ndet",&ndet);
  int acc;       t->SetBranchAddress("acc",&acc);
  float p[10];   t->SetBranchAddress("p[ndet]",&p);
  float dp[10];  t->SetBranchAddress("dp[ndet]",&dp);
  float x[10];   t->SetBranchAddress("x[ndet]",&x);
  float y[10];   t->SetBranchAddress("y[ndet]",&y);
  int nrec;      t->SetBranchAddress("nrec",&nrec);
  float pr[10];   t->SetBranchAddress("pr[nrec]",&pr);
  float xr[10];   t->SetBranchAddress("xr[nrec]",&xr);
  float yr[10];   t->SetBranchAddress("yr[nrec]",&yr);
  float sr[10];   t->SetBranchAddress("rstat[nrec]",&sr);

  const int idx_pr[nmeth] = {0,3};
  const int idx_xyr_dia[nmeth] = {1,4};
  const int idx_xyr_had[nmeth] = {2,5};

  int nent = t->GetEntries();
  //nent = 10000;
  cout << "Nent=  " << nent << endl;
  for (int ient=0; ient<nent; ++ient) {

    t->GetEntry(ient);
    int idel = floor(dp[0]+0.5) + 5;
    if (ient%10000==0) {
      //if (idel==5)
      cout << "ient= " << ient << "/"<< nent << " p[0]= " << p[0] << " pr[0]= " << pr[0] << " dp[0] = " << dp[0] << " idel= " << idel << endl;
    }

    if (idel<0 || idel>ndel) cout << "WTF! idel = " << idel << endl;
    //cout << "p[0]= " << p[0]
    if (sr[0] != 0 ) {
      cout << "status of fit bad " << endl;
      continue;
    }

    if (acc!=1) continue;

    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      hpres[imeth][idel]->Fill(pr[idx_pr[imeth]]);
      hxres_dia[imeth][idel]->Fill(xr[idx_xyr_dia[imeth]]-x[idx_dia]);
      hyres_dia[imeth][idel]->Fill(yr[idx_xyr_dia[imeth]]-y[idx_dia]);
      hxres_had[imeth][idel]->Fill(xr[idx_xyr_had[imeth]]-x[idx_had]);
      hyres_had[imeth][idel]->Fill(yr[idx_xyr_had[imeth]]-y[idx_had]);
    }
  }

  TCanvas *tc_xyres[nmeth];
  TCanvas *tc_xres_dia[nmeth];
  TCanvas *tc_yres_dia[nmeth];
  TCanvas *tc_xres_had[nmeth];
  TCanvas *tc_yres_had[nmeth];
  double d_xres_dia[nmeth][ndel] = {{0.0}}, d_xres_dia_e[nmeth][ndel] = {{0.0}};
  double d_yres_dia[nmeth][ndel] = {{0.0}}, d_yres_dia_e[nmeth][ndel] = {{0.0}};
  double d_xres_had[nmeth][ndel] = {{0.0}}, d_xres_had_e[nmeth][ndel] = {{0.0}};
  double d_yres_had[nmeth][ndel] = {{0.0}}, d_yres_had_e[nmeth][ndel] = {{0.0}};

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.08);
  tt->SetTextColor(4);
  tt->SetNDC(true);

  double max_xyres = 0.0;
  for (int imeth=0; imeth<nmeth_max; ++imeth) {

    tc_xyres[imeth] = new TCanvas(Form("tc_xres_dia_%s",cmeth[imeth]),Form("tc_xres_dia_%s",cmeth[imeth]),1600,1000);
    tc_xyres[imeth]->Divide(8,6);
    //tc_xres_dia[imeth] = new TCanvas(Form("tc_xres_dia_%s",cmeth[imeth]),Form("tc_xres_dia_%s",cmeth[imeth]));
    //tc_xres_dia[imeth]->Divide(4,3);
    //tc_yres_dia[imeth] = new TCanvas(Form("tc_yres_dia_%s",cmeth[imeth]),Form("tc_yres_dia_%s",cmeth[imeth]));
    //tc_yres_dia[imeth]->Divide(4,3);
    //tc_xres_had[imeth] = new TCanvas(Form("tc_xres_had_%s",cmeth[imeth]),Form("tc_xres_had_%s",cmeth[imeth]));
    //tc_xres_had[imeth]->Divide(4,3);
    //tc_yres_had[imeth] = new TCanvas(Form("tc_yres_had_%s",cmeth[imeth]),Form("tc_yres_had_%s",cmeth[imeth]));
    //tc_yres_had[imeth]->Divide(4,3);

    for (int idel=0; idel<ndel; ++idel) {

      //fxres_dia[imeth][idel] = (TF1*) fit_xy_res(d_xres_dia[imeth][idel], d_xres_dia_e[imeth][idel], max_xyres, do_fit, idel, tc_xres_dia[imeth], hxres_dia[imeth][idel], cmeth[imeth], "xres", "dia", "X_{rec}-X", "Diam. Det.");
      //fyres_dia[imeth][idel] = (TF1*) fit_xy_res(d_yres_dia[imeth][idel], d_yres_dia_e[imeth][idel], max_xyres, do_fit, idel, tc_yres_dia[imeth], hyres_dia[imeth][idel], cmeth[imeth], "yres", "dia", "Y_{rec}-Y", "Diam. Det." );
      //fxres_had[imeth][idel] = (TF1*) fit_xy_res(d_xres_had[imeth][idel], d_xres_had_e[imeth][idel], max_xyres, do_fit, idel, tc_xres_had[imeth], hxres_had[imeth][idel], cmeth[imeth], "xres", "had", "X_{rec}-X", "HADES" );
      //fyres_had[imeth][idel] = (TF1*) fit_xy_res(d_yres_had[imeth][idel], d_yres_had_e[imeth][idel], max_xyres, do_fit, idel, tc_yres_had[imeth], hyres_had[imeth][idel], cmeth[imeth], "yres", "had", "Y_{rec}-Y", "HADES" );

      const int ii = ((idel/4)*8) + (idel%4);
      cout << "idel= " << idel << " ii= " << ii<< endl;

      fxres_dia[imeth][idel] = (TF1*)
	fit_xy_res(d_xres_dia[imeth][idel], d_xres_dia_e[imeth][idel], max_xyres, do_fit, idel, ii,
		   tc_xyres[imeth], hxres_dia[imeth][idel], cmeth[imeth], "xres", "dia", "X_{rec}-X", "Diam. Det.");

      fyres_dia[imeth][idel] = (TF1*)
	fit_xy_res(d_yres_dia[imeth][idel], d_yres_dia_e[imeth][idel], max_xyres, do_fit, idel, 4+ii,
		   tc_xyres[imeth], hyres_dia[imeth][idel], cmeth[imeth], "xres", "dia", "Y_{rec}-Y", "Diam. Det." );

      fxres_had[imeth][idel] = (TF1*)
	fit_xy_res(d_xres_had[imeth][idel], d_xres_had_e[imeth][idel], max_xyres, do_fit, idel, 24+ii, tc_xyres[imeth],
		   hxres_had[imeth][idel], cmeth[imeth], "xres", "had", "X_{rec}-X", "HADES" );

      fyres_had[imeth][idel] = (TF1*)
	fit_xy_res(d_yres_had[imeth][idel], d_yres_had_e[imeth][idel], max_xyres, do_fit, idel, 28+ii, tc_xyres[imeth],
		   hyres_had[imeth][idel], cmeth[imeth], "xres", "had", "Y_{rec}-Y", "HADES" );

    }

    //tc_xres_dia[imeth]->Print(Form("slides/v4/%s_case%s.eps",tc_xres_dia[imeth]->GetName(),s_tag.c_str()));
    //tc_yres_dia[imeth]->Print(Form("slides/v4/%s_case%s.eps",tc_yres_dia[imeth]->GetName(),s_tag.c_str()));
    //tc_xres_had[imeth]->Print(Form("slides/v4/%s_case%s.eps",tc_xres_had[imeth]->GetName(),s_tag.c_str()));
    //tc_yres_had[imeth]->Print(Form("slides/v4/%s_case%s.eps",tc_yres_had[imeth]->GetName(),s_tag.c_str()));

  }

  save_hist_limits(fName);

  TF1 *fpres[nmeth][ndel];
  double d_pres[nmeth][ndel] = {{0.0}}, d_pres_e[nmeth][ndel] = {{0.0}};
  if (do_fit) {
    for (int idel=0; idel<ndel; ++idel) {
      for (int imeth=0; imeth<nmeth_max; ++imeth) {
	//const char  *cmeth[imeth] = (imeth==0?"tdr":"min");
	double fit_min= hpres[imeth][idel]->GetMean()-1.5*hpres[imeth][idel]->GetRMS();
	double fit_max= hpres[imeth][idel]->GetMean()+1.5*hpres[imeth][idel]->GetRMS();
	fpres[imeth][idel] = new TF1(Form("fpres_%d_%s",idel,cmeth[imeth]), "gaus", fit_min, fit_max);
	fpres[imeth][idel]->SetLineColor(idel!=9?idel+1:42);
      }
    }
    TCanvas *tc_tmp = new TCanvas("tmp","tmp");
    tc_tmp->cd();
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      for (int idel=0; idel<ndel; ++idel) {
	hpres[imeth][idel]->Fit(fpres[imeth][idel],"R+");

	double fwhm = fpres[imeth][idel]->GetParameter(2);
	double mean = fpres[imeth][idel]->GetParameter(1);
	double fwhm_re = fpres[imeth][idel]->GetParError(2)/fwhm;
	double mean_re = fpres[imeth][idel]->GetParError(1)/mean;
	d_pres[imeth][idel] = 100*fwhm/mean;
	d_pres_e[imeth][idel] = d_pres[imeth][idel]*(TMath::Hypot(fwhm_re,mean_re));

	//double fwhm = fpres[imeth][idel]->GetParameter(2);
	//double fwhm_re = fpres[imeth][idel]->GetParError(2);
	//d_pres[imeth][idel] = 100*fwhm;
	//d_pres_e[imeth][idel] = 100*fwhm_re;

      }
    }
  }

  double max[nmeth] = {0.0};
  int imaxdel[nmeth] = {0};
  for (int imeth=0; imeth<nmeth_max; imeth++) {
    for (int idel=0; idel<ndel; ++idel) {
      if (hpres[imeth][idel]->GetMaximum() > max[imeth]) {
	max[imeth] = hpres[imeth][idel]->GetMaximum();
	imaxdel[imeth] = idel;
      }
    }
  }

  TCanvas *tc_pres[nmeth];
  for (int imeth=0; imeth<nmeth_max; imeth++) {
    //const char  *cmeth[imeth] = (imeth==0?"tdr":"min");
    tc_pres[imeth] = new TCanvas(Form("tc_pres_%s",cmeth[imeth]),Form("tc_pres_%s",cmeth[imeth]));
    tc_pres[imeth]->cd();
    hpres[imeth][imaxdel[imeth]]->Draw();
    hpres[imeth][imaxdel[imeth]]->SetTitle(Form("Reconstructed Momentum Distribution (%s); Reconstructed Momentum (p_{rec}[GeV/c]); dN/dp_{rec}",cmeth[imeth]));
    for (int idel=0; idel<ndel; ++idel) {
      if (idel==imaxdel[imeth]) continue;
      hpres[imeth][idel]->Draw("same");
    }
    tc_pres[imeth]->Print(Form("slides/v4/%s_case%s.eps",tc_pres[imeth]->GetName(), s_tag.c_str()));
  }

  if (do_fit) {

    TGraphErrors *tge_pres[nmeth];
    TGraphErrors *tge_xres_dia[nmeth];
    TGraphErrors *tge_yres_dia[nmeth];
    TGraphErrors *tge_xres_had[nmeth];
    TGraphErrors *tge_yres_had[nmeth];

    int meth_color[2] = {2, 4};
    //int meth_color[2] = {30, 46};
    //int meth_marker[2] = {20, 22};
    int det_marker[2] = {33, 20};
    double det_marker_size[2] = {1.5, 1};

    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      tge_pres[imeth]= new TGraphErrors(ndel, d_del, d_pres[imeth], d_del_e, d_pres_e[imeth]);
      style_tgraph(tge_pres[imeth], Form("tge_pres_%s",cmeth[imeth]),
		   Form("Momentum Resolution (%s); mom offset(%%); mom res (%%)",cmeth[imeth]), meth_color[imeth], 20, 1.);

      tge_xres_dia[imeth]= new TGraphErrors(ndel, d_del, d_xres_dia[imeth], d_del_e, d_xres_dia_e[imeth]);
      style_tgraph(tge_xres_dia[imeth], Form("tge_xres_dia_%s",cmeth[imeth]),
		   "X Pos Resolution (mm); mom offset(%%); Position resolution (mm)", meth_color[imeth], det_marker[0], det_marker_size[0]);

      tge_yres_dia[imeth] = new TGraphErrors(ndel, d_del, d_yres_dia[imeth], d_del_e, d_yres_dia_e[imeth]);
      style_tgraph(tge_yres_dia[imeth], Form("tge_yres_dia_%s",cmeth[imeth]),
		   "Y Pos Resolution (mm); mom offset(%%); Position resolution (mm)", meth_color[imeth], det_marker[0], det_marker_size[0]);

      tge_xres_had[imeth] = new TGraphErrors(ndel, d_del, d_xres_had[imeth], d_del_e, d_xres_had_e[imeth]);
      style_tgraph(tge_xres_had[imeth], Form("tge_xres_had_%s",cmeth[imeth]),
		   "X Pos Resolution (mm); mom offset(%%); Position resolution (mm)", meth_color[imeth], det_marker[1], det_marker_size[1]);

      tge_yres_had[imeth] = new TGraphErrors(ndel, d_del, d_yres_had[imeth], d_del_e, d_yres_had_e[imeth]);
      style_tgraph(tge_yres_had[imeth], Form("tge_yres_had_%s",cmeth[imeth]),
		   "Y Pos Resolution (mm); mom offset(%%); Position resolution (mm)", meth_color[imeth], det_marker[1], det_marker_size[1]);
    }

    TMultiGraph *tmg_xres = new TMultiGraph();
    tmg_xres->SetTitle(Form("X - Position Resolution (mm); mom offset (%%); Position Reoslution (mm)"));
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      //tmg_xres->Add(tge_xres_dia[imeth],"p");
      tmg_xres->Add(tge_xres_had[imeth],"p");
    }

    TLegend *tl_xres = new TLegend(0.3,0.5,1.4,0.9);
    tl_xres->SetFillStyle(0);
    tl_xres->SetBorderSize(0);
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      //tl_xres->AddEntry(tge_xres_dia[imeth],Form("Diamond, %s",cmeth[imeth]),"pl");
      tl_xres->AddEntry(tge_xres_had[imeth],Form("HADES, %s",cmeth[imeth]),"pl");
    }

    TMultiGraph *tmg_yres = new TMultiGraph();
    tmg_yres->SetTitle(Form("Y - Position Resolution (mm); mom offset (%%); Position Reoslution (mm)"));
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      //tmg_yres->Add(tge_yres_dia[imeth],"p");
      tmg_yres->Add(tge_yres_had[imeth],"p");
    }

    TLegend *tl_yres = new TLegend(0.0,0.55,0.99,0.9);
    tl_yres->SetFillStyle(0);
    tl_yres->SetBorderSize(0);
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      //tl_yres->AddEntry(tge_yres_dia[0],Form("Diamond, %s",cmeth[0]),"pl");
      tl_yres->AddEntry(tge_yres_had[0],Form("HADES, %s",cmeth[0]),"pl");
    }

    TMultiGraph *tmg_pres = new TMultiGraph();
    tmg_pres->SetTitle("Momentum Resolution (%); mom offset (%); Momentum Reoslution (%)");
    for (int imeth=0; imeth<nmeth_max; ++imeth) tmg_pres->Add(tge_pres[imeth],"p");

    TLegend *tl_pres = new TLegend(0.5,0.65,0.9,0.85);
    tl_pres->SetFillStyle(0);
    tl_pres->SetBorderSize(0);
    for (int imeth=0; imeth<nmeth_max; ++imeth) tl_pres->AddEntry(tge_pres[imeth],cmeth[imeth],"pl");

    //TCanvas *tc_xyres = new TCanvas("tc_xyres","tc_xyres");
    //tmg_xyres->Draw("a");
    //tl_xyres->Draw();

//    double max_xyres = 0.0;
//    for (int imeth=0; imeth<nmeth_max; ++imeth) {
//      max_xyres = TMath::Max(max_xyres, tge_xres_dia[imeth]->GetMaximum() );
//      max_xyres = TMath::Max(max_xyres, tge_yres_dia[imeth]->GetMaximum() );
//      max_xyres = TMath::Max(max_xyres, tge_xres_had[imeth]->GetMaximum() );
//      max_xyres = TMath::Max(max_xyres, tge_yres_had[imeth]->GetMaximum() );
//    }

    cout << "Calling for x"<< endl;
    double max_xres = get_max(nmeth_max, tge_xres_had);
    cout << "max_xres = " << max_xres << endl;
    cout << "Calling for y"<< endl;
    double max_yres = get_max(nmeth_max, tge_yres_had);
    cout << "max_yres = " << max_yres << endl;
    double max_pres = get_max(nmeth_max, tge_pres);

    TCanvas *tc_res = new TCanvas("tc_res","tc_res",1200,600);
    tc_res->Divide(3,1);

    tc_res->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    tmg_pres->Draw("a");
    tmg_pres->SetMinimum(0.0);
    //tmg_pres->SetMaximum(0.6);
    tmg_pres->SetMaximum(max_pres*1.1);
    tl_pres->Draw();

    tc_res->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    tmg_xres->Draw("a");
    tmg_xres->SetMinimum(0.0);
    tmg_xres->SetMaximum(max_xres*1.1);
    //tmg_xres->SetMaximum(8);
    //tl_xres->Draw();

    tc_res->cd(3);
    gPad->SetGridx();
    gPad->SetGridy();
    tmg_yres->Draw("a");
    tmg_yres->SetMinimum(0.0);
    tmg_yres->SetMaximum(max_yres*1.1);
    //tmg_yres->SetMaximum(8);
    //tl_yres->Draw();

    tc_res->ls();
    double xl=0, xu=0, yl=0, yu=0;
    TPad *pad1 = (TPad*) tc_res->FindObject("tc_res_1");
    pad1->GetPadPar(xl,yl,xu,yu);
    cout << "pad1: xl= " << xl << "xu= " << xu << "yl= " << yl  << "yu= " << yu  << endl;
    xu= 0.44;
    pad1->SetPad(xl,yl,xu,yu);
    TPad *pad2 = (TPad*) tc_res->FindObject("tc_res_2");
    pad2->GetPadPar(xl,yl,xu,yu);
    cout << "pad2: xl= " << xl << "xu= " << xu << "yl= " << yl  << "yu= " << yu  << endl;
    xl= 0.45; xu=0.735;
    pad2->SetPad(xl,yl,xu,yu);
    TPad *pad3 = (TPad*) tc_res->FindObject("tc_res_3");
    pad3->GetPadPar(xl,yl,xu,yu);
    cout << "pad3: xl= " << xl << "xu= " << xu << "yl= " << yl  << "yu= " << yu  << endl;
    xl= 0.735; xu=0.99;
    pad3->SetPad(xl,yl,xu,yu);

    double eps = 1e-9;
    pad3->SetLeftMargin(0.1);
    pad3->SetRightMargin(0.1);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(eps);
    pad1->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.15);
    tc_res->Print(Form("slides/v4/%s_case%s.eps",tc_res->GetName(),s_tag.c_str()) );

    TFile *rootfout = TFile::Open(Form("%s_graphs.root", fName.c_str()),"RECREATE");
    rootfout->cd();
    for (int imeth=0; imeth<nmeth_max; ++imeth) {
      tge_pres[imeth]->Write();
      tge_xres_dia[imeth]->Write();
      tge_yres_dia[imeth]->Write();
      tge_xres_had[imeth]->Write();
      tge_yres_had[imeth]->Write();
    }
    rootfout->Close();
  }

}

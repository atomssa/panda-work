#include "draw_pad1.C"
#include "draw_pad2.C"
#include "draw_pad3.C"
#include "draw_pad4.C"
#include "draw_pad5.C"
#include "draw_pad6.C"
#include "draw_pad7.C"
#include "pad_setup.C"

void draw_pad1(TCanvas*, int, const char*, const char*, double, double, double);
void draw_pad2(int, TCanvas*, int, const char*, const char*, double, double, double);
void draw_pad3(TCanvas*, int, const char*, const char*, double, double, double);
void draw_pad4(TCanvas*, int, const char*, const char*, double, double, double);
TPad* draw_pad5(int, TPad*, const char*, const char*, double, double, double, TF1*, TF1*, double *);
void draw_pad6(TCanvas*, int, const char*, const char*, double, double, double);
void draw_pad7(TCanvas*, int, const char*, const char*, double, double, double);
TPad* padSetup3x2(TCanvas*, int);
TPad* padSetup3x3(TCanvas*, int);

const int nconfig = 27;
Double_t config[nconfig][3] = {
  /*0*/{0.15, 30.0, 45.0}, /*1*/{0.25, 30.0, 45.0}, /*2*/{0.5, 30.0, 45.0}, /*3*/{0.5, 45.0, 60.0}, /*4*/{0.5, 60.0, 75.0},
  /*5*/{0.5, 75.0, 90.0}, /*6*/{0.5, 90.0, 105.0}, /*7*/{0.5, 105.0, 120.0}, /*8*/{1.0, 30.0, 45.0}, /*9*/{1.0, 45.0, 60.0},
  /*10*/{1.0, 60.0, 75.0}, /*11*/{1.0, 75.0, 90.0}, /*12*/{1.5, 30.0, 45.0}, /*13*/{1.5, 45.0, 60.0}, /*14*/{2.0, 30.0, 45.0},
  /*15*/{0.15, 10.0, 20.0}, /*16*/{0.25, 10.0, 20.0}, /*17*/{0.5, 10.0, 20.0}, /*18*/{1.0, 10.0, 20.0}, /*19*/{1.5, 10.0, 20.0},
  /*20*/{2.0, 10.0, 20.0}, /*21*/{2.5, 10.0, 20.0}, /*22*/{0.2, 30.0, 45.0}, /*23*/{0.2, 90.0, 105.0}, /*24*/{0.2, 10.0, 20.0},
  /* manually added below*/
  /*25*/{1.5, 10.0, 15.0}, /*26*/{1.5, 15.0, 20.0}
};

const int nbs=8;
const double epbs_binSize[nbs]={2.25,1,1.5,1.8,2.0,2.5,3.0,4.0};

// low energy- top: barrel, bot: endcap
const bool draw_a = false;
const int ncfg_a = 6;
int icfg_a[ncfg_a] = {0, 22, 1, 15, 24, 16};

// polar angle dependence, 500MeV from 10-105
const bool draw_b = false;
const int ncfg_b = 6;
int icfg_b[ncfg_b] = {17, 2, 3, 4, 5, 6};

// Energy dependence Endcap
const bool draw_c = false;
const int ncfg_c = 6;
int icfg_c[ncfg_c] = {24, 16, 17, 18, 19, 20};
//const int ncfg_c = 7;
//int icfg_c[ncfg_c] = {24, 16, 17, 18, 19, 20, 15};

// Energy dependence Barrel
const bool draw_d = false;
const int ncfg_d = 6;
int icfg_d[ncfg_d] = {22, 1, 2, 8, 12, 14};
//const int ncfg_d = 7;
//int icfg_d[ncfg_d] = {22, 1, 2, 8, 12, 14, 0};

// Electron muon comparison
const bool draw_e = false;
const int ncfg_e = 2;
int icfg_e[ncfg_e] = {8, 14};

// binning comparison
const bool draw_f = false;
const int ncfg_f = 2;
int icfg_f[ncfg_f] = {12, 19};

// weighting comparison
const bool draw_g = true;
const int ncfg_g = 2;
int icfg_g[ncfg_g] = {2, 17};

// 10-15 vs 15-20 comp in fwd. 1.5GeV/c
const bool draw_h = false;
const int ncfg_h = 2;
int icfg_h[ncfg_h] = {25, 26};

void set_style(TMultiGraph* tmg) {
  TH1F* h = tmg->GetHistogram();
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  h->GetXaxis()->SetNdivisions(610);
}

void draw_graphs(int icfg, TCanvas *tc, const char *name, const char *title,
		 double min, double max, int start, int num, double* x, double* xe,
		 double* rec, double* rece, double* cor, double* core) {

  TGraphErrors *tge_rec = new TGraphErrors(num, &x[start], &rec[start], &xe[start], &rece[start]);
  tge_rec->SetMarkerStyle(20);
  tge_rec->SetMarkerSize(1.5);
  tge_rec->SetMarkerColor(1);
  tge_rec->SetLineColor(1);
  tge_rec->SetLineWidth(2);
  TGraphErrors *tge_cor = new TGraphErrors(num, &x[start], &cor[start], &xe[start], &core[start]);
  tge_cor->SetMarkerStyle(20);
  tge_cor->SetMarkerSize(1.5);
  tge_cor->SetMarkerColor(2);
  tge_cor->SetLineColor(1);
  tge_cor->SetLineWidth(2);

  TMultiGraph *tmg = new TMultiGraph(name,title);
  tmg->Add(tge_rec,"p");
  tmg->Add(tge_cor,"p");

  tc->cd();
  gPad->SetGridy();

  tmg->SetMinimum(min);
  tmg->SetMaximum(max);
  tmg->Draw("a");
  set_style(tmg);
  tc->Update();

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);
  tt->SetTextSize(0.06);

  double pt=config[icfg][0], thmin=config[icfg][1], thmax=config[icfg][2];
  tt->SetTextColor(1);
  tt->DrawLatex(0.58,0.24,"Uncorrected");
  tt->SetTextColor(2);
  tt->DrawLatex(0.61,0.17,"Corrected");
  tt->SetTextColor(4);
  if (draw_b) {
    tt->DrawLatex(0.15,0.17,Form("p_{T}= %4.2g GeV/c",pt));
  } else if (draw_c || draw_d ) {
    tt->DrawLatex(0.15,0.24,Form("%3.0f#circ < #theta < %3.0f#circ",thmin,thmax));
  }
  if (!draw_b) {
    bool fwd = thmin < 15.0 && 15.0<thmax;
    tt->DrawLatex(fwd?0.15:0.2,0.17,fwd?"(Fwd Endcap)":"(Barrel)");
  }
}

void draw_graphs(TCanvas *tc, const char* name, const char *title,
		 double min, double max, int start1, int num1, int start2, int num2,
		 double* x, double* xe, double* rec, double* rece, double* cor, double* core) {

  TGraphErrors *tge_rec1 = new TGraphErrors(num1, &x[start1], &rec[start1], &xe[start1], &rece[start1]);
  tge_rec1->SetMarkerStyle(20);
  tge_rec1->SetMarkerColor(1);
  TGraphErrors *tge_cor1 = new TGraphErrors(num1, &x[start1], &cor[start1], &xe[start1], &core[start1]);
  tge_cor1->SetMarkerStyle(20);
  tge_cor1->SetMarkerColor(2);

  TGraphErrors *tge_rec2 = new TGraphErrors(num2, &x[start2], &rec[start2], &xe[start2], &rece[start2]);
  tge_rec2->SetMarkerStyle(21);
  tge_rec2->SetMarkerColor(1);
  TGraphErrors *tge_cor2 = new TGraphErrors(num2, &x[start2], &cor[start2], &xe[start2], &core[start2]);
  tge_cor2->SetMarkerStyle(21);
  tge_cor2->SetMarkerColor(2);

  TMultiGraph *tmg = new TMultiGraph(name,title);
  tmg->Add(tge_rec1,"p");
  tmg->Add(tge_cor1,"p");
  tmg->Add(tge_rec2,"p");
  tmg->Add(tge_cor2,"p");

  tc->cd();
  gPad->SetGridy();

  tmg->SetMinimum(min);
  tmg->SetMaximum(max);
  tmg->Draw("a");
  set_style(tmg);
  tc->Update();

}

void normalize_canv(TCanvas *canv, int n, const char* in, const char* out) {
  double max = 0;
  for (int ic=0; ic< n; ++ic) {
    TPad* pad = (TPad*) canv->GetPrimitive(Form("pads%d",ic));
    TH1F* htmp = (TH1F*) pad->GetPrimitive(in);
    if (htmp->GetMaximum()>max) max = htmp->GetMaximum();
  }
  for (int ic=0; ic< n; ++ic) {
    TPad* pad = (TPad*) canv->GetPrimitive(Form("pads%d",ic));
    TH1F* htmp = (TH1F*) pad->GetPrimitive(out);
    htmp->SetMaximum(max*1.1);
  }
}

void res(int ibs=0) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  TF1 *rec_gaus = new TF1("rec_gaus","gaus",-0.02,0.02);
  TF1 *cor_gaus = new TF1("cor_gaus","gaus",-0.02,0.02);
  double x[nconfig];
  double xe[nconfig];
  // gc= gaussian center, f2s= fraction in 2sigma
  double rec_gc[nconfig], rec_gce[nconfig];
  double rec_f2s[nconfig], rec_f2se[nconfig];
  double cor_gc[nconfig], cor_gce[nconfig];
  double cor_f2s[nconfig], cor_f2se[nconfig];

  const char *a = "../grid.out/";
  //const char *b = "_oct14_binsong_configs/bremcorr.all.ibs.cfg.";
  const char *b = "_oct14_binsong_configs/all.ibs.xr-0.1_0.3/bremcorr.all.ibs.cfg.";

  // low energy- top: barrel, bot: endcap
  if (draw_a) {
    TCanvas *tc_a = new TCanvas("tc_a","tc_a",1000,1500);
    //TCanvas *tc_a_gc = new TCanvas("tc_a_gc","tc_a_gc",720,10,700,700);
    //TCanvas *tc_a_f2s = new TCanvas("tc_a_f2s","tc_a_f2s",1440,10,700,700);

    for (int ic = 0; ic < ncfg_a; ++ic) {
      double pt=config[icfg_a[ic]][0], thmin=config[icfg_a[ic]][1], thmax=config[icfg_a[ic]][2];
      double tmp[4];
      draw_pad5(ic, padSetup2x3vert(tc_a, ic), a, Form("%s%d_hists.root", b, icfg_b[ic]), pt, thmin, thmax, rec_gaus, cor_gaus, tmp);
      x[ic] = pt; xe[ic] = 0;
      rec_gc[ic] = rec_gaus->GetParameter(1); rec_gce[ic] = rec_gaus->GetParError(1);
      rec_f2s[ic] = tmp[0]; rec_f2se[ic] = tmp[1];
      cor_gc[ic] = cor_gaus->GetParameter(1); cor_gce[ic] = cor_gaus->GetParError(1);
      cor_f2s[ic] = tmp[2]; cor_f2se[ic] = tmp[3];
    }

    //draw_graphs(icfg_a[0], tc_a_gc, "tmg_gc_a", ";p_{T}[GeV/c];Peak position", -0.02, 0.02, 0, 3, 3, 3, x, xe, rec_gc, rec_gce, cor_gc, cor_gce);
    //draw_graphs(icfg_a[0], tc_a_f2s, "tmg_f2s_a", ";p_{T}[GeV/0];N_{2#sigma}/N_{tot}", 0.0, 1.0, 0, 3, 3, 3, x, xe, rec_f2s, rec_f2se, cor_f2s, cor_f2se);

    //normalize_canv(tc_a, ncfg_a, "h_rec_mu", "h_rec_mu");

    tc_a->Print("note/perf_low_pt_style6.pdf");

  }

  if (draw_b) {

    TCanvas *tc_b = new TCanvas("tc_b","tc_b",1000,1500);
    //TCanvas *tc_b_gc = new TCanvas("tc_b_gc","tc_b_gc",720,10,700,700);
    //TCanvas *tc_b_f2s = new TCanvas("tc_b_f2s","tc_b_f2s",1440,10,700,700);

    for (int ic = 0; ic < ncfg_b; ++ic) {
      //TString mufn = Form("%s/mum%s%d_hists.root",a, b, icfg_b[ic]),
      //elfn = Form("%s/esim%s%d_hists.root",a, b, icfg_b[ic]);
      double pt=config[icfg_b[ic]][0], thmin=config[icfg_b[ic]][1], thmax=config[icfg_b[ic]][2];
      double tmp[4];
      draw_pad5(ic, padSetup2x3vert(tc_b, ic), a, Form("%s%d_hists.root", b, icfg_b[ic]), pt, thmin, thmax, rec_gaus, cor_gaus, tmp);
      x[ic] = (thmin+thmax)/2; xe[ic] = 0;
      rec_gc[ic] = rec_gaus->GetParameter(1); rec_gce[ic] = rec_gaus->GetParError(1);
      rec_f2s[ic] = tmp[0]; rec_f2se[ic] = tmp[1];
      cor_gc[ic] = cor_gaus->GetParameter(1); cor_gce[ic] = cor_gaus->GetParError(1);
      cor_f2s[ic] = tmp[2]; cor_f2se[ic] = tmp[3];
    }

    //draw_graphs(icfg_b[0], tc_b_gc, "tmg_gc_b", ";#theta(#circ);Peak position", -0.004, 0.006, 0, 6, x, xe, rec_gc, rec_gce, cor_gc, cor_gce);
    //draw_graphs(icfg_b[0], tc_b_f2s, "tmg_f2s_b", ";#theta(#circ);N_{2#sigma}/N_{tot}", 0.0, 1.0, 0, 6, x, xe, rec_f2s, rec_f2se, cor_f2s, cor_f2se);

    //normalize_canv(tc_b, ncfg_b, "h_rec_mu", "h_rec_mu");

    tc_b->Print("note/perf_angular_dep_style6.pdf");
    //tc_b_gc->Print("note/perf_angular_dep_graph_gc.pdf");
    //tc_b_f2s->Print("note/perf_angular_dep_graph_f2s.pdf");
  }

  if (draw_c) {
    TCanvas *tc_c = new TCanvas("tc_c","tc_c",1000,1500);
    //TCanvas *tc_c_gc = new TCanvas("tc_c_gc","tc_c_gc",700,10,700,700);
    //TCanvas *tc_c_f2s = new TCanvas("tc_c_f2s","tc_c_f2s",1400,10,700,700);

    for (int ic = 0; ic < ncfg_c; ++ic) {
      double pt=config[icfg_c[ic]][0], thmin=config[icfg_c[ic]][1], thmax=config[icfg_c[ic]][2];
      double tmp[4];
      draw_pad5(ic, padSetup2x3vert(tc_c, ic), a, Form("%s%d_hists.root", b, icfg_c[ic]), pt, thmin, thmax, rec_gaus, cor_gaus, tmp);
      //draw_pad5(ic, padSetup3x3(tc_c, ic), a, Form("%s%d_hists.root", b, icfg_c[ic]), pt, thmin, thmax, rec_gaus, cor_gaus, tmp);
      x[ic] = pt; xe[ic] = 0;
      rec_gc[ic] = rec_gaus->GetParameter(1); rec_gce[ic] = rec_gaus->GetParError(1);
      rec_f2s[ic] = tmp[0]; rec_f2se[ic] = tmp[1];
      cor_gc[ic] = cor_gaus->GetParameter(1); cor_gce[ic] = cor_gaus->GetParError(1);
      cor_f2s[ic] = tmp[2]; cor_f2se[ic] = tmp[3];
    }

    //draw_graphs(icfg_c[0], tc_c_gc, "tmg_gc_c", ";p_{T}[GeV/c];Peak position", -0.005, 0.01, 0, 6, x, xe, rec_gc, rec_gce, cor_gc, cor_gce);
    //draw_graphs(icfg_c[0], tc_c_f2s, "tmg_f2s_c", ";p_{T}[GeV/c];N_{2#sigma}/N_{tot}", 0.0, 1.0, 0, 6, x, xe, rec_f2s, rec_f2se, cor_f2s, cor_f2se);

    //normalize_canv(tc_c, ncfg_c, "h_rec_mu", "h_rec_mu");

    tc_c->Print("note/perf_ene_dep_fwd_style6.pdf");
    //tc_c_gc->Print("note/perf_ene_dep_fwd_grph_gc.pdf");
    //tc_c_f2s->Print("note/perf_ene_dep_fwd_grph_gf2s.pdf");
  }

  // Energy dependence Barrel
  if (draw_d) {
    cout << "drawing d" << endl;
    TCanvas *tc_d = new TCanvas("tc_d","tc_d",1000,1500);
    //TCanvas *tc_d_gc = new TCanvas("tc_d_gc","tc_d_gc",700,10,700,700);
    //TCanvas *tc_d_f2s = new TCanvas("tc_d_f2s","tc_d_f2s",1400,10,700,700);

    for (int ic = 0; ic < ncfg_d; ++ic) {
      double pt=config[icfg_d[ic]][0], thmin=config[icfg_d[ic]][1], thmax=config[icfg_d[ic]][2];
      double tmp[4];
      draw_pad5(ic, padSetup2x3vert(tc_d, ic), a, Form("%s%d_hists.root", b, icfg_d[ic]), pt, thmin, thmax, rec_gaus, cor_gaus, tmp);
      //draw_pad5(ic, padSetup3x3(tc_d, ic), a, Form("%s%d_hists.root", b, icfg_d[ic]), pt, thmin, thmax, rec_gaus, cor_gaus, tmp);
      x[ic] = pt; xe[ic] = 0;
      rec_gc[ic] = rec_gaus->GetParameter(1); rec_gce[ic] = rec_gaus->GetParError(1);
      rec_f2s[ic] = tmp[0]; rec_f2se[ic] = tmp[1];
      cor_gc[ic] = cor_gaus->GetParameter(1); cor_gce[ic] = cor_gaus->GetParError(1);
      cor_f2s[ic] = tmp[2]; cor_f2se[ic] = tmp[3];
    }

    //draw_graphs(icfg_d[0], tc_d_gc, "tmg_gc_d", ";p_{T}[GeV/c];Peak position", -0.005, 0.008, 0, 6, x, xe, rec_gc, rec_gce, cor_gc, cor_gce);
    //draw_graphs(icfg_d[0], tc_d_f2s, "tmg_f2s_d", ";p_{T}[GeV/c];N_{2#sigma}/N_{tot}", 0.0, 1.0, 0, 6, x, xe, rec_f2s, rec_f2se, cor_f2s, cor_f2se);

    //normalize_canv(tc_d, ncfg_d, "h_rec_mu", "h_rec_mu");

    tc_d->Print("note/perf_ene_dep_brl_style6.pdf");
    //tc_d_gc->Print("note/perf_ene_dep_brl_grph_gc.pdf");
    //tc_d_f2s->Print("note/perf_ene_dep_brl_grph_gf2s.pdf");
  }

  // e-mu comp
  if (draw_e) {
    cout << "drawing e" << endl;
    TCanvas *tc_e = new TCanvas("tc_e","tc_e",1000,500);
    tc_e->Divide(2,1);
    for (int ic = 0; ic < ncfg_e; ++ic) {
      draw_pad4(tc_e, ic,
		Form("%s/mum%s%d_hists.root",a, b, icfg_e[ic]),
		Form("%s/esim%s%d_hists.root",a, b, icfg_e[ic]),
		config[icfg_e[ic]][0], config[icfg_e[ic]][1], config[icfg_e[ic]][2]);
    }
  }

  // binning comparison
  if (draw_f) {
    cout << "drawing f" << endl;
    TCanvas *tc_f = new TCanvas("tc_f","tc_f",1000,500);
    tc_f->Divide(2,1);
    for (int ic = 0; ic < ncfg_f; ++ic) {
      b = "_oct14_binsong_configs/all.ibs.xr-0.1_0.1/bremcorr.all.ibs.cfg.";
      draw_pad3(tc_f, ic,
		Form("%s/mum%s%d_hists.root",a, b, icfg_f[ic]),
		Form("%s/esim%s%d_hists.root",a, b, icfg_f[ic]),
		config[icfg_f[ic]][0], config[icfg_f[ic]][1], config[icfg_f[ic]][2]);
    }
  }

  // weighting comparison
  if (draw_g) {
    cout << "drawing g" << endl;
    TCanvas *tc_g = new TCanvas("tc_g","tc_g",1000,500);
    tc_g->Divide(2,1);
    for (int ic = 0; ic < ncfg_g; ++ic) {
      b = "_oct14_binsong_configs/all.ibs.xr-0.1_0.1/bremcorr.all.ibs.cfg.";
      draw_pad6(tc_g, ic,
		Form("%s/mum%s%d_hists.root",a, b, icfg_g[ic]),
		Form("%s/esim%s%d_hists.root",a, b, icfg_g[ic]),
		config[icfg_g[ic]][0], config[icfg_g[ic]][1], config[icfg_g[ic]][2]);
    }
    tc_g->Print("note/perf_weighting.pdf");
  }

  // Fwd, 10-15 vs 15-20 comparison
  if (draw_h) {
    cout << "drawing h" << endl;
    TCanvas *tc_h = new TCanvas("tc_h","tc_h",1000,500);
    tc_h->Divide(2,1);
    for (int ic = 0; ic < ncfg_h; ++ic) {
      b = "_oct14_binsong_configs/all.ibs.xr-0.1_0.3/bremcorr.all.ibs.cfg.";
      draw_pad7(tc_h, ic,
		Form("%s/mum%s%d_hists.root",a, b, icfg_h[ic]),
		Form("%s/esim%s%d_hists.root",a, b, icfg_h[ic]),
		config[icfg_h[ic]][0], config[icfg_h[ic]][1], config[icfg_h[ic]][2]);
    }
  }

}

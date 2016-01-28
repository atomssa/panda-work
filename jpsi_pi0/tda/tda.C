#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TGAxis.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TBox.h"
#include "TLegend.h"
#include "TText.h"

// constants for prob. calculation
//_Mp = 0.938;
//_Mpi0 = 0.0; // This is required to reproduce x-sect dists. on PLB paper
//_Mjpsi2 = 9.59;
//pi = 3.14159;
//fpi = 0.093;// GeV
//gpiNN = 13;
//M0 = 8776.88; //GeV
//fphi = 0.413; //GeV
//fN = 0.005;//GeV
//alpha_s = 0.252062;
//M = 3.;//GeV
//C = pow(4.*pi*alpha_s,3)*(fN*fN*fphi/fpi)*(10./81.);

static const double pi = TMath::Pi();
static const double Mp = 0.938;
static const double Mp_sq = Mp*Mp;
static const double Mpi0 = 0.0; //0.135;
static const double Mpi0_sq = Mpi0*Mpi0;
static const double Mjpsi = 3.097;
static const double Mjpsi2 = Mjpsi*Mjpsi;
static const double msqTot = 2*Mp_sq + Mjpsi2 + Mpi0_sq;
static const double fpi = 0.093; // # GeV;
static const double fpsi = 0.413; //  #GeV;
static const double fN = 0.005; // #GeV;
static const double gpiNN = 13.0;
static const double alpha_s = 0.252062;
static const double Mbar = 3.0; // #GeV;
static const double C = pow(4.0*pi*alpha_s,3)*fN*fN*fpsi*10.0/fpi/81.0;
static const double M0 = 8776.88; //GeV

//static const double jpsi_pbarp_width = 1.2e-6;
//static const double M0 = sqrt(jpsi_pbarp_width*243.0*pi*pow(Mbar,5)/pow(pi*alpha_s,6)/1280.0/pow(fpsi,2)/pow(fN,4));
//print "M0 = %4.2f" % M0

//const double pbarn = 0.38941e9;
static const double pbarn = 0.38941*pow(10.0,9);

static const double Lumi = 2e3;
static const double Br = 5.94e-2;

//const double pbarn = 4e6;
double _ssss = 10.0;
const double _s[3] = {12.25, 16.87, 24.35};
const double _tmin[3] = {-0.45, -2.76, -6.5};
const double _tmax[3] = {0.6, 0.46, 0.32};
//const double _val_min[3] = {-0.092, -1.3, -2.85};
const double _val_min[3] = {-0.092, -1.0, -1.0};
const double _val_max[3] = {0.59, 0.43, 0.3};

//funciton definitions
double _lambda_sq(double x, double y, double z){
  return (x*x)+(y*y)+(z*z)-(2*x*y)-(2*x*z)-(2*y*z);
}

double _xi_a(double s){
  return Mbar*Mbar/(2.0*s - Mbar*Mbar);
}

double _xi(double Del_sq, double s){
  return (Mjpsi*Mjpsi - Del_sq - Mp_sq)/(2.0*s - Mjpsi*Mjpsi - Del_sq - 3.0*Mp_sq);
}

double _DelT_sq(double Del_sq, double s){
  return ((1-_xi(Del_sq,s))/(1+_xi(Del_sq,s))) * (Del_sq - 2.0*_xi(Del_sq,s)*( (Mp_sq/(1+_xi(Del_sq,s))) - (Mpi0_sq/(1-_xi(Del_sq,s))) ) );
}

double _DelT_sq_a(double Del_sq, double s){
  return ((1-_xi_a(s))/(1+_xi_a(s))) * (Del_sq - 2.0*_xi_a(s)*( (Mp_sq/(1+_xi_a(s))) - (Mpi0_sq/(1-_xi_a(s))) ) );
}

double _Del_sq_max(double s){
  return 2.0*_xi_a(s)*(Mp_sq*(_xi_a(s)-1) + Mpi0_sq*(_xi_a(s)+1))/(pow(_xi_a(s),2)-1);
}

double _dsig_dDel_sq(double Del_sq, double s){
  const double DelT_sq = _DelT_sq_a(Del_sq,s);
  const double F1 = 1/16.0/pi/_lambda_sq(s,Mp_sq,Mp_sq);
  const double F2 = C*C*2.0*(1+_xi(Del_sq,s))/4.0/_xi(Del_sq,s)/pow(Mbar, 8);
  const double I = fpi*gpiNN*Mp*(1-_xi(Del_sq,s))*M0/(Del_sq - Mp_sq)/(1+_xi(Del_sq,s));
  const double Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq);
  return Lumi*Br * pbarn * F1 * F2 * ( pow(I,2) - DelT_sq*pow(Iprim,2)/Mp_sq);
}

double _dsig_dDel_sq_a(double Del_sq, double s){
  const double DelT_sq = _DelT_sq_a(Del_sq,s);
  const double F1 = 1/16.0/pi/_lambda_sq(s,Mp_sq,Mp_sq);
  const double F2 = C*C*2.0*(1+_xi_a(s))/4.0/_xi_a(s)/pow(Mbar, 8);
  const double I = fpi*gpiNN*Mp*(1-_xi_a(s))*M0/(Del_sq - Mp_sq)/(1+_xi_a(s));
  const double Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq);
  //return Lumi*Br*pbarn*F1 * F2 * ( pow(I,2) - DelT_sq*pow(Iprim,2)/Mp_sq);
  return pbarn*F1 * F2 * ( pow(I,2) - DelT_sq*pow(Iprim,2)/Mp_sq);
}

//double anatda_func(double mand){
//  double _deltaT2 = ((1-_xi)/(1+_xi)) * (mand - 2.*_xi * ( (pow(Mp,2)/(1+_xi)) - (pow(Mpi0,2)/(1-_xi)) ));
//  double I = M0*fpi*gpiNN*Mp*(1-_xi)/(mand-Mp*Mp)/(1+_xi);
//  double Iprim = M0*fpi*gpiNN*Mp/(mand-Mp*Mp);
//  double MT2 = 0.25*C*C*(2*(_xi+1)/_xi/pow(M,8)) *  (  pow(I,2) - _deltaT2 * pow(Iprim,2)/pow(Mp,2) );
//  double lamda2 = pow(_s,2)-4*_s*pow(Mp,4);
//  return MT2/(16*pi*lamda2)*0.38941*pow(10,9);
//}

double func_Del_sq(double *x, double* par) {
  return _dsig_dDel_sq_a(x[0],par[0]);
}

double mirror(double t) {
  return msqTot - _ssss - t;
}

double func_Del_sq_mirror(double *x, double *par) {
  return _dsig_dDel_sq_a(mirror(x[0]), par[0]);
}

double func_s(double *x, double* par) {
  return _dsig_dDel_sq(par[0],x[0]);
}

void tda(int ip=0) {
  TF1* fDel_sq =  new TF1("Del_sq",func_Del_sq, _tmin[ip], _tmax[ip], 1);
  fDel_sq->SetParameter(0,_s[ip]);
  fDel_sq->Draw();
}

TF1* get_func(int ip, double _min, double _max) {
  TF1* fDel_sq =  new TF1(Form("Del_sq_%d",ip),func_Del_sq, _min, _max, 1);
  fDel_sq->SetParameter(0,_s[ip]);
  fDel_sq->SetLineWidth(3);
  fDel_sq->SetLineColor(1);
  return fDel_sq;
}

TF1* get_func(int ip) {
  TF1* fDel_sq =  new TF1(Form("Del_sq_%d",ip),func_Del_sq, _tmin[ip], _tmax[ip], 1);
  fDel_sq->SetParameter(0,_s[ip]);
  fDel_sq->SetLineWidth(3);
  return fDel_sq;
}

TF1* get_func_mirror(int ip=0) {
  TF1* fDel_sq_mir =  new TF1(Form("Del_sq_%d_mir",ip),func_Del_sq_mirror, 2*_tmin[ip]-_tmax[ip], _tmin[ip], 1);
  fDel_sq_mir->SetParameter(0,_s[ip]);
  fDel_sq_mir->SetLineWidth(3);
  return fDel_sq_mir;
}

void set_style(TH1* h, int col) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);

  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);

  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void tda_func(int ip=0, TPad *tp=0) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  //TString file_name[3] = {
  //  "/Users/tujuba/panda/work/jpsi_pi0/hists/note.aug.2015/eid90pct/pass0/anav2_jpsi_brem_plab5.5.root",
  //  "/Users/tujuba/panda/work/jpsi_pi0/hists/note.aug.2015/eid90pct/pass0/anav2_jpsi_brem_plab8.0.root",
  //  "/Users/tujuba/panda/work/jpsi_pi0/hists/note.aug.2015/eid90pct/pass0/anav2_jpsi_brem_plab12.0.root" };

  TString file_name[3] = {
    "/Users/tujuba/panda/work/jpsi_pi0/hists/pcm.jul.2015//anav2_jpsi_brem_plab5.5.root",
    "/Users/tujuba/panda/work/jpsi_pi0/hists/pcm.jul.2015//anav2_jpsi_brem_plab8.0.root",
    "/Users/tujuba/panda/work/jpsi_pi0/hists/pcm.jul.2015//anav2_jpsi_brem_plab12.0.root" };


  TCanvas *tc_gen_comp = new TCanvas(Form("gen_comp_p%d",ip),Form("gen_comp_p%d",ip),1200,500);

  _ssss = _s[ip];

  TFile *f = new TFile(file_name[ip]);

  TH1F* h = (TH1F*) f->Get("tu/httrumc");
  //h->SetTitle(Form("Expected signal rates, s=%5.2f GeV^{2};t[GeV^{2}]; dN_{sig}/dt [Counts/GeV^{2}]",_s[ip]));
  h->SetTitle(Form("s = %5.2f GeV^{2};t[GeV^{2}]; dN_{sig}/dt [Counts/GeV^{2}]",_s[ip]));
  set_style(h,1);

  if (ip==0) {h->GetXaxis()->SetRangeUser(-1.6,0.7);}

  double bw = h->GetXaxis()->GetBinWidth(3);
  //double evtCorr[3] = {32779.32/28040.0,   50142.17/43013.0,  51860.49/44425.};
  cout << "BW = " << bw << " Lumi= " << Lumi << " BW*Lumi*Br = " << bw*Lumi*Br << endl;

  if (tp!=0) tp->cd();

  h->Scale(1.0/bw/100); // WTF? 100
  //h->Scale(1.0/bw); // Why this is not working anymore?

  h->Draw();

  TBox *b_fwd = new TBox(_val_min[ip], 0, _val_max[ip], 1.05*h->GetMaximum());
  b_fwd->SetFillColor(kRed-10);
  b_fwd->Draw();

  TBox *b_bwd = new TBox(mirror(_val_min[ip]), 0, mirror(_val_max[ip]), 1.05*h->GetMaximum());
  b_bwd->SetFillColor(kCyan-10);
  b_bwd->Draw();

  h->Draw("same");

  TF1* func =  get_func(ip);
  func->DrawCopy("same");

  TF1* func_mirror =  get_func_mirror(ip);
  func_mirror->DrawCopy("same");

  TLegend *tl = new TLegend(0.4,0.55,0.7,0.85);
  tl->SetTextSize(0.05);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  tl->SetHeader("             L_{int} = 2 fb^{-1}");
  tl->AddEntry(h, "Generator Output","l");
  tl->AddEntry(func, "Model prediction","l");
  tl->Draw();

  TText *tt = new TText();
  tt->DrawText(_val_min[ip]+0.1, 0.1*h->GetMaximum(),"Fwd Kin");
  tt->DrawText(_val_min[ip]+0.1, 0.02*h->GetMaximum(),"Valid. range");

  //TText *tt = new TText();
  tt->DrawText(mirror(_val_max[ip])+0.1, 0.1*h->GetMaximum(),"Bwd Kin");
  tt->DrawText(mirror(_val_max[ip])+0.1, 0.02*h->GetMaximum(),"Valid. range");

  tc_gen_comp->Print(Form("gen_comp_%d.pdf",ip));

}

void draw_all() {
  TCanvas *tc = new TCanvas("tc","tc",1500,500);
  tc->Divide(3,1);
  tda_func(0, (TPad*)tc->GetPad(1));
  tda_func(1, (TPad*)tc->GetPad(2));
  tda_func(2, (TPad*)tc->GetPad(3));
}

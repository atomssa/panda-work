#include "TMath.h"

static const double pi = TMath::Pi();
static const double Mp = 0.938;
static const double Mp_sq = Mp*Mp;
static const double Mpi0 = 0.135;
static const double Mpi0_sq = Mpi0*Mpi0;
static const double Mjpsi = 3.097;
static const double Mjpsi2 = Mjpsi*Mjpsi;
static const double fpi = 0.093; // # GeV;
static const double fpsi = 0.413; //  #GeV;
static const double fN = 0.005; // #GeV;
static const double gpiNN = 13.0;
static const double alpha_s = 0.25;
static const double Mbar = 3.0; // #GeV;
static const double C = pow(4.0*pi*alpha_s,3)*fN*fN*fpsi*10.0/fpi/81.0;
//static const double M0 = 8776.88; //GeV
static const double jpsi_pbarp_width = 1.2e-6;
static const double M0 = sqrt(jpsi_pbarp_width*243.0*pi*pow(Mbar,5)/pow(pi*alpha_s,6)/1280.0/pow(fpsi,2)/pow(fN,4));
//print "M0 = %4.2f" % M0

//funciton definitions
double _lambda_sq(double x, double y, double z){
  return (x*x)+(y*y)+(z*z)-(2*x*y)-(2*x*z)-(2*y*z);
}

double _xi_a(double s){
  return Mbar*Mbar/(2.0*pow(s,2) - Mbar*Mbar);
}

double _xi(double Del_sq, double s){
  return (Mjpsi*Mjpsi - Del_sq - Mp_sq)/(2.0*pow(s,2) - Mjpsi*Mjpsi - Del_sq - 3.0*Mp_sq);
}

double _DelT_sq(double Del_sq, double s){
  return ((1-_xi(Del_sq,s))/(1+_xi(Del_sq,s))) * (Del_sq - 2.0*_xi(Del_sq,s)*( (Mp_sq/(1+_xi(Del_sq,s))) - (Mpi0_sq/(1-_xi(Del_sq,s))) ) );
}

double _Del_sq_max(double s){
  return 2.0*_xi_a(s)*(Mp_sq*(_xi_a(s)-1) + Mpi0_sq*(_xi_a(s)+1))/(pow(_xi_a(s),2)-1);
}

double _dsig_dDel_sq(double Del_sq,double s){
  const double DelT_sq = _DelT_sq(Del_sq,s);
  const double F1 = 1/16.0/pi/_lambda_sq(s,Mp_sq,Mp_sq);
  const double F2 = C*C*2.0*(1+_xi(Del_sq,s))/4.0/_xi(Del_sq,s)/pow(Mbar, 6);
  const double I = fpi*gpiNN*Mp*(1-_xi(Del_sq,s))*M0/(Del_sq - Mp_sq)/(1+_xi(Del_sq,s));
  const double Iprim = fpi*gpiNN*Mp*M0/(Del_sq - Mp_sq);
  return 4e6 * F1 * F2 * ( pow(I,2) - DelT_sq*pow(Iprim,2)/Mp_sq);
}

double func_Del_sq(double *x, double* par) {
  return _dsig_dDel_sq(x[0],par[0]);
}

double func_s(double *x, double* par) {
  return _dsig_dDel_sq(par[0],x[0]);
}

void tda() {
  TF1* fDel_sq =  new TF1("Del_sq",func_Del_sq, -0.5, 0.1, 1);
  fDel_sq->SetParameter(0,15.0);
  fDel_sq->Draw();
}

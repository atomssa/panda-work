//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: EvtPpbarJpsiPi0.cpp,v 1.2 2008/01/08 17:15:11 steinke Exp $
//
// Description:
//            Generator of pbar p -> jpsi + pi0
//
// Author List:
//
//
//------------------------------------------------------------------------
#include "EvtGenModels/EvtPpbarJpsiPi0.hh"

#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cmath> // for hypot
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtStdHep.hh"

//using std::fstream;
using std::string;
using namespace std;

EvtPpbarJpsiPi0::~EvtPpbarJpsiPi0()
{
}

std::string EvtPpbarJpsiPi0::getName()
{
  return "PpbarJpsiPi0";
}

EvtDecayBase* EvtPpbarJpsiPi0::clone()
{
  return new EvtPpbarJpsiPi0;
}

void EvtPpbarJpsiPi0::initProbMax()
{
  setProbMax(prob_max);
}

void EvtPpbarJpsiPi0::init()
{
  checkNDaug(2);
  checkNArg(1);
  _s = getArg(0);

  // constants for prob. calculation
  _Mp = 0.938;
  _Mpi0 = 0.0; // This is required to reproduce x-sect dists. on PLB paper
  _Mjpsi2 = 9.59;
  pi = 3.14159;
  fpi = 0.093;// GeV
  gpiNN = 13;
  M0 = 8776.88; //GeV
  fphi = 0.413; //GeV
  fN = 0.005;//GeV
  alpha_s = 0.252062;
  M = 3.;//GeV
  C = pow(4.*pi*alpha_s,3)*(fN*fN*fphi/fpi)*(10./81.);

  _Ebeam = (_s -2.*pow(_Mp,2))/(2.*_Mp);//5.59(s=12.25)   9.72(s=20)
  _Pbeam = sqrt(pow(_Ebeam,2)-pow(_Mp,2));//5.51(s=12.25)    9.68(s=20)

  // These only depend on _s, no need to calculate for every event
  //double _qsit = (_Mjpsi2 - _t - pow(_Mp,2))/(2.*_s-_Mjpsi2+_t-3.*pow(_Mp,2));
  _xi = M*M/(2*_s-M*M);
  _max = 2*_xi*(_Mp*_Mp*(_xi-1)+_Mpi0*_Mpi0*(_xi+1))/(_xi*_xi-1);

  boost_to_cm.set(0., 0., -_Pbeam/(_Ebeam + _Mp));

  prob_max = prob(0.7);
  cout << "EvtPpbarJpsiPi0::init() prob_max = " << prob_max << endl;
}

double EvtPpbarJpsiPi0::prob(double mand) {
  //if (mand > -1.0 and mand<_max) {
    double _deltaT2 = ((1-_xi)/(1+_xi)) * (mand - 2.*_xi * ( (pow(_Mp,2)/(1+_xi)) - (pow(_Mpi0,2)/(1-_xi)) ));
    double I = M0*fpi*gpiNN*_Mp*(1-_xi)/(mand-_Mp*_Mp)/(1+_xi);
    double Iprim = M0*fpi*gpiNN*_Mp/(mand-_Mp*_Mp);
    double MT2 = 0.25*C*C*(2*(_xi+1)/_xi/pow(M,8)) *  (  pow(I,2) - _deltaT2 * pow(Iprim,2)/pow(_Mp,2) );
    double lamda2 = pow(_s,2)-4*_s*pow(_Mp,4);
    return MT2/(16*pi*lamda2)*0.38941*pow(10,9);
  //} else {
  //  return 0;
  //}
}

void EvtPpbarJpsiPi0::decay(EvtParticle* p)
{
  p->initializePhaseSpace(getNDaug(),getDaugs());

  static const EvtId pi0ID = EvtPDL::getId("pi0");
  EvtParticle* pi0 =  p->getDaug(0);
  EvtParticle* jpsi =  p->getDaug(1);

  if(pi0->getId() != pi0ID ) {
    pi0= p->getDaug(1);
    jpsi= p->getDaug(0);
  }

  if(pi0->getId() != pi0ID )
    cout << "EvtPpbarJpsiPi0::decay():\n wrong id of produced particles!"<<endl;

  // calculate t and u of this decay
  EvtVector4R p4pi0 = pi0->getP4Lab();
  double _t = pow((_Ebeam-p4pi0.get(0)),2)-pow((0-p4pi0.get(1)),2)-pow((0-p4pi0.get(2)),2)-pow((_Pbeam-p4pi0.get(3)),2);
  EvtVector4R p4jpsi = jpsi->getP4Lab();
  double _u = pow((_Ebeam-p4jpsi.get(0)),2)-pow((0-p4jpsi.get(1)),2)-pow((0-p4jpsi.get(2)),2)-pow((_Pbeam-p4jpsi.get(3)),2);

  EvtVector4R p4pi0_cm(p4pi0);
  p4pi0_cm.applyBoostTo(boost_to_cm);

  double _prob = 0;
  if ( p4pi0_cm.get(3) > 0.0 ) {
    _prob = prob(_t);
  } else {
    _prob = prob(_u);
  }
  setProb(_prob);
  return;

}

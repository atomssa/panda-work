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
  setProbMax(500);
}

void EvtPpbarJpsiPi0::init()
{
  checkNDaug(2);
  checkNArg(1);
  _s = getArg(0); 
}

void EvtPpbarJpsiPi0::decay(EvtParticle* p)
{
  p->initializePhaseSpace(getNDaug(),getDaugs());
   
  static const EvtId pi0ID = EvtPDL::getId("pi0");
  EvtParticle* pi0 =  p->getDaug(0);
  
  if(pi0->getId() != pi0ID )
    pi0= p->getDaug(1);
  if(pi0->getId() != pi0ID )
    cout << "EvtPpbarJpsiPi0::decay():\n wrong id of produced particles!"<<endl;
  

  double _Mp = 0.938;
  double _Mpi0 = 0.135;
  double _Mjpsi2 = 9.59;

  
  double pi = 3.14159;
  double fpi = 0.093;// GeV
  double gpiNN = 13;
  double M0 = 8776.88; //GeV
  double fphi = 0.413; //GeV
  double fN = 0.005;//GeV
  double alpha_s = 0.25;
  double M = 3.;//GeV
  double C = pow(4.*pi*alpha_s,3)*(fN*fN*fphi/fpi)*(10./81.);

  
  // calculate t of this decay
  EvtVector4R p4 = pi0->getP4Lab();
  double _Ebeam = (_s -2.*pow(_Mp,2))/(2.*_Mp);//5.59(s=12.25)   9.72(s=20)
  double _Pbeam = sqrt(pow(_Ebeam,2)-pow(_Mp,2));//5.51(s=12.25)    9.68(s=20)
  double _t = pow((_Ebeam-p4.get(0)),2)-pow((0-p4.get(1)),2)-pow((0-p4.get(2)),2)-pow((_Pbeam-p4.get(3)),2); // t of the event

  //double _qsit = (_Mjpsi2 - _t - pow(_Mp,2))/(2.*_s-_Mjpsi2+_t-3.*pow(_Mp,2));
  double _qsi = M*M/(2*_s-M*M);
  double _deltaT2 = ((1-_qsi)/(1+_qsi))*(_t - 2.*_qsi*(pow(_Mp,2)/(1+_qsi)-pow(_Mpi0,2)/(1-_qsi)));
  double t_max = 2*_qsi*(_Mp*_Mp*(_qsi-1)+_Mpi0*_Mpi0*(_qsi+1))/(_qsi*_qsi-1);
  
  double prob;
  if (_t < -0.5 || _t >= t_max)
    prob =0;
  else
    {
      double I = M0*fpi*gpiNN*_Mp*(1-_qsi)/(_t-_Mp*_Mp)/(1+_qsi);
      double Iprim = M0*fpi*gpiNN*_Mp/(_t-_Mp*_Mp);
      double MT2 = 0.25*C*C*(2*(_qsi+1)/_qsi/pow(M,8))*(I*I-_deltaT2*Iprim*Iprim/pow(_Mp,2));
      double lamda2 = pow(_s,2)-4*_s*pow(_Mp,2);
      prob = MT2/(16*pi*lamda2)*0.38941*pow(10,9);
    }
  if (prob!=0)
    cout << "EvtPpbarJpsiPi0: t= " << _t << " prob= " << prob << endl;
  setProb(prob);
  return;

}


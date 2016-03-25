//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: EvtPpbarPi0Pi0.hh,v 1.2 2008/01/08 17:15:11 steinke Exp $
//
// Description:
//            Generator of pbar p -> Jpsi + Pi0
//
// Author List:
//
//------------------------------------------------------------------------
#ifndef EVTPPABRJPSIPI0_HH
#define EVTPPABRJPSIPI0_HH

#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtVector3R.hh"

class EvtParticle;
class EvtStdHep;

class EvtPpbarJpsiPi0:public EvtDecayProb{

public:
  EvtPpbarJpsiPi0() {};
  virtual ~EvtPpbarJpsiPi0();
  std::string getName();
  EvtDecayBase* clone();
  void init();
  void initProbMax();
  void decay(EvtParticle *p);

  double prob(double);
private:

  double _qsit;
  double _s;


  double _Mp;
  double _Mpi0;
  double _Mjpsi2;
  double pi;
  double fpi;
  double gpiNN;
  double M0;
  double fphi;
  double fN;
  double alpha_s;
  double M;
  double C;

  double _xi;
  double _max;

  double _Ebeam;
  double _Pbeam;

  double prob_max;

  EvtVector3R boost_to_cm;

};

#endif

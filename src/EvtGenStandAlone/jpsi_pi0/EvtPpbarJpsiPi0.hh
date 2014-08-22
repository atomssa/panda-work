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

class EvtParticle;

class EvtPpbarJpsiPi0:public EvtDecayProb{

public:
  EvtPpbarJpsiPi0() {};
  virtual ~EvtPpbarJpsiPi0();
  std::string getName();
  EvtDecayBase* clone();
  void init();
  void initProbMax();
  void decay(EvtParticle *p);

private:
  double _qsit;
  double _s;
};

#endif


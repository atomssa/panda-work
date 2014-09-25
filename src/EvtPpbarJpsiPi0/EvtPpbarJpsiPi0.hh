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

#include "TMath.h"

class EvtParticle;
class EvtStdHep;

class TTree;
class TParticle;
class TFile;
class TClonesArray;
class TLorentzVector;

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

  void set_p4_v4(EvtStdHep&, const int&, TLorentzVector&, TLorentzVector&);
  inline double _pow(double base, double exp) { return TMath::Power(base,exp); }

  double _qsit;
  double _s;

  EvtStdHep fEvtStdHep;  //! The decay tree
  TFile* fFile;
  TTree* fTree;
  TClonesArray* fEvt;
  int fNpart;


};

#endif

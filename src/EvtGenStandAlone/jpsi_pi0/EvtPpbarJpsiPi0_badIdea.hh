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

class EvtParticle;
class EvtStdHep;

// Temporary to store all events regardless of acceptance 
class TFile;
class TTree;
class TClonesArray;
class TParticle;
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
  
  double _qsit;
  double _s;

  // Temporary to store all events regardless of acceptance 
  TFile* fFile;
  TTree* fTree;
  TClonesArray* fEvt;
  TClonesArray* fEvtTlv;  
  float fProb;
  EvtStdHep fEvtStdHep;  //! The decay tree  

};

#endif


// -------------------------------------------------------------------------
// -----                 PndEvtGenStandAlone header file                  -----
// -----               Created 11/04/08  by M.Al-Turany              -----
// -------------------------------------------------------------------------

/** PndEvtGenStandAlone.h
 *@author M.Al-Turany <m.al-turany@gsi.de>
 *
 The PndEvtGenStandAlone generates EVT event using EvtGen
 and saves the event in a TTree of TClonesArray of TParticle

 Derived from FairTask so that it can be run in a fair task
**/

#ifndef PND_EVTSTANDALONE_H
#define PND_EVTSTANDALONE_H

#include "FairTask.h"

#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtId.hh"

#include <vector>

class EvtId;
class EvtStdHep;
class EvtParticle;

class EvtGen;
class TTree;
class TParticle;
class TFile;
class TClonesArray;
class TLorentzVector;

class PndEvtGenStandAlone : public FairTask
{

 public:

  /** Default constructor (should not be used) **/
  PndEvtGenStandAlone();

  PndEvtGenStandAlone(TString particle,TString decfile="",Double_t Mom=0, Long_t Seed=-1,TString defDECAY="",TString defPDL=""); // Mom>0 -> pbar Momentum; Mom<0 -> cms Energy

  /** Destructor **/
  virtual ~PndEvtGenStandAlone();

  /* virtual initializer inherited from superclass*/
  virtual InitStatus Init();
  // well no way to run as FairAna without fucking par file, running directly in macro
  void Initialize();

  /* virtual executor inherited from superclass*/
  virtual void Exec(Option_t* option);
  // well no way to run as FairAna without fucking par file, running directly in macro
  void Exec();
  void Exec(TClonesArray *);

  /* virtual finisher inherited from superclass*/
  virtual void Finish();
  void Finalize();
  void close_root_file();

  inline void SetVerbose(int v=1){fVerb=v;};

  virtual void set_evt_topo(const std::vector<Int_t> &arg) { ref_topo=arg; }

 private:

  /** Generate one event using EVT
   ** @param primGen  pointer to the FairPrimaryGenerator
   **/
  //Bool_t generate_event();
  Bool_t generate_event(TClonesArray *);

  void init_root_tree();
  void set_p4_v4(EvtStdHep&, const int &, TLorentzVector&, TLorentzVector&);
  void print_detail(EvtParticle *part);
  Bool_t verb_flag();
  Bool_t compare();

  TFile* fFile;
  TTree* fTree;
  TTree *fTreeUnfilt;
  TClonesArray* fEvt;
  Int_t fNpart;

  Int_t iEvt;

  Int_t fVerb;

  std::vector<Int_t> ref_topo;

  /**
   * P_lab(GeV/c)
   */
  Double_t fEnergy;	//! Energy of System
  Double_t fMomentum; 	//! Momentum of System

  EvtGen *fGenerator;	//! Pointer to the actual EvtGen
  EvtStdHep fEvtStdHep;  //! The decay tree
  EvtId PART;		    //! The mother particle

  ClassDef(PndEvtGenStandAlone,1);

};

#endif

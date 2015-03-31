#ifndef PNDPIDBREMCORRECTED4MOM_H
#define PNDPIDBREMCORRECTED4MOM_H
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// PndPidBremCorrected4Mom	                                                //
//                                                                      //
// Definition of the Panda pid candidate.	                        //
//                                                                      //
// Author: Klaus Goetzen, GSI, 12.06.08		                        //
// Copyright (C) 2008, GSI Darmstadt.					//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

#include <assert.h>

#include "FairRecoCandidate.h"
#include "TArrayI.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class PndPidBremCorrected4Mom : public FairMultiLinkedData
{

 public:

  PndPidBremCorrected4Mom();
  PndPidBremCorrected4Mom(TLorentzVector &p4);

  ~PndPidBremCorrected4Mom();

  TVector3  GetMomentum() const { return TVector3(fXmomentum, fYmomentum, fZmomentum); }
  Double_t  GetEnergy()   const { return fEnergy; }
  Int_t GetPidCandIdx() const { return fPidCandIdx; }
  const std::vector<Int_t> &GetPhiBumpList() {return fPhiBumpList; }
  const std::vector<Int_t> &GetSepBumpList() {return fSepBumpList; }

  void	SetMomentum(TVector3 &mom) { fXmomentum=mom.X(); fYmomentum=mom.Y(); fZmomentum=mom.Z(); }
  void	SetEnergy(Double_t en)     { fEnergy=(Float_t) en;}
  void  AddToPhiBumpList(Int_t idx) { fPhiBumpList.push_back(idx); }
  void  AddToSepBumpList(Int_t idx) { fSepBumpList.push_back(idx); }
  void  SetPidCandIdx(Int_t idx) { fPidCandIdx = idx; }

 protected:

  Double_t fXmomentum;		// The momentum in x
  Double_t fYmomentum;		// The momentum in y
  Double_t fZmomentum;		// The momentum in z
  Double_t fEnergy;

  Int_t fPidCandIdx;
  std::vector<Int_t> fPhiBumpList;
  std::vector<Int_t> fSepBumpList;

  ClassDef(PndPidBremCorrected4Mom,1)

};


#endif

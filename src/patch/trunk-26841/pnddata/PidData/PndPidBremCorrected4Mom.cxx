//////////////////////////////////////////////////////////////////////////
//                                                                      //
// PndPidBremCorrected4Mom	                                        //
//                                                                      //
// Container for Bremstrahlung radiaton corrected momentum              //
//                                                                      //
// Author: Klaus Goetzen, GSI, 12.06.08		                        //
// Copyright (C) 2008, GSI Darmstadt.		         		//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "PndPidBremCorrected4Mom.h"

PndPidBremCorrected4Mom::PndPidBremCorrected4Mom():
  fXmomentum(0.),
  fYmomentum(0.),
  fZmomentum(0.),
  fEnergy(0.),
  fPhiBumpList(),
  fSepBumpList()
{
  fPhiBumpList.clear();
  fSepBumpList.clear();
}

PndPidBremCorrected4Mom::PndPidBremCorrected4Mom(TLorentzVector &p4) :
  fXmomentum(p4.X()),
  fYmomentum(p4.Y()),
  fZmomentum(p4.Z()),
  fEnergy(p4.E()),
  fPhiBumpList(),
  fSepBumpList()
{
  fPhiBumpList.clear();
  fSepBumpList.clear();
}

PndPidBremCorrected4Mom::~PndPidBremCorrected4Mom() {

}

ClassImp(PndPidBremCorrected4Mom)

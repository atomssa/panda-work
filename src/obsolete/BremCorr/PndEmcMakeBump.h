//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:$
//
// Description:
//	Class Template
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
// Adapted for the PANDA experiment at GSI
//
// Author List:
//	Xiaorong Shi            Lawrence Livermore National Lab
//	Steve Playfer           University of Edinburgh
//	Stephen Gowdy           University of Edinburgh
//      Phil Strother           Imperial College
//
// Dima Melnychuk, adaption for PANDA
// Modified:
// M. Babai
//------------------------------------------------------------------------
#ifndef PNDEMCMAKEBUMP_H
#define PNDEMCMAKEBUMP_H

#include "FairTask.h"

class TClonesArray;
class TObjectArray;

class PndEmcCluster;
class PndEmcDigi;
class PndEmcSharedDigi;
class PndEmcBump;
class PndEmcTwoCoordIndex;

class PndEmcMakeBump  : public FairTask
{
	
	
public:

  // Constructors
  PndEmcMakeBump(Int_t verbose=0, Bool_t storebumps=kTRUE);

  // Destructor
  virtual ~PndEmcMakeBump( );
  
  /** Virtual method Init **/
  virtual InitStatus Init();
  

  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt);
  
  void SetStorageOfData(Bool_t val); // Method to specify whether bumps are stored or not.
  
 protected:
  
 private:
	
  /** Get parameter containers **/
  virtual void SetParContainers();
  
  /** Verbosity level **/
  Int_t fVerbose;
  
  Bool_t fPersistance;
  
  static Int_t fEventCounter;
  
  ClassDef(PndEmcMakeBump,1);
};
#endif //PNDEMCMAKEBUMP_HH

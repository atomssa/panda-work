//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:$
//
// Description:
//      This module takes Clusters (Connected Regions) and slits them
//      up into Bumps. There are defined by local maxima
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
// Adapted for the PANDA experiment at GSI
//
// Author List:
//	Stephen J. Gowdy           University of Edinburgh
//      Phil Strother              Imperial College
//
// Dima Melnychuk, adaption for PANDA
// Modified:
// M. Babai
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "PndEmcMakeBump.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "PndEmcDataTypes.h"

#include "PndEmcBump.h"
#include "PndEmcCluster.h"
#include "PndEmcDigi.h"
#include "PndEmcSharedDigi.h"
#include "PndEmc2DLocMaxFinder.h"
#include "PndEmcExpClusterSplitter.h"
#include "PndEmcPhiBumpSplitter.h"
#include "PndEmcTwoCoordIndex.h"

#include "TClonesArray.h"
#include "TROOT.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>

using std::endl;
using std::cout;

Int_t PndEmcMakeBump::fEventCounter=1;

//----------------
// Constructors --
//----------------
PndEmcMakeBump::PndEmcMakeBump(Int_t verbose, Bool_t persistance):
FairTask("EMC Bump splitting Task"), fVerbose(verbose), fPersistance(persistance)
{
  this->Add(new PndEmc2DLocMaxFinder());
  this->Add(new PndEmcExpClusterSplitter());
  for (int i=0; i<8; ++i)
    this->Add(new PndEmcPhiBumpSplitter(i));

  TList* thistasks = this->GetListOfTasks();
  for(Int_t i=0;i<thistasks->GetEntries();i++)
  {
    ((FairTask*)thistasks->At(i))->SetVerbose(fVerbose);
  }

	SetStorageOfData(fPersistance);
}

void PndEmcMakeBump::SetStorageOfData(Bool_t val)
{
  fPersistance=val;
  TList* thistasks = this->GetListOfTasks();
  ((PndEmc2DLocMaxFinder*)thistasks->At(0))->SetStorageOfData(fPersistance);
  ((PndEmcExpClusterSplitter*)thistasks->At(1))->SetStorageOfData(fPersistance);
  for (int i=0; i<8; ++i) {
    ((PndEmcPhiBumpSplitter*)thistasks->At(2+i))->SetStorageOfData(fPersistance);
  }
  return;
}

//--------------
// Destructor --
//--------------
PndEmcMakeBump::~PndEmcMakeBump()
{
}

// -----   Public method Init   -------------------------------
InitStatus PndEmcMakeBump::Init() {
  return kSUCCESS;
}

void PndEmcMakeBump::Exec(Option_t* opt)
{
	if (fVerbose>0)
		std::cout<<"***************** PndEmcMakeBump, event: "
		<<fEventCounter<<" **************"<<endl;
	fEventCounter++;
	return;
}

void PndEmcMakeBump::SetParContainers() {}

ClassImp(PndEmcMakeBump)

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:$
//
// Description:
//
// Environment:
//
//
//
// Author List:
//
//
// Copyright Information:
//	Copyright (C) 1997               Imperial College
//
// Modified:
//
//------------------------------------------------------------------------
#pragma once
#ifndef PNDEMCPHIBUMPSPLITTER_H
#define PNDEMCPHIBUMPSPLITTER_H

// Path of file:
// ----- $pandaroot/emc/EmcReco

//---------------
// C++ Headers --
//---------------
#include <vector>
#include <map>

#include "FairTask.h"
#include "TObject.h"
#include "PndEmcDataTypes.h"
//#include "PndEmcDigiCalibrator.h"
//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class PndEmcCluster;
class PndEmcBump;
class PndEmcDigi;
class PndEmcTwoCoordIndex;

class PndEmcGeoPar;
class PndEmcDigiPar;
class PndEmcRecoPar;

class PndEmcSharedDigi;

//		---------------------
// 		-- Class Interface --
//		---------------------

class PndEmcPhiBumpSplitter: public FairTask
{

 public:

	PndEmcPhiBumpSplitter(Int_t verbose=0);

	// Destructor
	virtual ~PndEmcPhiBumpSplitter( );

	// Methods
	/** Virtual method Init **/
	virtual InitStatus Init();

	/** Virtual method Exec **/
	virtual void Exec(Option_t* opt);

	void SetStorageOfData(Bool_t p = kTRUE) {fPersistance=p;};
	PndEmcBump* AddPhiBump();
	//PndEmcSharedDigi* AddPhiBumpSharedDigi(PndEmcDigi*, Double_t);

	virtual void FinishTask();

 private:
	/** Input array of PndEmcClusters **/
	TClonesArray* fDigiArray;
	TClonesArray* fClusterArray;

	/** Output array of PndEmcBumps **/
	TClonesArray* fPhiBumpArray;
	//TClonesArray* fPhiBumpSharedDigiArray;


	PndEmcGeoPar*     fGeoPar;       /** Geometry parameter container **/
	PndEmcDigiPar*    fDigiPar;      /** Digitisation parameter container **/
	PndEmcRecoPar*    fRecoPar;      /** Reconstruction parameter container **/
	/** Get parameter containers **/
	virtual void SetParContainers();

	std::vector<Double_t> fClusterPosParam;

	Bool_t fPersistance; // switch to turn on/off storing the arrays to a file
	// Data members

	/** Verbosity level **/
	Int_t fVerbose;

        PndEmcPhiBumpSplitter(const  PndEmcPhiBumpSplitter& L);
        PndEmcPhiBumpSplitter& operator= (const  PndEmcPhiBumpSplitter&) {return *this;};

	ClassDef(PndEmcPhiBumpSplitter,1);

//added for time information

	//PndEmcDigiCalibrator digiCalibrator;
	//Int_t HowManyDidis;

};
#endif

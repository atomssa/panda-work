//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:$
//
// Description:
//	Class PndEmcExpClusterSplitter.
//      Implementation of ClusterSplitter which splits
//      on the basis of exponential distance from the bump centroid.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Adapted for the PANDA experiment at GSI
//
// Author List:
//      Phil Strother               
//
// Copyright Information:
//	Copyright (C) 1997               Imperial College
//
// Modified:
// M. Babai
//------------------------------------------------------------------------


//-----------------------
// This Class's Header --
//-----------------------
#include "PndEmcPhiBumpSplitter.h"

//---------------
// C++ Headers --
//---------------
//#include <vector>
//#include <set>
//#include <map>
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

#include "PndEmcClusterProperties.h"
#include "PndEmcXClMoments.h"

#include "PndEmcRecoPar.h"
#include "PndEmcGeoPar.h"
#include "PndEmcDigiPar.h"
#include "PndEmcStructure.h"
#include "PndEmcMapper.h"

#include "PndDetectorList.h"
#include "PndEmcTwoCoordIndex.h"
#include "PndEmcBump.h"
#include "PndEmcCluster.h"
#include "PndEmcDigi.h"
#include "PndEmcSharedDigi.h"
#include "PndEmcXtal.h"
#include "PndEmcDataTypes.h"
using std::endl;

//----------------
// Constructors --
//----------------

PndEmcPhiBumpSplitter::PndEmcPhiBumpSplitter(Int_t verbose):
fDigiArray(0), fClusterArray(0), fPhiBumpArray(0), fPhiBumpSharedDigiArray(0), fGeoPar(new PndEmcGeoPar()), fDigiPar(new PndEmcDigiPar()), fRecoPar(new PndEmcRecoPar()), fPersistance(kTRUE), fMaxIterations(0), fCentroidShift(0), fClusterPosParam(), fVerbose(verbose)
{
  fClusterPosParam.clear();
}

//--------------
// Destructor --
//--------------

PndEmcPhiBumpSplitter::~PndEmcPhiBumpSplitter()
{
// 	delete fGeoPar;
// 	delete fDigiPar;
// 	delete fRecoPar;
}

InitStatus PndEmcPhiBumpSplitter::Init() {
  
  // Get RootManager
  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- PndEmcMakeBump::Init: "
	 << "RootManager not instantiated!" << endl;
    return kFATAL;
  }
	
  // Geometry loading
  fGeoPar->InitEmcMapper();
  PndEmcStructure::Instance();
	
	// Get input array
	fDigiArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcDigi"));
	if ( ! fDigiArray ) {
		cout << "-W- PndEmcMakeCluster::Init: "
		<< "No PndEmcDigi array!" << endl;
		return kERROR;
	}
  
  fClusterArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcCluster"));
  if ( ! fClusterArray ) {
    cout << "-W- PndEmcMakeBump::Init: "
	 << "No PndEmcCluster array!" << endl;
    return kERROR;
  }



  // Set minimum SharedDigi energy to 20keV.
	
	if (!strcmp(fRecoPar->GetEmcClusterPosMethod(),"lilo"))
	{
		cout<<"Lilo cluster position method"<<endl;
		fClusterPosParam.push_back(fRecoPar->GetOffsetParmA());
		fClusterPosParam.push_back(fRecoPar->GetOffsetParmB());
		fClusterPosParam.push_back(fRecoPar->GetOffsetParmC());
	}	

  // Create and register output array
  fPhiBumpArray = new TClonesArray("PndEmcPhiBump");
  ioman->Register("EmcPhiBump","Emc",fPhiBumpArray,fPersistance);
  
  fPhiBumpSharedDigiArray = new TClonesArray("PndEmcPhiBumpSharedDigi");
  ioman->Register("EmcPhiBumpSharedDigi","Emc",fPhiBumpSharedDigiArray,fPersistance);


	HowManyDidis = 0;

  cout << "-I- PndEmcExpClusterSplitter: Intialization successfull" << endl;

}

void PndEmcPhiBumpSplitter::Exec(Option_t* opt) 
{

	PndEmcMapper *fEmcMap=PndEmcMapper::Instance();
	// Reset output array
	if ( ! fPhiBumpArray ) Fatal("Exec", "No Phi-Bump Array");
	fPhiBumpArray->Delete();

	int nClusters = fClusterArray->GetEntriesFast();
	
	for (Int_t iCluster = 0; iCluster < nClusters; iCluster++){
		PndEmcCluster* theCluster = (PndEmcCluster*) fClusterArray->At(iCluster);

		int numberOfPhiBumps = -1;
		numberOfBumps = (theCluster->LocalMaxMap()).size();
	
		if (numberOfBumps<=1 || numberOfBumps >= fMaxBumps){
			// Limit the max number of bumps in the cluster to 8 (default)
			// in this case, we clearly have a cluster, but no bumps to speak of.  
			// Make 1 bump with weights all equal to 1
		
			std::map<Int_t, Int_t>::const_iterator theDigiIterator;
			
			PndEmcBump* theNewBump = AddBump(); 
			theNewBump->MadeFrom(iCluster);
			theNewBump->SetLink(FairLink("EmcCluster", iCluster));

			for(theDigiIterator = theCluster->MemberDigiMap().begin();
			theDigiIterator != theCluster->MemberDigiMap().end(); ++theDigiIterator){
				PndEmcDigi *theDigi = (PndEmcDigi *) fDigiArray->At(theDigiIterator->second);
				PndEmcSharedDigi* sharedDigi=AddSharedDigi(theDigi, 1.0);
				Int_t iSharedDigi=fSharedDigiArray->GetEntriesFast()-1;
				theNewBump->addDigi(fSharedDigiArray,iSharedDigi);
			}
		}
		else {
			std::map<Int_t,Int_t> theMaximaDigis=theCluster->LocalMaxMap();
			std::map<Int_t,Int_t>::iterator theMaximaDigisIterator;
			std::map<PndEmcTwoCoordIndex*, TVector3*> theCentroidPoints;
			std::map<PndEmcTwoCoordIndex*, TVector3*> theMaximaPoints;
			std::map<PndEmcTwoCoordIndex*, PndEmcBump*> theIndexedBumps;
			std::map<PndEmcTwoCoordIndex*, TVector3*> theAllDigiPoints;
			
			std::map<Int_t, Int_t> theDigiDict = theCluster->MemberDigiMap();
			
			double totalEnergy=0;
			
			for(theMaximaDigisIterator = theMaximaDigis.begin(); 
			theMaximaDigisIterator != theMaximaDigis.end();
			++theMaximaDigisIterator){
				PndEmcDigi *theMaxDigi = (PndEmcDigi *) fDigiArray->At(theMaximaDigisIterator->second);
				
				Int_t detId = theMaximaDigisIterator->first;
				PndEmcTwoCoordIndex *theTCI = fEmcMap->GetTCI(detId);
				
				totalEnergy += theMaxDigi->GetEnergy();
				
				TVector3 *digiLocation = new TVector3(theMaxDigi->where());
				TVector3 *sameLocation = new TVector3(theMaxDigi->where());
				
				theMaximaPoints.insert(std::map<PndEmcTwoCoordIndex*, TVector3*>::value_type(theTCI, digiLocation));
				theCentroidPoints.insert(std::map<PndEmcTwoCoordIndex*, TVector3*>::value_type(theTCI, sameLocation));
			}
	
			std::map<Int_t,Int_t>::iterator theAllDigisIterator;		
			// This loop works out the location of all the digis in the cluster 
		
			for( theAllDigisIterator = theDigiDict.begin(); theAllDigisIterator != theDigiDict.end();
					++theAllDigisIterator){
				Int_t detId = theAllDigisIterator->first;
				PndEmcTwoCoordIndex *theTCI = fEmcMap->GetTCI(detId);
				PndEmcDigi *theDigi = (PndEmcDigi *) fDigiArray->At(theAllDigisIterator->second);
				TVector3 *digiLocation = new TVector3(theDigi->where());
				theAllDigiPoints.insert(std::map<PndEmcTwoCoordIndex*, TVector3*>::value_type( theTCI, digiLocation ));
			}

			theMaximaDigisIterator = theMaximaDigis.begin();

			// Now we can create the EmcBumps

			// The algorithm is as follows: We will index each bump by its
			// maximum digi's PndEmcTwoCoordIndex.  We will set up a list of
			// bump centroids which to start with will be synonymous with the
			// location of the maxima.  We then apportion a weight to each
			// digi, according to its distance from the centroids.  We then
			// construct the bumps according to these weights, which will
			// presumably give a different set of centroids.  This is repeated
			// until the centroids are static within tolerance, or we reach
			// the maximum number of iterations.

			Int_t iterations = 0;

			Double_t averageCentroidShift;

			do {
				if (fVerbose>=3){
					std::cout<<"iteration No "<<iterations<<std::endl;
				}
				averageCentroidShift=0.0;

				// First clean up the old bumps
				std::map<PndEmcTwoCoordIndex*, PndEmcBump*>::iterator theBumpKiller = theIndexedBumps.begin();
				while(theBumpKiller != theIndexedBumps.end()){
					PndEmcBump *theBump = theBumpKiller->second;
					delete theBump;
					++theBumpKiller;
				}
				theIndexedBumps.clear();

				// Then loop over all the maxima and assign weights accordingly
				for (theMaximaDigisIterator = theMaximaDigis.begin();
						theMaximaDigisIterator != theMaximaDigis.end();
						++theMaximaDigisIterator) {
					Int_t detId = theMaximaDigisIterator->first;
					PndEmcTwoCoordIndex *theCurrentMaximaTCI = fEmcMap->GetTCI(detId);

					if (fVerbose>=3){
						std::cout<<"***************** current maximum: theta = "<<theCurrentMaximaTCI->XCoord()
							<<", phi = "<<theCurrentMaximaTCI->YCoord()<<"*********"<<std::endl;
					}

					// Create the bump which will correspond to this digi maxima
					PndEmcBump* theNewBump = new PndEmcBump();
					theNewBump->MadeFrom(iCluster);
					theIndexedBumps.insert(std::map<PndEmcTwoCoordIndex*,
							PndEmcBump*>::value_type(theCurrentMaximaTCI, theNewBump));


					// Now we will look over all the digis and add each of them
					// to this Bump with an appropriate weight

					for (theAllDigisIterator = theDigiDict.begin();
							theAllDigisIterator != theDigiDict.end();++theAllDigisIterator) {
						PndEmcDigi *theCurrentDigi = (PndEmcDigi *) fDigiArray->At(theAllDigisIterator->second);
						PndEmcTwoCoordIndex *theCurrentTCI = theCurrentDigi->GetTCI();

						Double_t weight;

						// We are on the first pass and the digi is not a local max, or we are not on the
						// first pass.   Assign a weight according to the distance from the centroid position.

						Double_t myEnergy = 0;
						Double_t myDistance = 0;

						// Now share the digi out according to its distance from the maxima, 
						// and the maxima energies

						Double_t totalDistanceEnergy=0;

						//Moliere Radius for Shashlyk is different
						Double_t MoliereRadius;
						if(theCurrentDigi->GetModule() == 5)
							MoliereRadius = fMoliereRadiusShashlyk;
						else
							MoliereRadius = fMoliereRadius;

						std::map<PndEmcTwoCoordIndex*,TVector3*>::iterator theMaxPointsIterator;

						for(theMaxPointsIterator = theCentroidPoints.begin(); theMaxPointsIterator != theCentroidPoints.end();++theMaxPointsIterator){
							PndEmcTwoCoordIndex *theMaxPointsTCI =theMaxPointsIterator->first; 

							TVector3 *theMaxPoint = theMaxPointsIterator->second;
							TVector3 *theCurrentDigiPoint = theAllDigiPoints.find(theCurrentTCI)->second;

							Double_t theDistance;

							// This next bit just checks to see if the maxima point in
							// hand is the same as the crystal from which we are
							// trying to find distance - just an FP trap really.

							if ((*theCurrentTCI)==(*theMaxPointsTCI)){
								theDistance=0.0;
							} else {
								TVector3 distance( *theMaxPoint - *theCurrentDigiPoint);

								theDistance = distance.Mag();
							}

							if (*theCurrentMaximaTCI == *(theMaxPointsTCI)){
								// i.e. the maximum we are trying to find the distance from is 
								// the one for which we are currently trying to make a bump
								myDistance = theDistance;
								Int_t iCurentMaxDigi = (theDigiDict.find(theMaxPointsTCI->Index()))->second;
								myEnergy = ((PndEmcDigi *) fDigiArray->At(iCurentMaxDigi))->GetEnergy();
							}

							Int_t iMaxPoint = (theDigiDict.find(theMaxPointsTCI->Index()))->second;
							totalDistanceEnergy += ((PndEmcDigi *) fDigiArray->At(iMaxPoint))->GetEnergy() *
								exp(-fExponentialConstant * theDistance/MoliereRadius);
						}

						if(totalDistanceEnergy > 0.0)
							weight = myEnergy*exp(-fExponentialConstant* 
									myDistance/MoliereRadius) / ( totalDistanceEnergy);
						else 
							weight=0;

						if (fVerbose>=3){
							std::cout<<"\t digi theta = "<<theCurrentDigi->GetTCI()->XCoord()
								<<", phi = "<<theCurrentDigi->GetTCI()->YCoord()<<std::endl;
							std::cout<<"energy = "<<theCurrentDigi->GetEnergy()<<", weight = "<< weight<<std::endl;
						}
						PndEmcSharedDigi* sharedDigi= AddSharedDigi( theCurrentDigi, weight );

						if (fVerbose>=3){
							std::cout<<"shared digi energy = "<<sharedDigi->GetEnergy()<<std::endl;
						}

						Int_t iSharedDigi=fSharedDigiArray->GetEntriesFast()-1;
						if( sharedDigi->GetEnergy() > fMinDigiEnergy){
							theNewBump->addDigi( fSharedDigiArray, iSharedDigi );
						} else {
							fSharedDigiArray->RemoveAt(iSharedDigi);
							fSharedDigiArray->Compress();
						}
					}

					// Compute the shift of the centroid we have just calculated
					TVector3 *theOldCentroid = theCentroidPoints.find(theCurrentMaximaTCI)->second;

					PndEmcClusterProperties clusterProperties(*theNewBump, fSharedDigiArray);

					TVector3 newbumppos =  clusterProperties.Where(fRecoPar->GetEmcClusterPosMethod(), fClusterPosParam);
					theNewBump->SetPosition(newbumppos);
					TVector3 centroidShift(*theOldCentroid - newbumppos);
					averageCentroidShift+=centroidShift.Mag();
				}

				averageCentroidShift/=(Double_t)numberOfBumps;

				// Put the new centroids in the list of centroid points,
				// remembering to delete the old ones.
				std::map<PndEmcTwoCoordIndex*, TVector3*>::iterator 
					theCentroidPointsIterator = theCentroidPoints.begin();
				for(theCentroidPointsIterator = theCentroidPoints.begin();
						theCentroidPointsIterator != theCentroidPoints.end(); 
						++theCentroidPointsIterator) {
					delete theCentroidPointsIterator->second;
				}
				theCentroidPoints.clear();

				std::map<PndEmcTwoCoordIndex*, PndEmcBump*>::iterator theIndexedBumpsIterator;
				for(theIndexedBumpsIterator = theIndexedBumps.begin();
						theIndexedBumpsIterator != theIndexedBumps.end(); ++theIndexedBumpsIterator){
					TVector3 *theNewCentroid = new TVector3((theIndexedBumpsIterator->second)->where());
					theCentroidPoints.insert(std::map<PndEmcTwoCoordIndex*, TVector3*>::value_type(theIndexedBumpsIterator->first,theNewCentroid));
				}

				iterations++;

			} while (iterations < fMaxIterations && averageCentroidShift > fCentroidShift);
			//End of do loop

			// Finally append the new bumps to the TClonesArray.
			std::map<PndEmcTwoCoordIndex*, PndEmcBump*>::iterator theBumpsIterator;
			PndEmcBump *theBump;

			for(theBumpsIterator = theIndexedBumps.begin(); theBumpsIterator != theIndexedBumps.end();++theBumpsIterator){
				theBump = theBumpsIterator->second;
				Int_t size_ba = fBumpArray->GetEntriesFast();
				PndEmcBump* theNextBump = new((*fBumpArray)[size_ba]) PndEmcBump(*(theBump));
				if (fVerbose>0)
					std::cout << "Bump Created!" << std::endl;
				theNextBump->SetLink(FairLink("EmcCluster", iCluster));
			}

			std::map<PndEmcTwoCoordIndex*,TVector3*>::iterator theGrimReaper = theMaximaPoints.begin();
			for(theGrimReaper = theMaximaPoints.begin();
					theGrimReaper != theMaximaPoints.end();
					++theGrimReaper){
				delete theGrimReaper->second;
			}
			theMaximaPoints.clear();

			for(theGrimReaper  = theAllDigiPoints.begin();
					theGrimReaper != theAllDigiPoints.end();
					++theGrimReaper){
				delete theGrimReaper->second;
			}
			theAllDigiPoints.clear();

			for(theGrimReaper = theCentroidPoints.begin(); 
					theGrimReaper != theCentroidPoints.end();
					++theGrimReaper){
				delete theGrimReaper->second;
			}
			theCentroidPoints.clear();
		}

		Int_t nBumps = (theCluster->LocalMaxMap()).size();
		theCluster->SetNBumps(nBumps);
	}

	// At that moment internal state fEnergy and fWhere of Clusters are
	// not initialized, the following make it possible to see energy and
	// position from output root file
	Int_t nBump = fBumpArray->GetEntriesFast();

	Double_t CalibTimeOfaDigi, fTimeError;
	Double_t WeightedFactor1(0.), NormWeightedFactor1(0.), AverageTime1(0.);
	//Double_t WeightedFactor2(0.), NormWeightedFactor2(0.), AverageTime2(0.);
	//Double_t WeightedFactor3(0.), NormWeightedFactor3(0.), AverageTime3(0.);
	for (Int_t i=0; i<nBump; i++){
		PndEmcBump *tmpbump = (PndEmcBump*) fBumpArray->At(i);
		PndEmcClusterProperties clusterProperties(*tmpbump, fSharedDigiArray);

		if (!tmpbump->IsEnergyValid())
			tmpbump->SetEnergy(clusterProperties.Energy());
		if (!tmpbump->IsPositionValid())
			tmpbump->SetPosition(clusterProperties.Where(fRecoPar->GetEmcClusterPosMethod(), fClusterPosParam));
		PndEmcXClMoments xClMoments(*tmpbump, fSharedDigiArray);
		tmpbump->SetZ20(xClMoments.AbsZernikeMoment(2, 0, 15));
		tmpbump->SetZ53(xClMoments.AbsZernikeMoment(5, 3, 15));
		tmpbump->SetLatMom(xClMoments.Lat());

		//for time information
		WeightedFactor1 = 0.;//= WeightedFactor2 = WeightedFactor3 = 0.;
		AverageTime1 = 0.;//AverageTime2 = AverageTime3 = 0.;
		Double_t fMaxDigiEnergy = -1.;
		const std::vector<Int_t>& listOfDigi = tmpbump->DigiList();
		for(Int_t id=0;id <listOfDigi.size();++id){
			PndEmcDigi* theDigi = (PndEmcDigi*)fSharedDigiArray->At(listOfDigi[id]);
			//CalibTimeOfaDigi = digiCalibrator.CalibrationEvtTimeByDigi(theDigi, kFALSE);
			fTimeError = digiCalibrator.GetTimeResolutionOfDigi(theDigi);
			WeightedFactor1 += 1./fTimeError/fTimeError;
			//WeightedFactor2 += theDigi->GetEnergy()*1./fTimeError/fTimeError;
			//WeightedFactor3 += theDigi->GetEnergy();
			if(theDigi->GetEnergy() > fMaxDigiEnergy){
				fMaxDigiEnergy = theDigi->GetEnergy();
				tmpbump->SetEventNo(theDigi->fEvtNo);
				//tmpbump->SetTimeStamp1(theDigi->GetTimeStamp());
				//tmpbump->fSeedPosition = theDigi->where();
			}
		}
		for(Int_t id=0;id <listOfDigi.size();++id){
			PndEmcDigi* theDigi = (PndEmcDigi*)fSharedDigiArray->At(listOfDigi[id]);
			CalibTimeOfaDigi = digiCalibrator.CalibrationEvtTimeByDigi(theDigi, kFALSE);
			fTimeError = digiCalibrator.GetTimeResolutionOfDigi(theDigi);
			NormWeightedFactor1 = 1./fTimeError/fTimeError;
			NormWeightedFactor1 /= WeightedFactor1;
			AverageTime1 += NormWeightedFactor1*theDigi->GetTimeStamp();

			//NormWeightedFactor2 = theDigi->GetEnergy()*1./fTimeError/fTimeError;
			//NormWeightedFactor2 /= WeightedFactor2;
			//AverageTime2 += NormWeightedFactor2*theDigi->GetTimeStamp();

			//NormWeightedFactor3 = theDigi->GetEnergy();
			//NormWeightedFactor3 /= WeightedFactor3;
			//AverageTime3 += NormWeightedFactor3*theDigi->GetTimeStamp();
		}
		tmpbump->SetTimeStamp(AverageTime1);
		HowManyDidis += tmpbump->NumberOfDigis();
		//tmpbump->SetTimeStamp1(AverageTime1);
		//tmpbump->SetTimeStamp2(AverageTime2);
		//tmpbump->SetTimeStamp3(AverageTime3);
		//end for time information
	}

	if (fVerbose>=1){
		std::cout<<"PndEmcExpClusterSplitter:: Number of clusters = "<<nClusters<<std::endl;
		std::cout<<"PndEmcExpClusterSplitter:: Number of bumps = "<<nBump<<std::endl;
	}

}

PndEmcBump* PndEmcExpClusterSplitter::AddBump(){
	TClonesArray& clref = *fBumpArray;
	Int_t size = clref.GetEntriesFast();
	return new(clref[size]) PndEmcBump();
}

PndEmcSharedDigi* PndEmcExpClusterSplitter::AddSharedDigi(PndEmcDigi* digi, Double_t weight){
	TClonesArray& clref = *fSharedDigiArray;
	Int_t size = clref.GetEntriesFast();
	return new(clref[size]) PndEmcSharedDigi(*digi, weight);
}

void PndEmcExpClusterSplitter::FinishTask()
{
	cout<<"================================================="<<endl;
	cout<<"PndEmcExpClusterSplitter::FinishTask"<<endl;
	cout<<"================================================="<<endl;
	cout<<"read digis #"<<HowManyDidis<<endl;
}

void PndEmcExpClusterSplitter::SetParContainers() {

	// Get run and runtime database
	FairRun* run = FairRun::Instance();
	if ( ! run ) Fatal("SetParContainers", "No analysis run");

	FairRuntimeDb* db = run->GetRuntimeDb();
	if ( ! db ) Fatal("SetParContainers", "No runtime database");
	// Get Emc digitisation parameter container
	fGeoPar = (PndEmcGeoPar*) db->getContainer("PndEmcGeoPar");
	// Get Emc digitisation parameter container
	fDigiPar = (PndEmcDigiPar*) db->getContainer("PndEmcDigiPar");
	// Get Emc reconstruction parameter container
	fRecoPar = (PndEmcRecoPar*) db->getContainer("PndEmcRecoPar");
}

ClassImp(PndEmcExpClusterSplitter)

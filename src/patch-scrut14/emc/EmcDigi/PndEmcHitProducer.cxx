/////////////////////////////////////////////////////////////
//
//  PndEmcHitProducer
//
//  Filler of PndEmcHit
//
//  Created 14/08/06  by S.Spataro
//
/////////////////////////////////////////////////////////////// 

#include "PndEmcHitProducer.h"

#include "PndEmcStructure.h"
#include "PndEmcHit.h"
#include "PndEmcPoint.h"
#include "PndEmcGeoPar.h"
#include "PndEmcDigiPar.h"		
#include "PndEmcDigiNonuniformityPar.h"		
#include "PndMCTrack.h"

#include "PndEmcXtal.h"

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairDetector.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

#include "TClonesArray.h"
#include "TROOT.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TString.h"

static Int_t HowManyPoints = 0;
static Int_t HowManyHitsAll = 0;
static Int_t HowManyHitsAboveThreshold = 0;

// -----   Default constructor   -------------------------------------------
PndEmcHitProducer::PndEmcHitProducer() :
	FairTask("Ideal EMC hit Producer"),
	fUse_nonuniformity(0), fNonuniformityFile(""), fPointArray(), fMCTrackArray(), fHitArray(), fVolumeArray(), fMapVersion(0), emcX(), emcY(), emcZ(), fEmcStr(), fMapper(), fDigiPar(), fGeoPar(), fNonuniformityPar(), fStoreHits(kTRUE), fEnergyThreshold(0)
{
	fNonuniformityFile=gSystem->Getenv("VMCWORKDIR");
	fNonuniformityFile+="/input/EmcDigiNoniformityPars.root";
}
// -------------------------------------------------------------------------

PndEmcHitProducer::PndEmcHitProducer(Bool_t val) :
	FairTask("Ideal EMC hit Producer"),
	fUse_nonuniformity(0), fNonuniformityFile(""), fPointArray(), fMCTrackArray(), fHitArray(), fVolumeArray(), fMapVersion(0), emcX(), emcY(), emcZ(), fEmcStr(), fMapper(), fDigiPar(), fGeoPar(), fNonuniformityPar(), fStoreHits(val), fEnergyThreshold(0)
{ 
	fNonuniformityFile=gSystem->Getenv("VMCWORKDIR");
	fNonuniformityFile+="/input/EmcDigiNoniformityPars.root";
}

// -----   Destructor   ----------------------------------------------------
PndEmcHitProducer::~PndEmcHitProducer() { delete fEmcStr;}
// -------------------------------------------------------------------------


// -----   Public method Init   --------------------------------------------
InitStatus PndEmcHitProducer::Init(){

	cout << " INITIALIZATION *********************" << endl;

	//FairDetector::Initialize();
	//FairRun* sim = FairRun::Instance();
	//FairRuntimeDb* rtdb=sim->GetRuntimeDb();

	// Get RootManager
	FairRootManager* ioman = FairRootManager::Instance();
	if (! ioman ){
		cout << "-E- PndEmcHitProducer::Init: "
			<< "RootManager not instantiated!" << endl;
		return kFATAL;
	}

	// Get input array
	fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
	if (! fPointArray ){
		cout << "-W- PndEmcHitProducer::Init: "
			<< "No EmcPoint array!" << endl;
		return kERROR;
	}

	// Get input array
	fMCTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
	if (! fMCTrackArray ){
		cout << "-W- PndEmcMakeCluster::Init: "
			<< "No MCTrack array! Needed for MC Truth" << endl;
		//return kERROR;
	}

	// Create and register output array
	fHitArray = new TClonesArray("PndEmcHit");

	ioman->Register("EmcHit","Emc",fHitArray,fStoreHits);

	fGeoPar->InitEmcMapper();
	fMapper=PndEmcMapper::Instance();
	fEmcStr=PndEmcStructure::Instance();

	emcX=fEmcStr->GetEmcX();
	emcY=fEmcStr->GetEmcY();
	emcZ=fEmcStr->GetEmcZ();;

	fEnergyThreshold =fDigiPar->GetEnergyHitThreshold();
	fUse_nonuniformity = fDigiPar->GetUse_nonuniformity();

	if(fUse_nonuniformity){
		cout << "-I- PndEmcHitProducer: Using nonuniform lightoutput" << endl;
	}
	if(fUse_nonuniformity && fNonuniformityFile.Length()>0){
		TFile *nonuniformityfile = new TFile(fNonuniformityFile);
		if(nonuniformityfile==NULL){
			cout << "-E- PndEmcHitProducer: Could not open file " << fNonuniformityFile.Data() << " for Nonuniformity Information" << endl;
		} else {
			PndEmcDigiNonuniParObject *parObject;
			nonuniformityfile->GetObject("PndEmcDigiNonuniParObject",parObject);
			if(parObject == NULL){
				cout << "-E- PndEmcHitProducer: Could not get Nonuniformity information from file " << fNonuniformityFile.Data() << endl;
			} else {
				fNonuniformityPar->SetNonuniParObject(parObject);
			}
		}
	}


	printf("HitProducer has EnergyHitThreshold of %f GeV and Use_nonuniformity %i\n", fEnergyThreshold, fUse_nonuniformity);
	cout << "-I- PndEmcHitProducer: Intialization successfull" << endl;

	return kSUCCESS;
}

void PndEmcHitProducer::SetParContainers(){

	// Get run and runtime database
	FairRun* run = FairRun::Instance();
	if (! run ) Fatal("SetParContainers", "No analysis run");

	FairRuntimeDb* db = run->GetRuntimeDb();
	if (! db ) Fatal("SetParContainers", "No runtime database");

	// Get Emc geometry parameter container
	fGeoPar = (PndEmcGeoPar*) db->getContainer("PndEmcGeoPar");

	// Get Emc digitisation parameter container
	fDigiPar = (PndEmcDigiPar*) db->getContainer("PndEmcDigiPar");

	fNonuniformityPar = (PndEmcDigiNonuniformityPar*) db->getContainer("PndEmcDigiNonuniformityPar");

	fDigiPar->setChanged();
	fDigiPar->setInputVersion(run->GetRunId(),1); 

	fNonuniformityPar->setChanged();
	fNonuniformityPar->setInputVersion(run->GetRunId(),1);
}

// -------------------------------------------------------------------------

// Helper function, does not depend on class, identical to the one in PndEmcMakeCluster
void PndEmcHitProducer::cleansortmclist( std::vector <Int_t> &newlist,TClonesArray* mcTrackArray)
{
  std::vector <Int_t> tmplist;
  std::vector <Int_t> tmplist2;
  // Sort list...
  std::sort( newlist.begin(), newlist.end());
  // and copy every id only once (even so it might be in the list several times)
  std::unique_copy( newlist.begin(), newlist.end(), std::back_inserter( tmplist ) );

  newlist.clear();
  std::unique_copy( tmplist.begin(), tmplist.end(), std::back_inserter( newlist) );
}

//// Helper function, does not depend on class, identical to the one in PndEmcMakeCluster
//void PndEmcHitProducer::cleansortmclist( std::vector <Int_t> &newlist,TClonesArray* mcTrackArray)
//{
//	std::vector <Int_t> tmplist;
//	std::vector <Int_t> tmplist2;
//	// Sort list...
//	std::sort( newlist.begin(), newlist.end());
//	// and copy every id only once (even so it might be in the list several times)
//	std::unique_copy( newlist.begin(), newlist.end(), std::back_inserter( tmplist ) );
//
//	// Now check if mother or (grand)^x-mother are already in the list
//	// (which means i am a secondary)... if so, remove myself
//	for(Int_t j=tmplist.size()-1; j>=0; j--){
//		bool flag = false;
//		PndMCTrack *pt;
//		//pt=((PndMCTrack*)mcTrackArray->At(tmplist[j]));
//		//if(pt->GetMotherID()<0) { 
//		//	tmplist2.push_back(tmplist[j]);
//		//	continue;
//		//}
//		Int_t id = tmplist[j];
//		if(id < 0) {
//			tmplist2.push_back(id);
//			continue;
//		}
//		while(!flag){
//			pt=((PndMCTrack*)mcTrackArray->At(id));
//			//id=pt->GetMotherID();
//			if(pt->GetMotherID()<0) {
//				tmplist2.push_back(id);
//				break;
//			}
//			id = pt->GetMotherID();
//			//pt=(PndMCTrack*)mcTrackArray->At(id);
//
//			for(Int_t k=j-1; k>=0; k--){
//				if(tmplist[k]==id){
//					tmplist.erase(tmplist.begin()+j);
//					flag=true;
//					break;
//				}
//			}
//		}
//	}
//	newlist.clear();
//	std::unique_copy( tmplist2.begin(), tmplist2.end(), std::back_inserter( newlist) );
//}

// -----   Public method Exec   --------------------------------------------
void PndEmcHitProducer::Exec(Option_t* opt)
{  
	cout << " POINT EXECUTION *********************" << endl;
	// Reset output array
	if (! fHitArray ) Fatal("Exec", "No DigiArray");

	fHitArray->Delete();

	// Declare some variables
	//PndEmcPoint* point  = NULL;
	Int_t DetId;

	fTrackEnergy.clear();
	fTrackTime.clear();
	fTrackMcTruth.clear();
	fPointMatch.clear();

	map<Int_t, Float_t>::const_iterator p;

	std::vector<PndEmcPoint*> fPointList;// to pass to EmcHit
	const PndEmcTciXtalMap &XtalMap = fEmcStr->GetTciXtalMap();
	TVector3 frontvec;
	TVector3 normvec;
	TVector3 pointvec;
	TVector3 distvec;
	Double_t zpos;
	Double_t energyscalefactor=1.0;
	Double_t c[3];
	PndEmcXtal *tmpXtal;
	PndEmcTwoCoordIndex *tmpTCI;
	// Loop over EmcPoints
	Int_t nPoints = fPointArray->GetEntriesFast();

	Double_t point_time = 0.00;
	//------- init containers --- 

	for (Int_t iPoint = 0; iPoint < nPoints; iPoint++){
		PndEmcPoint* point  = (PndEmcPoint*) fPointArray->At(iPoint);
		fTrackEnergy[point->GetDetectorID()] = 0.00;
		fTrackTime  [point->GetDetectorID()] = std::numeric_limits<float>::max();
	}

	//----------------------------

	//Int_t Counter[3] = {0};
	HowManyPoints += nPoints;

	for (Int_t iPoint=0; iPoint<nPoints; iPoint++)
	{
		PndEmcPoint* point  = (PndEmcPoint*) fPointArray->At(iPoint);
		DetId = point->GetDetectorID();

		if(point->GetEnergyLoss() == 0 ) continue;
		//if(point->GetTrackID() < 0)
		//	std::cout<<"negative track id #"<<point->GetTrackID()<<std::endl;

		//std::cout<<"point belongs to track #"<<point->GetTrackID()<<std::endl;
		//if(point->GetTrackID() == 0) ++Counter[0];
		//if(point->GetTrackID() == 1) ++Counter[1];
		//if(point->GetTrackID() == 2) ++Counter[2];

		if(fUse_nonuniformity !=0 ){
			//light output is z-dependent, so calculate z

			tmpTCI = fMapper->GetTCI(DetId);
			if(tmpTCI == NULL){
				printf("no TCI found for DetectorID %d\n",DetId);
				continue;
			}
			tmpXtal =XtalMap.find(tmpTCI)->second;
			point->Position(pointvec);
			frontvec = tmpXtal->frontCentre();
			normvec = tmpXtal->normalToFrontFace();
			distvec = pointvec-frontvec;
			zpos = distvec.Dot(normvec);
			fNonuniformityPar->GetNonuniformityParameters(DetId,c);
			energyscalefactor=c[0]+zpos*(c[1]+zpos*c[2]);
			fTrackEnergy[DetId] += point->GetEnergyLoss() * energyscalefactor;
			fPointMatch[DetId].push_back(iPoint);
			//        printf("point with detID %d has z Position %f and energyloss %f scaled with %f\n",DetId,zpos, point->GetEnergyLoss(),energyscalefactor);	
			//        printf("front is at x: %f y: %f z: %f\n", frontvec.X(),frontvec.Y(),frontvec.Z());
		} else {
			fTrackEnergy[DetId] += point->GetEnergyLoss();
			fPointMatch[DetId].push_back(iPoint);
			//        printf("point with detID %d has z Position %f and energyloss %f not scaled\n",DetId,zpos, point->GetEnergyLoss());	
		}
		point_time=point->GetTime();

		if (point_time < fTrackTime[point->GetDetectorID()]){
			fTrackTime[point->GetDetectorID()] = point_time;
		}

		// Check and save MC truth information
		// Eloss==0 tracks are only stored in point, if track is entering detector from outside
		// and thats what we are interested in...
		//std::cout<<"track id #"<<point->GetTrackID()<<", Energyloss #"<<point->GetEnergyLoss()<<endl;
		if(point->GetEnergyLoss() >0)
			(fTrackMcTruth[point->GetDetectorID()]).push_back(point->GetTrackID());
		//if(point->GetEnergyLoss() == 0 ){
		//	cout << "ELoss== 0 : " <<point->GetEnergyLoss()<<", ID "<<
		//		point->GetTrackID()<<","<<point->GetDetectorID()<<","<<point->GetXPad()<<","<<point->GetYPad()<<endl;
		//}else{
		//	cout << "ELoss>0 : " <<point->GetEnergyLoss()<<", ID "<<
		//		point->GetTrackID()<<","<<point->GetDetectorID()<<","<<point->GetXPad()<<","<<point->GetYPad()<<endl;
		//}
	}

	// Loop over EmcPoint

	// Loop to register EmcHit
	Int_t idx = 0;
	for( p = fTrackEnergy.begin(); p != fTrackEnergy.end(); ++p){
		++HowManyHitsAll;
		++idx;
		if ((*p).second > fEnergyThreshold){
			++ HowManyHitsAboveThreshold;
			// Check and save MC truth information B.S.
			// remove MC Truth particles which are not needed (eg grand^x-daugherts)
			if( fMCTrackArray){
				std::vector<Int_t>& plist = fTrackMcTruth[(*p).first];
				cleansortmclist(fTrackMcTruth[(*p).first],fMCTrackArray);
				//std::cout<<"The "<<(idx)<<" hit produced by track #";
				//for(Int_t ip=0;ip<plist.size();++ip){
				//	cout<<plist[ip]<<", ";
				//}
				//cout<<endl;
			}
			//std::vector<Int_t>& track = fTrackMcTruth[(*p).first];
			//std::cout<<track.size()<<" tracks contributes energy to this hit, ";
			//for(size_t i=0;i<track.size();++i){
			//	if(i < track.size() - 1){
			//		std::cout<<track[i]<<", ";
			//	}else{
			//		std::cout<<track[i];
			//	}
			//}
			//std::cout<<std::endl;
			PndEmcHit* myHit = AddHit(1, (*p).first, (*p).second, fTrackTime[(*p).first], fTrackMcTruth[(*p).first]);
			myHit->AddLinks(FairMultiLinkedData("EmcPoint", fPointMatch[p->first]));
		}
	}
	//check
	//Int_t nTrack = fMCTrackArray->GetEntriesFast();
	//for(Int_t itrack = 0; itrack < nTrack; ++itrack){
	//	PndMCTrack* pt1 =((PndMCTrack*)fMCTrackArray->At(itrack));
	//	if(pt1->IsGeneratorCreated()){
	//		pt1->Print(itrack);
	//	}
	//}


}

// -------------------------------------------------------------------------
// -----   Private method AddDigi   --------------------------------------------
PndEmcHit* PndEmcHitProducer::AddHit(Int_t trackID,Int_t detID, Float_t energy,
		Float_t time, std::vector <Int_t> &mctruth)
{
	// It fills the PndEmcHit category

	//cout << "PndEmcHitProducer: track " << trackID << " evt " << eventID
	//<< " sec " << sec << " plane " << pla << " strip " << strip << "box
	//" << box << " tube " << tub << endl;
	TClonesArray& clref = *fHitArray;
	Int_t size = clref.GetEntriesFast();
	return new(clref[size]) PndEmcHit(trackID, detID, energy, time, emcX[detID], 
			emcY[detID], emcZ[detID], mctruth);
}
// ----

void PndEmcHitProducer::SetStorageOfData(Bool_t val)
{
	fStoreHits=val;
	return;
}
void PndEmcHitProducer::FinishTask()
{
	std::cout<<"========================================================="<<std::endl;
	std::cout<<"PndEmcHitProducer::FinishTask"<<std::endl;
	std::cout<<"*********************************************************"<<std::endl;
	std::cout<<"Read points # "<<HowManyPoints<<std::endl;
	std::cout<<"Produc hits# "<<HowManyHitsAll<<", threshold# "<<fEnergyThreshold<<std::endl;
	std::cout<<"Hits above threshhod#"<<HowManyHitsAboveThreshold<<std::endl;
	std::cout<<"*********************************************************"<<std::endl;
}

ClassImp(PndEmcHitProducer)

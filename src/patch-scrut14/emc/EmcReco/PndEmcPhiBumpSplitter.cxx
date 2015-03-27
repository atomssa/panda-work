//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:$
//
// Description:
//	Class PndEmcPhiBumpSplitter
//      Implementation of PhiBumpSplitter which splits clusters based on
//      local maxima in the Phi Direction for use with Bremstrahlung correction
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
// Binsong Ma, Ermias Atomssa
//------------------------------------------------------------------------

// Path of file:
// ----- $pandaroot/emc/EmcReco

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
  fDigiArray(0), fClusterArray(0), fPhiBumpArray(0), fGeoPar(new PndEmcGeoPar()), fDigiPar(new PndEmcDigiPar()), fRecoPar(new PndEmcRecoPar()), fPersistance(kTRUE), fClusterPosParam(), fVerbose(verbose)
{
  fClusterPosParam.clear();
}

//--------------
// Destructor --
//--------------

PndEmcPhiBumpSplitter::~PndEmcPhiBumpSplitter()
{
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
  fPhiBumpArray = new TClonesArray("PndEmcBump");
  ioman->Register("EmcPhiBump","Emc",fPhiBumpArray,fPersistance);

  cout << "-I- PndEmcPhiBumpSplitter: Intialization successfull" << endl;

}

void PndEmcPhiBumpSplitter::Exec(Option_t* opt)
{

  PndEmcMapper *fEmcMap=PndEmcMapper::Instance();

  // Reset output array
  if ( ! fPhiBumpArray ) Fatal("Exec", "No Phi-Bump Array");
  fPhiBumpArray->Delete();

  // loop on each cluster. For each cluster there can be any number of phi_bumps (atleast one)
  int nClusters = fClusterArray->GetEntriesFast();
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++){

    PndEmcCluster* theCluster = (PndEmcCluster*) fClusterArray->At(iCluster);

    const Int_t TotNumOfHitPhi = 160;

    std::vector<double> phi_bump(TotNumOfHitPhi,0);

    Int_t digiSize = theCluster->DigiList().size();
    std::vector<Int_t> digiList = theCluster->DigiList();
    for (Int_t i_digi = 0; i_digi<digiSize; ++i_digi)
      {
	PndEmcDigi *emcDigi = (PndEmcDigi *) fDigiArray->At(digiList[i_digi]);
	Double_t emcDigiPhi = emcDigi->GetPhi()*TMath::RadToDeg();
	Double_t emcDigiEnergy = emcDigi->GetEnergy();
	if ( fabs(emcDigiPhi)<=180 ) {
	  Int_t iEmcDigiPhi = int( (emcDigiPhi + 180.) * TotNumOfHitPhi / 360. );
	  phi_bump.at(iEmcDigiPhi) += emcDigiEnergy;
	}
      }

    std::vector<double> vPhiList;
    std::vector<double> vDepoEnergyList;
    std::vector<int> vGapSizeList;

    int i_phi_prev = 0;
    for (Int_t i_phi=0; i_phi < TotNumOfHitPhi; ++i_phi) {
      Double_t BinValue = phi_bump.at(i_phi);
      if (BinValue != 0)
	{
	  vPhiList.push_back(-180. + (160.*(0.5+i_phi)/360.));
	  vDepoEnergyList.push_back(BinValue);
	  vGapSizeList.push_back(i_phi - i_phi_prev );
	  i_phi_prev = i_phi;
	}
    }

    // Find start bin number for the phi projection of cluster.
    Int_t StartIndex = 0;
    for (Int_t i = 0;i < vGapSizeList.size();i++)
      {
	if (vGapSizeList.at(i) > vGapSizeList.at(StartIndex)) StartIndex = i;
      }

    // Rotate to the begining of the cluster in case cluster folds over to index=0
    std::rotate(vDepoEnergyList.begin(),vDepoEnergyList.begin()+StartIndex,vDepoEnergyList.end());
    vDepoEnergyList.push_back(0);
    vDepoEnergyList.push_back(0);
    std::rotate(vDepoEnergyList.begin(),vDepoEnergyList.begin()+(vDepoEnergyList.size()-1),vDepoEnergyList.end());
    // Do the same for the phi positions
    std::rotate(vPhiList.begin(),vPhiList.begin()+StartIndex,vPhiList.end());
    vPhiList.push_back(0);
    vPhiList.push_back(0);
    std::rotate(vPhiList.begin(),vPhiList.begin()+(vPhiList.size()-1),vPhiList.end());


    // Loop through deposited energy vector and classify bins
    std::vector<int> Type;
    std::vector<double> Weight;
    Weight.push_back(0);
    Type.push_back(-3);
    for (Int_t n_sel = 1; n_sel < vDepoEnergyList.size()-1; n_sel++)
      {
	if (vDepoEnergyList.at(n_sel-1) < vDepoEnergyList.at(n_sel) && vDepoEnergyList.at(n_sel) < vDepoEnergyList.at(n_sel+1) ) Type.push_back(1);
	else if(vDepoEnergyList.at(n_sel-1) < vDepoEnergyList.at(n_sel) && vDepoEnergyList.at(n_sel) > vDepoEnergyList.at(n_sel+1) )
	  {
	    Type.push_back(0);
	    Weight.push_back(vDepoEnergyList.at(n_sel));
	  }
	else if(vDepoEnergyList.at(n_sel-1) > vDepoEnergyList.at(n_sel) && vDepoEnergyList.at(n_sel) > vDepoEnergyList.at(n_sel+1) ) Type.push_back(-1);
	else if(vDepoEnergyList.at(n_sel-1) > vDepoEnergyList.at(n_sel) && vDepoEnergyList.at(n_sel) < vDepoEnergyList.at(n_sel+1) ) Type.push_back(-2);
      }
    Weight.push_back(0);

    std::vector<double> enePhiBump, phiPhiBump;
    int ValleyIndex = 0;
    int iWeight = 0;
    for (Int_t n_sel = 1; n_sel < vDepoEnergyList.size()-1; n_sel++)
      {
	if (Type.at(n_sel) == -2 || n_sel == vDepoEnergyList.size()-2)
	  {
	    iWeight++;

	    const double _eneRightEdge = vDepoEnergyList.at(ValleyIndex)*(Weight.at(iWeight)/(Weight.at(iWeight)+Weight.at(iWeight-1)));
	    double _enePhiBump = _eneRightEdge;
	    double _phiPhiBump = vPhiList.at(ValleyIndex)*_eneRightEdge;
	    for(Int_t p = ValleyIndex+1;p < n_sel;p++) {
	      _enePhiBump += vDepoEnergyList.at(p);
	      _phiPhiBump += vPhiList.at(p)*vDepoEnergyList.at(p);
	    }
	    const double _eneLeftEdge = vDepoEnergyList.at(n_sel)*(Weight.at(iWeight)/(Weight.at(iWeight)+Weight.at(iWeight+1) ));
	    _enePhiBump += _eneLeftEdge;
	    _phiPhiBump += vPhiList.at(n_sel)*_eneLeftEdge;
	    _phiPhiBump /= _enePhiBump;
	    enePhiBump.push_back(_enePhiBump);
	    phiPhiBump.push_back(_phiPhiBump);
	    ValleyIndex = n_sel;
	  }
      }

    TVector3 posClust = theCluster->position();
    for (int i_phibump=0; i_phibump<enePhiBump.size(); ++i_phibump) {
      PndEmcBump* theNewPhiBump = AddPhiBump();
      theNewPhiBump->MadeFrom(iCluster);
      theNewPhiBump->SetLink(FairLink("EmcCluster", iCluster));
      theNewPhiBump->SetEnergy(enePhiBump.at(i_phibump));
      TVector3 posPhiBump;
      posPhiBump.SetMagThetaPhi(posClust.Mag(),posClust.Theta(),phiPhiBump.at(i_phibump));
      theNewPhiBump->SetPosition(posPhiBump);
    }

  }

}

PndEmcBump* PndEmcPhiBumpSplitter::AddPhiBump(){
  TClonesArray& clref = *fPhiBumpArray;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) PndEmcBump();
}

void PndEmcPhiBumpSplitter::FinishTask()
{
  cout<<"================================================="<<endl;
  cout<<"PndEmcPhiBumpSplitter::FinishTask"<<endl;
  cout<<"================================================="<<endl;
}

void PndEmcPhiBumpSplitter::SetParContainers() {

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

ClassImp(PndEmcPhiBumpSplitter)

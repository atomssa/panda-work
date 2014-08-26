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

#include "TH1.h"

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

  h_phi_bump = new TH1F("h_phi_bump","phi bumps",160,-180,180);
  
  HowManyDidis = 0;

  cout << "-I- PndEmcPhiBumpSplitter: Intialization successfull" << endl;

}

void PndEmcPhiBumpSplitter::Exec(Option_t* opt) 
{

  PndEmcMapper *fEmcMap=PndEmcMapper::Instance();

  // Reset output array
  if ( ! fPhiBumpArray ) Fatal("Exec", "No Phi-Bump Array");
  fPhiBumpArray->Delete();

  Double_t DepoEnergyList[100];
  Int_t DepoENoOfBinList[100];
  Int_t GapSizeList[100];
  Int_t DepoEListIndex = 1;
  
  // loop on each cluster. For each cluster there can be any number of phi_bumps (atleast one) 
  int nClusters = fClusterArray->GetEntriesFast();
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++){  

    PndEmcCluster* theCluster = (PndEmcCluster*) fClusterArray->At(iCluster);
        
    Int_t digiSize = theCluster->DigiList().size();
    std::vector<Int_t> digiList = theCluster->DigiList();

    for (Int_t i_digi = 0; i_digi<digiSize; ++i_digi)
      {
	PndEmcDigi *emcDigi = (PndEmcDigi *) fDigiArray->At(digiList[i_digi]);
	Double_t emcDigiPhi = emcDigi->GetPhi()*TMath::RadToDeg();
	Double_t emcDigiEnergy = emcDigi->GetEnergy();
	h_phi_bump->Fill(EleEmcDigiPhi,EleEmcDigiEnergy); 
      }

    Int_t TotNumOfHitPhi = h_phi_bump->GetNbinsX();

    Double_t DepoEnergyList[100] = {0.};
    Int_t DepoENoOfBinList[100] = {0};
    Int_t GapSizeList[100] = {0};
    Int_t DepoEListIndex = 1;

    DepoEnergyList[0] = 0;
    DepoENoOfBinList[0] = -10;
    GapSizeList[0] = -1;
     
    for (Int_t i_phi = 1; i_phi <= TotNumOfHitPhi; i_phi++)
      {
	Double_t BinValue = h_phi_bump->GetBinContent(i_phi);
	if (BinValue != 0) 
	  {
	    DepoEnergyList[DepoEListIndex] = BinValue;
	    DepoENoOfBinList[DepoEListIndex] = i_phi;
	    GapSizeList[DepoEListIndex] = DepoENoOfBinList[DepoEListIndex] - DepoENoOfBinList[DepoEListIndex-1];
	    DepoEListIndex++;
	  }
      }
    DepoEnergyList[DepoEListIndex] = 0;
    DepoENoOfBinList[DepoEListIndex] = -1;
    GapSizeList[DepoEListIndex] = -1;

    Int_t StartIndice = 0;
    for (Int_t i = 1;i < DepoEListIndex;i++)
      {
	if (GapSizeList[i] > GapSizeList[StartIndice]) StartIndice = i;     
      }

    if(StartIndice > 1)
      {   
	Int_t i_corr = 0;
	Double_t DepoEnergyListCorr[100];
	DepoEnergyListCorr[0] = DepoEnergyList[0];
	for (Int_t i2 = StartIndice; i2 < DepoEListIndex;i2++)
	  {
	    i_corr++;
	    DepoEnergyListCorr[i_corr] = DepoEnergyList[i2];
	  }
	for(Int_t i2 = 1;i2 < StartIndice;i2++)
	  {
	    i_corr++;
	    DepoEnergyListCorr[i_corr] = DepoEnergyList[i2];
	  }
	DepoEnergyListCorr[DepoEListIndex] = DepoEnergyList[DepoEListIndex];
      } //StartIndice
    else if (StartIndice = 1)
      {
	for (Int_t i2 = 0; i2 <= DepoEListIndex ; i2++) DepoEnergyListCorr[i2] = DepoEnergyList[i2];
      }

     
    Int_t Case[100];
    Double_t enePhiBump[100] = {0.};
    Double_t Poid[100] = {0.};


    //std::vector<int> case;
    //std::vector<double> enePhiBump;
    //std::vector<double> phiPhiBump;
    //std::vector<double> weightPhiBump;
    //Int_t IndiceVally = 0;

     
    Int_t PhiBumpIndex = 0, PoidIndex = 0;
    Case[0] = -3;
    Poid[PoidIndex] = 0;
    Int_t IndiceVally = 0;
     
    for (Int_t n_sel = 1; n_sel < DepoEListIndex; n_sel++)
      { 
	if (DepoEnergyListCorr[n_sel-1] < DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] < DepoEnergyListCorr[n_sel+1] ) Case[n_sel] = 1;
	else if(DepoEnergyListCorr[n_sel-1] < DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] > DepoEnergyListCorr[n_sel+1] )
	  {
	    PoidIndex++;
	    Case[n_sel] = 0;
	    Poid[PoidIndex] = DepoEnergyListCorr[n_sel];
	  }
	else if(DepoEnergyListCorr[n_sel-1] > DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] > DepoEnergyListCorr[n_sel+1] ) Case[n_sel] = -1;
	else if(DepoEnergyListCorr[n_sel-1] > DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] < DepoEnergyListCorr[n_sel+1] ) Case[n_sel] = -2;
      }
    Poid[PoidIndex+1] = 0;

    Int_t iPoid = 0;

    for (Int_t n_sel = 1; n_sel < DepoEListIndex; n_sel++)
      {
	if (Case[n_sel] == -2 || n_sel == DepoEListIndex-1)
	  {
	    iPoid++;
	    enePhiBump[PhiBumpIndex] = DepoEnergyListCorr[IndiceVally]*(Poid[iPoid]/(Poid[iPoid]+Poid[iPoid-1]));
	    for(Int_t p = IndiceVally;p < n_sel;p++) enePhiBump[PhiBumpIndex] += DepoEnergyListCorr[p];
	    enePhiBump[PhiBumpIndex] += DepoEnergyListCorr[n_sel]*(Poid[iPoid]/(Poid[iPoid]+Poid[iPoid+1]));
	    PhiBumpIndex++;
	    IndiceVally = n_sel;
	  }
      }

    //for (int i=0; i<enePhiBump.size(); ++i) {
    //  PndEmcBump* theNewPhiBump = AddPhiBump(); 
    //  theNewPhiBump->MadeFrom(iCluster);
    //  theNewPhiBump->SetLink(FairLink("EmcCluster", iCluster));  
    //  theNewBump->SetEnergy(enePhiBump.at(i));
    //  TVector3 pos = calcPosition(phiPhiBump.at(i));
    //  theNewBump->SetPhi(phiPhiBump.at(i));       
    //}
       
    for (int i_phibump=0; i_phibump<=PhiBumpIndex; ++i_phibump) {
      PndEmcBump* theNewPhiBump = AddPhiBump(); 
      theNewPhiBump->MadeFrom(iCluster);
      theNewPhiBump->SetLink(FairLink("EmcCluster", iCluster));  
      theNewPhiBump->SetEnergy(enePhiBump[i_phibump]);
      //  TVector3 pos = calcPosition(phiPhiBump.at(i));
      //  theNewBump->SetPhi(phiPhiBump.at(i));              
    }
    
  }
	
}

PndEmcBump* PndEmcPhiBumpSplitter::AddPhiBump(){
  TClonesArray& clref = *fPhiBumpArray;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) PndEmcBump();
}

PndEmcSharedDigi* PndEmcPhiBumpSplitter::AddPhiBumpSharedDigi(PndEmcDigi* digi, Double_t weight){
  TClonesArray& clref = *fPhiBumpSharedDigiArray;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) PndEmcSharedDigi(*digi, weight);
}

void PndEmcPhiBumpSplitter::FinishTask()
{
  cout<<"================================================="<<endl;
  cout<<"PndEmcPhiBumpSplitter::FinishTask"<<endl;
  cout<<"================================================="<<endl;
  cout<<"read digis #"<<HowManyDidis<<endl;
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

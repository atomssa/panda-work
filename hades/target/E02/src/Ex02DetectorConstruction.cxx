// $Id: Ex02DetectorConstruction.cxx 375 2008-12-23 09:50:46Z ivana $

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007, 2008 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: vmc@pcroot.cern.ch
//-------------------------------------------------

/// \file Ex02DetectorConstruction.cxx 
/// \brief Implementation of the Ex02DetectorConstruction class 
///
/// Geant4 ExampleN02 adapted to Virtual Monte Carlo \n
/// Id: Ex02DetectorConstruction.cc,v 1.11 2002/01/09 17:24:09 ranjard Exp \n
/// GEANT4 tag Name: geant4-04-00-patch-02 
///
/// \date 21/04/2002
/// \author I. Hrivnacova; IPN, Orsay
 
#include <iostream>

#include <TVirtualMC.h>
#include <TGeoManager.h>

#include "Ex02DetectorConstruction.h"
#include "Ex02ChamberParameterisation.h"

/// \cond CLASSIMP
ClassImp(Ex02DetectorConstruction)
/// \endcond

using namespace std;

//_____________________________________________________________________________
Ex02DetectorConstruction::Ex02DetectorConstruction()
  : TObject(),
    fWorldLength(0.),  
    fTargetLength(0.), 
    fTrackerLength(0.),
    fNofChambers(0) ,  
    fChamberWidth(0.),  
    fChamberSpacing(0.),
    fImedAir(0),
    fImedPb(0),
    fImedXe(0)
{
/// Default constuctor

  //fpMagField = new ExN02MagneticField();

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  fNofChambers = 5;
  fChamberWidth =   20;
  fChamberSpacing = 80;
  
  fTrackerLength = (fNofChambers+1)*fChamberSpacing; // Full length of Tracker
  fTargetLength  = 5.0;                              // Full length of Target
  fWorldLength= 1.2 *(fTargetLength+fTrackerLength);
}

//_____________________________________________________________________________
Ex02DetectorConstruction::~Ex02DetectorConstruction()
{
/// Destructor

  //delete fpMagField;
}

//_____________________________________________________________________________
void Ex02DetectorConstruction::ConstructMaterials()
{
/// Construct materials using TGeo modeller

//--------- Material definition ---------

  // Create Root geometry manager 
  new TGeoManager("E02_geometry", "E02 VMC example geometry");
   
  Double_t a;        // Mass of a mole in g/mole   
  Double_t z;        // Atomic number
  Double_t density;  // Material density in g/cm3
 
  // Elements
  
  TGeoElement* elN  = new TGeoElement("Nitrogen", "N", z= 7., a= 14.01);
  TGeoElement* elO  = new TGeoElement("Oxygen"  , "O", z= 8., a= 16.00);
  
  // Materials

  TGeoMixture* matAir
    = new TGeoMixture("Air", 2, density = 1.29e-03);
  matAir->AddElement(elN, 0.7); 
  matAir->AddElement(elO, 0.3); 

  TGeoMaterial* matLead 
    = new TGeoMaterial("Lead", a = 207.19, z = 82., density = 11.35); 


  TGeoMaterial* matXe 
    = new TGeoMaterial("XenonGas", a = 131.29, z = 54., density = 5.458e-03);
    
  Double_t atmosphere = 6.32421e+08;
  Double_t pressure   = 1.*atmosphere;
  Double_t temperature = 293.15;
  matXe->SetPressure(pressure);
  matXe->SetTemperature(temperature);
  matXe->SetState(TGeoMaterial::kMatStateGas);

  // Tracking media

  Double_t param[20];
  param[0] = 0;     // isvol  - Not used
  param[1] = 2;     // ifield - User defined magnetic field
  param[2] = 10.;   // fieldm - Maximum field value (in kiloGauss)
  param[3] = -20.;  // tmaxfd - Maximum angle due to field deflection 
  param[4] = -0.01; // stemax - Maximum displacement for multiple scat 
  param[5] = -.3;   // deemax - Maximum fractional energy loss, DLS 
  param[6] = .001;  // epsil - Tracking precision
  param[7] = -.8;   // stmin
  for ( Int_t i=8; i<20; ++i) param[i] = 0.;

  fImedAir = 1; 
  new TGeoMedium("Air", fImedAir, matAir, param);
  
  fImedPb = 2;
  new TGeoMedium("Lead", fImedPb, matLead, param);
  
  fImedXe = 3;
  new TGeoMedium("XenonGas", fImedXe, matXe, param);
}    

//_____________________________________________________________________________
void Ex02DetectorConstruction::ConstructGeometry()
{
/// Contruct volumes using TGeo modeller

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  Int_t targetMater  = fImedPb;
  Int_t chamberMater = fImedXe;
   
  Double_t targetSize  = 0.5*fTargetLength;    // Half length of the Target  
  Double_t trackerSize = 0.5*fTrackerLength;   // Half length of the Tracker
      
  Double_t* ubuf = 0;

//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 

  Double_t halfWorldLength = 0.5*fWorldLength;
  
  Double_t world[3];
  world[0] = halfWorldLength;
  world[1] = halfWorldLength;
  world[2] = halfWorldLength;
  TGeoVolume *top = gGeoManager->Volume("WRLD","BOX", fImedAir, world, 3);
  gGeoManager->SetTopVolume(top);

  //------------------------------ 
  // Target
  //------------------------------
  
  Double_t target[3];
  target[0] = targetSize;
  target[1] = targetSize;
  target[2] = targetSize;
  gGeoManager->Volume("TARG", "BOX", targetMater, target, 3);

  Double_t posX =  0.;
  Double_t posY =  0.;
  Double_t posZ =  -(targetSize+trackerSize);
  gGeoManager->Node("TARG", 1 ,"WRLD", posX, posY, posZ, 0, kTRUE, ubuf);

  cout << "Target is " << fTargetLength << " cm" 
       //<< " of " << TargetMater->GetName() 
       << " of lead "
       << endl;

  //------------------------------ 
  // Tracker
  //------------------------------
  
  Double_t tracker[3];
  tracker[0] = trackerSize;
  tracker[1] = trackerSize;
  tracker[2] = trackerSize;
  gGeoManager->Volume("TRAK", "BOX", fImedAir, tracker, 3);

  posX = 0.;
  posY = 0.;
  posZ = 0.;
  gGeoManager->Node("TRAK", 1 ,"WRLD", posX, posY, posZ, 0, kTRUE, ubuf);
 
  //------------------------------ 
  // Tracker segments
  //------------------------------
  // 
  // An example of Parameterised volumes
  // dummy values for G4Box -- modified by parameterised volume
  // - implemented using Gsposp

  Double_t chamber[3];
  chamber[0] = -1;
  chamber[1] = -1;
  chamber[2] = -1;
  gGeoManager->Volume("CHMB","BOX", chamberMater, chamber, 0);

  Double_t firstPosition = -trackerSize + 0.5*fChamberWidth;
  Double_t firstLength = fTrackerLength/10;
  Double_t lastLength  = fTrackerLength;
  Ex02ChamberParameterisation* chamberParam 
    = new Ex02ChamberParameterisation(  
			   fNofChambers,          // NoChambers 
			   firstPosition,         // Z of center of first 
			   fChamberSpacing,        // Z spacing of centers
			   fChamberWidth,          // Width Chamber 
			   firstLength,           // lengthInitial 
			   lastLength);           // lengthFinal
  for (Int_t i=0; i<fNofChambers; i++) {
    Double_t pos[3];
    Double_t dim[3];
    chamberParam->ComputeTransformation(i, pos);
    chamberParam->ComputeDimensions(i, dim);
    gGeoManager->Node("CHMB", i ,"TRAK", pos[0], pos[1], pos[2], 0, 
                      kTRUE, dim, 3); 
  }  

  cout << "There are " << fNofChambers 
       << " chambers in the tracker region. "
       << "The chambers are " << fChamberWidth*10. << " mm of " 
       //<< ChamberMater->GetName() 
       << "Xenon gas"
       << "\n The distance between chamber is "
       << fChamberSpacing << " cm" 
       << endl;
	 

  // close geometry
  gGeoManager->CloseGeometry();
    
  // notify VMC about Root geometry
  gMC->SetRootGeometry();

  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 
/*

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerChamberSDname = "ExN02/TrackerChamberSD";
  ExN02TrackerSD* aTrackerSD = new ExN02TrackerSD( trackerChamberSDname );
  SDman->AddNewDetector( aTrackerSD );
  logicChamber->SetSensitiveDetector( aTrackerSD );
*/  

//--------- Visualization attributes -------------------------------
/*
  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  logicTarget ->SetVisAttributes(BoxVisAtt);
  logicTracker->SetVisAttributes(BoxVisAtt);
  
  G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  logicChamber->SetVisAttributes(ChamberVisAtt);
*/  
//--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in N02PhysicsList how to setup the process
  // G4UserSpecialCuts).  
  // Sets a max Step length in the tracker region
  // Double_t maxStep = 0.5*fChamberWidth, maxLength = 2*fTrackerLength;
  // Double_t maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
}

//_____________________________________________________________________________
void Ex02DetectorConstruction::SetTargetMaterial(const TString& /*materialName*/)
{
/// Set target material (not available with VMC)

  Warning("SetTargetMaterial", "Not available in virtual Monte Carlo");
}
 
//_____________________________________________________________________________
void Ex02DetectorConstruction::SetChamberMaterial(const TString& /*materialName*/)
{
/// Set chamber material (not available with VMC)

  Warning("SetTargetMaterial", "Not available in virtual Monte Carlo");
}

// -------------------------------------------------------------------------
// -----                PndDpmGeneratorMod source file                  -----
// -------------------------------------------------------------------------


#include <iostream>
#include "TClonesArray.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "PndDpmGeneratorMod.h"
#include "FairPrimaryGenerator.h"

using namespace std;

// -----   Default constructor   ------------------------------------------
PndDpmGeneratorMod::PndDpmGeneratorMod() {
  iEvent     = 0;
  fInputFile = NULL;
  fInputTree = NULL;
}
// ------------------------------------------------------------------------



// -----   Standard constructor   -----------------------------------------
PndDpmGeneratorMod::PndDpmGeneratorMod(const Char_t* fileName) {
  iEvent     = 0;
  fFileName  = fileName;
  fInputFile = new TFile(fFileName);
  fInputTree = (TTree*) fInputFile->Get("data");
  fParticles = new TClonesArray("TParticle",100);
  fInputTree->SetBranchAddress("Particles", &fParticles);
}
// ------------------------------------------------------------------------



// -----   Destructor   ---------------------------------------------------
PndDpmGeneratorMod::~PndDpmGeneratorMod() {
  CloseInput();
}
// ------------------------------------------------------------------------



// -----   Public method ReadEvent   --------------------------------------
Bool_t PndDpmGeneratorMod::ReadEvent(FairPrimaryGenerator* primGen) {

  // Check for input file
  if ( ! fInputFile ) {
    cout << "-E PndDpmGeneratorMod: Input file nor open!" << endl;
    return kFALSE;
  }

  // Check for number of events in input file
  if ( iEvent > fInputTree->GetEntries() ) {
    cout << "-E PndDpmGeneratorMod: No more events in input file!" << endl;
    CloseInput();
    return kFALSE;
  }
  TFile  *g=gFile;
  fInputFile->cd();
  fInputTree->GetEntry(iEvent++);
  g->cd();

   // Get number of particle in TClonesrray
  Int_t nParts = fParticles->GetEntriesFast();

  // Loop over particles in TClonesArray
  for (Int_t iPart=0; iPart < nParts; iPart++) {
    TParticle* part = (TParticle*) fParticles->At(iPart);
    Int_t pdgType = part->GetPdgCode();

    // Check if particle type is known to database
    if ( ! pdgType ) {
      cout << "-W PndDpmGeneratorMod: Unknown type " << part->GetPdgCode()
	   << ", skipping particle." << endl;
      continue;
    }

    Double_t px = part->Px();
    Double_t py = part->Py();
    Double_t pz = part->Pz();
    Double_t e = part->Energy();

    Double_t vx = part->Vx();
    Double_t vy = part->Vy();
    Double_t vz = part->Vz();

    // Give track to PrimaryGenerator
    if (pdgType != 111){ // not pi0
      primGen->AddTrack(pdgType, px, py, pz, vx, vy, vz, -1, true);
    } else {
      primGen->AddTrack(pdgType, px, py, pz, vx, vy, vz, -1, false);
      TLorentzVector pi0, g1, g2;
      pi0.SetPxPyPzE(px,py,pz,e);
      decay_pi0(pi0,g1,g2);
      primGen->AddTrack(pdgType, g1.Vect().Px(), g1.Vect().Py(), g1.Vect().Pz(), vx, vy, vz, iPart, true);
      primGen->AddTrack(pdgType, g2.Vect().Px(), g2.Vect().Py(), g2.Vect().Pz(), vx, vy, vz, iPart, true);
    }
  }        //  Loop over particle in event

  return kTRUE;

}
// ------------------------------------------------------------------------

void PndDpmGeneratorMod::decay_pi0(const TLorentzVector &pi0, TLorentzVector &g1, TLorentzVector &g2) {

  const TVector3 boost_to_lab = pi0.BoostVector();
  const double E = pi0.M()/2., p=E;
  const double cost = (2.0*gRandom->Uniform())-1.0;
  const double pz = p * cost;
  const double pt = sqrt(p*p-pz*pz);
  const double phi = gRandom->Uniform();
  const double py = pt * sin(phi);
  const double px = pt * cos(phi);

  // GO FIGURE -
  // TODO still not quite there but close enough for now
  TVector3 g(px,py,pz);
  TVector3 dir = -pi0.Vect().Unit();
  g.RotateUz(dir);
  g1.SetPxPyPzE(g.X(),g.Y(),g.Z(),E);
  g1.Boost(boost_to_lab);
  g2.SetPxPyPzE(-g.X(),-g.Y(),-g.Z(),E);
  g2.Boost(boost_to_lab);

}

// -----   Private method CloseInput   ------------------------------------
void PndDpmGeneratorMod::CloseInput() {
  if ( fInputFile ) {
    cout << "-I PndDpmGeneratorMod: Closing input file " << fFileName
	 << endl;
    fInputFile->Close();
    delete fInputFile;
  }
  fInputFile = NULL;
}
// ------------------------------------------------------------------------



ClassImp(PndDpmGeneratorMod)

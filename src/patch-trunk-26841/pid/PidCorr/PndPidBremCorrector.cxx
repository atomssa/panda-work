//-----------------------
// This Class's Header --
//-----------------------
#include "PndPidBremCorrector.h"

//---------------
// C++ Headers --
//---------------
#include <vector>
//#include <set>
//#include <map>
#include <iostream>

// Path of file:
//  ----- $pandaroot/pid/PidCorr

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

#include "PndEmcBump.h"
#include "PndEmcCluster.h"

#include "PndPidCandidate.h"

#include "PndPidBremCorrected4Mom.h"

using std::cout;
using std::endl;

PndPidBremCorrector::PndPidBremCorrector():
   fClusterArray(0), fPhiBumpArray(0), fBumpArray(0), fChargedCandidateArray(0), fNeutralCandidateArray(0), fBremCorrected4MomArray(0),fRecMomOfEle(0), fRecThetaOfEle(0), fRecPhiOfEle(0), fCharge(0), fSepPhotonE(0.), fMergPhotonE(0.), fEmcPhiBumpList(), fPersistance(kTRUE)
{

}

PndPidBremCorrector::~PndPidBremCorrector()
{

}

InitStatus PndPidBremCorrector::Init() {

  // Get RootManager
  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- PndPidBremCorrector::Init: "
	 << "RootManager not instantiated!" << endl;
    return kFATAL;
  }


 fClusterArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcCluster"));
  if ( ! fClusterArray ) {
    cout << "-W- PndEmcMakeBump::Init: "
	 << "No PndEmcCluster array!" << endl;
    return kERROR;
  }

 fPhiBumpArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcPhiBump"));
  if ( ! fPhiBumpArray ) {
    cout << "-W- PndEmcMakePhiBump::Init: "
	 << "No PhiBumpArray array!" << endl;
    return kERROR;
  }

 fBumpArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcBump"));
  if ( ! fBumpArray ) {
    cout << "-W- PndEmcMakeBump::Init: "
	 << "No PndEmcBump array!" << endl;
    return kERROR;
  }

 fChargedCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidChargedCand"));
  if ( ! fChargedCandidateArray ) {
    cout << "-W- PndEmcMakeBump::Init: "
	 << "No PidChargedCand array!" << endl;
    return kERROR;
  }


 fNeutralCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidNeutralCand"));
  if ( ! fNeutralCandidateArray ) {
    cout << "-W- PndEmcMakeBump::Init: "
	 << "No PidNeutralCand array!" << endl;
    return kERROR;
  }


 fBremCorrected4MomArray = new TClonesArray("PndPidBremCorrected4Mom");
  ioman->Register("BremCorrected4Mom","Pid",fBremCorrected4MomArray,fPersistance);

}

void PndPidBremCorrector::Exec(Option_t* opt)
{

  // Reset output array
  if ( ! fBremCorrected4MomArray ) Fatal("Exec", "No BremCorrected4Mom Array");
  fBremCorrected4MomArray->Delete();

  int nChargedCand = fChargedCandidateArray->GetEntriesFast();

  for (int iCand = 0; iCand<nChargedCand; ++iCand){

    PndPidCandidate* theChargedCand = (PndPidCandidate*) fChargedCandidateArray->At(iCand);
    fRecMomOfEle = theChargedCand->GetMomentum().Mag();
    fRecThetaOfEle = theChargedCand->GetMomentum().Theta()*TMath::RadToDeg();
    fRecPhiOfEle = theChargedCand->GetMomentum().Phi()*TMath::RadToDeg();
    fCharge = theChargedCand->GetCharge();

    TVector3 mom = theChargedCand->GetMomentum();
    double ene = theChargedCand->GetEnergy();
    TLorentzVector m4 = TLorentzVector(mom,ene);
    //double mass = m4.M();
    double mass = 0.000511; // Electron mass hypothesis (GeV)

    fSepPhotonE = GetSepPhotonE(theChargedCand);
    fMergPhotonE = GetMergPhotonE(theChargedCand);
    double energy_gamma = fSepPhotonE + fMergPhotonE;

    TVector3 momCorr = ((energy_gamma+fRecMomOfEle)/fRecMomOfEle) * mom;
    double eneCorr = TMath::Hypot(mass, momCorr.Mag());

    PndPidBremCorrected4Mom *bremCorr = AddBremCorrected4Mom();
    bremCorr->SetMomentum(momCorr);
    bremCorr->SetEnergy(eneCorr);
    bremCorr->SetPidCandIdx(iCand);

  }

}

PndPidBremCorrected4Mom* PndPidBremCorrector::AddBremCorrected4Mom(){
  TClonesArray& clref = *fBremCorrected4MomArray;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) PndPidBremCorrected4Mom();
}

double PndPidBremCorrector::GetSepPhotonE(PndPidCandidate *ChargedCand){

  Float_t PhotonTotEnergySepWtd = 0;

  const int iTrkEmcIdx = ChargedCand->GetEmcIndex();
  if (iTrkEmcIdx<0) return 0;

  const int nBump = fBumpArray->GetEntriesFast();
  for(Int_t iBump = 0; iBump<nBump; ++iBump)
    {
      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(iBump);
      const Float_t PhotonEnergySep = PhotonBump->GetEnergyCorrected();

      const Int_t iSepClust = PhotonBump->GetClusterIndex();
      if ( iSepClust == iTrkEmcIdx ) continue;

      const Double_t PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
      const Double_t PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();

      const Bool_t fwd = fRecThetaOfEle <= 23.;
      const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
      const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
      const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;

      const Float_t RealDeltaPhi = fCharge<0?PhotonPhiSep-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiSep;
      const Float_t RealDeltaTheta = fCharge<0?PhotonThetaSep-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaSep;

      const Float_t rad_calc = 100.*TMath::Sin(RealDeltaPhi*TMath::DegToRad()/2.)*2.0*Pt/0.3/2.0; // B=2T
      const Float_t zed_calc = rad_calc/TMath::Tan(TMath::DegToRad()*fRecThetaOfEle);

      const Float_t wt = fwd ? 1.0/(1.+TMath::Exp((zed_calc-90.)/25.)) : 1.0/(1.+TMath::Exp((rad_calc-21.)/5.));
      const Float_t ThetaCutUp = 2.;
      const Float_t ThetaCutDown = -2.;
      const Float_t PhiCutUp = fwd ? DeltaPhiForward : DeltaPhiBarrel;
      const Float_t PhiCutDown = -1;

      const Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
      const Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;

      if (PhiCut && ThetaCut) PhotonTotEnergySepWtd += wt*PhotonEnergySep;

    }//loop neutralcand

  if (PhotonTotEnergySepWtd < fRecMomOfEle/100.) PhotonTotEnergySepWtd = 0;

  return PhotonTotEnergySepWtd;

}

double PndPidBremCorrector::GetMergPhotonE(PndPidCandidate *ChargedCand){

  Double_t PhotonTotEnergyMerg = 0.0;

  // no EMcal cluster associated with track ...
  if (ChargedCand->GetEmcIndex() < 0) return 0.0;

  PndEmcBump *EleBump = (PndEmcBump *) fBumpArray->At(ChargedCand->GetEmcIndex());
  Int_t EleRefCluster = EleBump->GetClusterIndex();

  if (EleRefCluster < 0) return 0.0;
  GetEmcPhiBumpList(EleRefCluster);

  Double_t EnergyCut = 0.15/TMath::Sin(fRecThetaOfEle*TMath::DegToRad());
  Double_t EleEnergy = 0;

  int iMax = 0;
  Float_t eMax = -1e9;
  for (Int_t ib = 0; ib < fEmcPhiBumpList.size(); ++ib) {
    if( fEmcPhiBumpList[ib]->energy() > eMax) {
      iMax = ib;
      eMax = fEmcPhiBumpList[ib]->energy();
    }
  }
  Int_t iS = fCharge<0?0:iMax+1;
  Int_t iE = fCharge<0?iMax-1:fEmcPhiBumpList.size()-1;
  for(Int_t r = iS; r<=iE; r++) PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();

  if(PhotonTotEnergyMerg > fRecMomOfEle/100.)  {
    return PhotonTotEnergyMerg;
  } else {
    return 0.0;
  }

}

void PndPidBremCorrector::GetEmcPhiBumpList(int iClust) {
  fEmcPhiBumpList.clear();
  int nPhiBump = fPhiBumpArray->GetEntriesFast();
  for (int ipb=0; ipb<nPhiBump; ++ipb) {
    PndEmcBump *phibump = (PndEmcBump*) fPhiBumpArray->At(ipb);
    if ( phibump->GetClusterIndex() == iClust ) {
      fEmcPhiBumpList.push_back(phibump);
    }
  }
}

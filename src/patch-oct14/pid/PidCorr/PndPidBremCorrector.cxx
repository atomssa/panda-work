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

  int nNeutralCand = fNeutralCandidateArray->GetEntriesFast();

  Float_t PhotonTotEnergySep = 0;
  for(Int_t iNeutralCand = 0; iNeutralCand<nNeutralCand; ++iNeutralCand)
    {
      Float_t PhotonEnergySep = 0;
      PndPidCandidate *PhotonCand = (PndPidCandidate *) fNeutralCandidateArray->At(iNeutralCand);
      PhotonEnergySep = PhotonCand->GetEmcCalEnergy();
      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(PhotonCand->GetEmcIndex());
      double PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
      double PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();

      Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
      Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
      Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;
      Float_t RealDeltaPhi = 0, RealDeltaTheta =0;

      if (fCharge < 0){
	RealDeltaPhi = PhotonPhiSep - fRecPhiOfEle;
	RealDeltaTheta = PhotonThetaSep - fRecThetaOfEle;
      }
      else {
	RealDeltaPhi = fRecPhiOfEle - PhotonPhiSep;
	RealDeltaTheta =fRecThetaOfEle - PhotonThetaSep;
      }

      Float_t PhiCutUp = 0, ThetaCutUp = 0 ,  PhiCutDown = 0, ThetaCutDown = 0;

      if (fRecThetaOfEle <= 23.)
	{
	  PhiCutUp = DeltaPhiForward;
	  PhiCutDown = -1.;
	  ThetaCutUp = 2.;
	  ThetaCutDown = -2.;
	}
      else if (fRecThetaOfEle > 23.)
	{
	  PhiCutUp = DeltaPhiBarrel;
	  PhiCutDown = -1.;
	  ThetaCutUp = 2.;
	  ThetaCutDown = -2.;
	}
      Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
      Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
      if (PhiCut && ThetaCut) PhotonTotEnergySep += PhotonEnergySep;

    }//loop neutralcand

  if(PhotonTotEnergySep > fRecMomOfEle/100.)
    {
      return PhotonTotEnergySep;
    }
  else return 0.;

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

  if (fCharge<0) {
    for (Int_t i_phibump = fEmcPhiBumpList.size()-1; i_phibump >= 0; --i_phibump) {
      if( fEmcPhiBumpList[i_phibump]->energy() > EnergyCut) {
	for(Int_t r = 0; r<i_phibump; r++) PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();
	i_phibump = -1;
      }
    }
  } else {
    for (Int_t i_phibump =0; i_phibump < fEmcPhiBumpList.size(); ++i_phibump) {
      if( fEmcPhiBumpList[i_phibump]->energy() > EnergyCut) {
	for(Int_t r = 0; r<i_phibump; r++) PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();
	i_phibump = fEmcPhiBumpList.size()+1;
      }
    }
  }

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

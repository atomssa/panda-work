//-----------------------
// This Class's Header --
//-----------------------
#include "PndPidBremCorrectorNT.h"

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
#include "FairMCEventHeader.h"

#include "PndEmcPoint.h"
#include "PndEmcHit.h"
#include "PndEmcDigi.h"
#include "PndEmcBump.h"
#include "PndEmcCluster.h"

#include "PndPidCandidate.h"

#include "PndMCTrack.h"

#include "PndPidBremCorrected4Mom.h"

using std::cout;
using std::endl;

PndPidBremCorrectorNT::PndPidBremCorrectorNT():
   fClusterArray(0), fPhiBumpArray(0), fBumpArray(0), fChargedCandidateArray(0), fNeutralCandidateArray(0), fBremCorrected4MomArray(0),fRecMomOfEle(0), fRecThetaOfEle(0), fRecPhiOfEle(0), fCharge(0), fSepPhotonE(0.), fMergPhotonE(0.), fEmcPhiBumpList(), fPersistance(kTRUE)
{
}

PndPidBremCorrectorNT::~PndPidBremCorrectorNT()
{

}

InitStatus PndPidBremCorrectorNT::Init() {

  cout <<"PndPidBremCorrectorNT::Init " << endl;

  // Get RootManager
  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- PndPidBremCorrectorNT::Init: "
	 << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  fMCHeader = dynamic_cast<FairMCEventHeader*>(ioman->GetObjectFromInTree("MCEventHeader."));
  if ( ! fMCHeader ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No fMCHeader array!" << endl;
    return kERROR;
  }
  ioman->Register("MCEventHeader.","Event",fMCHeader,kTRUE);

  fClusterArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcCluster"));
  if ( ! fClusterArray ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No PndEmcCluster array!" << endl;
    return kERROR;
  }
  ioman->Register("EmcCluster","Emc",fClusterArray,kTRUE);

  fPhiBumpArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcPhiBump"));
  if ( ! fPhiBumpArray ) {
    cout << "-W- PndEmcMakePhiBump::Init: "
	 << "No PhiBumpArray array!" << endl;
    return kERROR;
  }
  ioman->Register("EmcPhiBump","Emc",fPhiBumpArray,kTRUE);

  fBumpArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcBump"));
  if ( ! fBumpArray ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No PndEmcBump array!" << endl;
    return kERROR;
  }
  ioman->Register("EmcBump","Emc",fBumpArray,kTRUE);

  fDigiArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcDigi"));
  if ( ! fDigiArray ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No PndEmcDigi array!" << endl;
    return kERROR;
  }
  ioman->Register("EmcDigi","Emc",fDigiArray,kTRUE);

  fHitArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcHit"));
  if ( ! fHitArray ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No PndEmcHit array!" << endl;
    return kERROR;
  }
  ioman->Register("EmcHit","Emc",fHitArray,kTRUE);

  fChargedCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidChargedCand"));
  if ( ! fChargedCandidateArray ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No PidChargedCand array!" << endl;
    return kERROR;
  }

  fNeutralCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidNeutralCand"));
  if ( ! fNeutralCandidateArray ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No PidNeutralCand array!" << endl;
    return kERROR;
  }

  fMcArray = dynamic_cast<TClonesArray *> (ioman->GetObject("MCTrack"));
  if ( ! fMcArray ) {
    cout << "-W- PndPidBremCorrectorNT::Init: "
	 << "No McTrack array!" << endl;
    return kERROR;
  }
  ioman->Register("MCTrack","Stack",fMcArray,kTRUE);

  cout <<"PndPidBremCorrectorNT::Init Creating PndPidBremCorrected4Mom" << endl;

  fBremCorrected4MomArray = new TClonesArray("PndPidBremCorrected4Mom");
  ioman->Register("BremCorrected4Mom","Pid",fBremCorrected4MomArray,fPersistance);

}

void PndPidBremCorrectorNT::print_cands() {

  int clustId = -1;
  cout << "============================ printing event ===============================" <<endl;
  nChCand = fChargedCandidateArray->GetEntriesFast();
  cout << "------ nChCand= " << nChCand << " -------" << endl;
  for (int iChCand = 0; iChCand<nChCand; ++iChCand){
    PndPidCandidate* theChargedCand = (PndPidCandidate*) fChargedCandidateArray->At(iChCand);
    cout << "iChCand= " << iChCand;
    cout << " (" << (theChargedCand->GetCharge()>0?"+":"-") << ")";
    cout << " P= " << theChargedCand->GetMomentum().Mag();
    cout << " Ang= (" << theChargedCand->GetMomentum().Theta()*TMath::RadToDeg();
    cout << ", " << theChargedCand->GetMomentum().Phi()*TMath::RadToDeg() << ")";
    cout << " McIdx= " << theChargedCand->GetMcIndex();
    cout << " clust= " << theChargedCand->GetEmcIndex();
    clustId = theChargedCand->GetEmcIndex();
    cout << endl;
  }

  nNeutCand = fNeutralCandidateArray->GetEntriesFast();
  cout << "------ nNeutCand= " << nNeutCand << " -------" << endl;
  for (int iNeutCand = 0; iNeutCand<nNeutCand; ++iNeutCand){
    PndPidCandidate *PhotonCand = (PndPidCandidate *) fNeutralCandidateArray->At(iNeutCand);
    cout << "iNeutCand= " << iNeutCand;
    cout << " Mom= " << PhotonCand->GetMomentum().Mag();
    cout << " The= " << PhotonCand->GetMomentum().Theta()*TMath::RadToDeg();
    cout << " Phi= " << PhotonCand->GetMomentum().Phi()*TMath::RadToDeg();
    cout << " Ch= " << PhotonCand->GetCharge();
    cout << " McIdx= " << PhotonCand->GetMcIndex();
    cout << " EmcIdx= " << PhotonCand->GetEmcIndex();
    cout << endl;
  }

  int nClust = fClusterArray->GetEntriesFast();
  cout << "------ nClust= " << nClust << " -------" << endl;
  for (int ii = 0; ii < nClust; ++ii) {
    PndEmcCluster *clust = (PndEmcCluster*) fClusterArray->At(ii);
    if (clustId==ii) cout << "--> ";

    for (int jj=0; jj<fSepClustId.size(); ++jj) {
      if (ii==jj) cout << "*** ";
    }

    cout << "ClustId= " << ii;

    //int mcsize = clust->GetMcSize();
    //cout << " McSize= " << mcsize << " McList: (";
    //for (int i=0; i<mcsize; ++i) {
    //  cout << clust->GetMcIndex(i) << (i<(mcsize-1)?", ": "");
    //}
    //cout << ")";
    cout << " Ang=(" << clust->position().Theta()*TMath::RadToDeg();
    cout << ", " << clust->position().Phi()*TMath::RadToDeg() << ") ";

    cout << " E= " << clust->energy();

    int digisize = clust->DigiList().size();
    cout << " Digi= " << digisize << "(";
    for (int i=0; i<digisize; ++i) {
      //cout << ((PndEmcDigi*)fDigiArray->At(clust->DigiList()[i]))->GetTrackId() << (i<(digisize-1)?", ": "");
      cout << clust->DigiList()[i] << (i<(digisize-1)?", ": "");
    }
    cout << ")";
    cout << endl;
  }

  int nBump = fBumpArray->GetEntriesFast();
  cout << "------ nBump= " << nBump << " -------" << endl;
  for (int ii = 0; ii < nBump; ++ii) {
    PndEmcBump *bump = (PndEmcBump*) fBumpArray->At(ii);
    cout << "BumpId= " << ii;
    cout << " clust= " << bump->GetClusterIndex();
    cout << " Ang=(" << bump->position().Theta()*TMath::RadToDeg();
    cout << ", " << bump->position().Phi()*TMath::RadToDeg() << ") ";
    cout << " E= " << bump->energy();

    //int mcsize = bump->GetMcSize();
    //cout << " McSize= " << mcsize << " McList: (";
    //for (int i=0; i<mcsize; ++i) {
    //  cout << bump->GetMcIndex(i) << (i<(mcsize-1)?", ": "");
    //}
    //cout << ")";
    cout << endl;
  }

  int nDigi = fDigiArray->GetEntriesFast();
  cout << "------ nDigi= " << nDigi << " -------" << endl;
  for (int ii = 0; ii < nDigi; ++ii) {
    PndEmcDigi *digi = (PndEmcDigi*) fDigiArray->At(ii);
    cout << "DigiId= " << ii << " trkId= " << digi->GetTrackId() << " hitId= " << digi->GetHitIndex();
    PndEmcHit *hit = (PndEmcHit*) fHitArray->At(digi->GetHitIndex());

    int mcsize = hit->GetMcList().size();
    cout << " McSize= " << mcsize << " McList: (";
    for (int i=0; i<mcsize; ++i) {
      cout << hit->GetMcList()[i] << (i<(mcsize-1)?", ": "");
      if (i==10)  {
	cout << "... ";
	break;
      }

    }
    cout << ")";

    int ptsize = hit->GetPointList().size();
    cout << " PointSize= " << ptsize << " PointList: (";
    for (int i=0; i<ptsize; ++i) {
      cout << hit->GetPointList()[i] << (i<(ptsize-1)?", ": "");
    }
    cout << ")";

    cout << endl;
  }

  int nMcTrack = fMcArray->GetEntriesFast();
  cout << "------ nMcTrack= " << nMcTrack << " ------" << endl;
  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack){
    PndMCTrack *McTrack = (PndMCTrack *) fMcArray->At(iMcTrack);

    double mom = McTrack->GetMomentum().Mag();
    int motherid = McTrack->GetMotherID();
    int motherid2 = McTrack->GetSecondMotherID();
    double mother_mom = 0.;
    if (motherid==0) {
      PndMCTrack *McMother = (PndMCTrack*)fMcArray->At(motherid);
      mother_mom = McMother->GetMomentum().Mag();
    }

    double radiusT = McTrack->GetStartVertex().Pt();
    //if (motherid==-1||(motherid==0&&mom>0.1*mother_mom)) {
    //if (motherid==-1||(radiusT<50/*cm*/)) {
    //if (motherid==-1||(motherid==0&&mom>0.01*motherid&&radiusT<150)) {
    //if (motherid==-1||motherid==0) {
    if (radiusT<4) {
      cout << "i= " << iMcTrack;
      cout << " (" << McTrack->GetPdgCode();
      cout << " M= "<< McTrack->GetMotherID() << ")";
      //cout << " 2M= "<< McTrack->GetSecondMotherID() << ")";
      cout << " P=(" << McTrack->GetMomentum().Mag();
      cout << ", " << McTrack->GetMomentum().Theta()*TMath::RadToDeg();
      cout << "," << McTrack->GetMomentum().Phi()*TMath::RadToDeg() << ")";
      //cout << " V=(" << McTrack->GetStartVertex().X();
      //cout << ", " << McTrack->GetStartVertex().Y();
      //cout << ", " << McTrack->GetStartVertex().Z();
      cout << " R=" << radiusT;
      //cout << ") ";
      cout << endl;
    }
  }

}

void PndPidBremCorrectorNT::Exec(Option_t* opt)
{

  // Reset output array
  if ( ! fBremCorrected4MomArray ) Fatal("Exec", "No BremCorrected4Mom Array");
  fBremCorrected4MomArray->Delete();

  nChCand = fChargedCandidateArray->GetEntriesFast();
  nNeutCand = fNeutralCandidateArray->GetEntriesFast();

  for (int iCand = 0; iCand<nChCand; ++iCand){

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

    int nPhotSep = 0;
    fSepPhotonE = GetSepPhotonE(theChargedCand, nPhotSep);

    int nPhotMrg = 0;
    fMergPhotonE = GetMergPhotonE(theChargedCand, nPhotMrg);

    double energy_gamma = fSepPhotonE + fMergPhotonE;

    TVector3 momCorr = ((energy_gamma+fRecMomOfEle)/fRecMomOfEle) * mom;
    double eneCorr = TMath::Hypot(mass, momCorr.Mag());

    PndPidBremCorrected4Mom *bremCorr = AddBremCorrected4Mom();
    bremCorr->SetMomentum(momCorr);
    bremCorr->SetEnergy(eneCorr);
    bremCorr->SetPidCandIdx(iCand);

  }

  //print_cands();
}



PndPidBremCorrected4Mom* PndPidBremCorrectorNT::AddBremCorrected4Mom(){
  TClonesArray& clref = *fBremCorrected4MomArray;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) PndPidBremCorrected4Mom();
}



double PndPidBremCorrectorNT::GetSepPhotonE(PndPidCandidate *ChargedCand, int &nphotsep){

  fSepClustId.clear();

  Float_t PhotonTotEnergySep = 0;
  for(Int_t iNeutralCand = 0; iNeutralCand<nNeutCand; ++iNeutralCand)
    {
      Float_t PhotonEnergySep = 0;
      PndPidCandidate *PhotonCand = (PndPidCandidate *) fNeutralCandidateArray->At(iNeutralCand);
      PhotonEnergySep = PhotonCand->GetEmcCalEnergy();
      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(PhotonCand->GetEmcIndex());

      Int_t iSepClust = PhotonCand->GetEmcIndex();

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
      if (PhiCut && ThetaCut) {
	fSepClustId.push_back(iSepClust);
	PhotonTotEnergySep += PhotonEnergySep;
	nphotsep++;
      }

    }//loop neutralcand

  if(PhotonTotEnergySep > fRecMomOfEle/100.)
    {
      return PhotonTotEnergySep;
    }
  else return 0.;

}

double PndPidBremCorrectorNT::GetMergPhotonE(PndPidCandidate *ChargedCand, int &nphotmrg){

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
	for(Int_t r = 0; r<i_phibump; r++) {
	  PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();
	  nphotmrg++;
	}
	i_phibump = -1;
      }
    }
  } else {
    for (Int_t i_phibump =0; i_phibump < fEmcPhiBumpList.size(); ++i_phibump) {
      if( fEmcPhiBumpList[i_phibump]->energy() > EnergyCut) {
	for(Int_t r = 0; r<i_phibump; r++) {
	  PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();
	  nphotmrg++;
	}
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

void PndPidBremCorrectorNT::GetEmcPhiBumpList(int iClust) {
  fEmcPhiBumpList.clear();
  int nPhiBump = fPhiBumpArray->GetEntriesFast();
  for (int ipb=0; ipb<nPhiBump; ++ipb) {
    PndEmcBump *phibump = (PndEmcBump*) fPhiBumpArray->At(ipb);
    if ( phibump->GetClusterIndex() == iClust ) {
      fEmcPhiBumpList.push_back(phibump);
    }
  }
}

void PndPidBremCorrectorNT::FinishTask() {
}

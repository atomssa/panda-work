//-----------------------
// This Class's Header --
//-----------------------
#include "BremPidReader.h"
#include "TFile.h"
#include "TTree.h"

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

BremPidReader::BremPidReader():
   fClusterArray(0), fPhiBumpArray(0), fBumpArray(0), fChargedCandidateArray(0), fNeutralCandidateArray(0), fBremCorrected4MomArray(0),fRecMomOfEle(0), fRecThetaOfEle(0), fRecPhiOfEle(0), fCharge(0), fSepPhotonE(0.), fMergPhotonE(0.), fEmcPhiBumpList(), fPersistance(kTRUE)
{
  output_name = "bremcorr.root";
  nEvt = 0;
}

BremPidReader::~BremPidReader()
{

}

InitStatus BremPidReader::Init() {

  cout <<"BremPidReader::Init " << endl;

  // Get RootManager
  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- BremPidReader::Init: "
	 << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  fMCHeader = dynamic_cast<FairMCEventHeader*>(ioman->GetObjectFromInTree("MCEventHeader."));
  if ( ! fMCHeader ) {
    cout << "-W- BremPidReader::Init: "
	 << "No fMCHeader array!" << endl;
    return kERROR;
  }

  fClusterArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcCluster"));
  if ( ! fClusterArray ) {
    cout << "-W- BremPidReader::Init: "
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
    cout << "-W- BremPidReader::Init: "
	 << "No PndEmcBump array!" << endl;
    return kERROR;
  }

  fDigiArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcDigi"));
  if ( ! fDigiArray ) {
    cout << "-W- BremPidReader::Init: "
	 << "No PndEmcDigi array!" << endl;
    return kERROR;
  }

  fHitArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcHit"));
  if ( ! fHitArray ) {
    cout << "-W- BremPidReader::Init: "
	 << "No PndEmcHit array!" << endl;
    return kERROR;
  }

  fChargedCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidChargedCand"));
  if ( ! fChargedCandidateArray ) {
    cout << "-W- BremPidReader::Init: "
	 << "No PidChargedCand array!" << endl;
    return kERROR;
  }

  fNeutralCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidNeutralCand"));
  if ( ! fNeutralCandidateArray ) {
    cout << "-W- BremPidReader::Init: "
	 << "No PidNeutralCand array!" << endl;
    return kERROR;
  }

  fMcArray = dynamic_cast<TClonesArray *> (ioman->GetObject("MCTrack"));
  if ( ! fMcArray ) {
    cout << "-W- BremPidReader::Init: "
	 << "No McTrack array!" << endl;
    return kERROR;
  }

  cout <<"BremPidReader::Init Creating PndPidBremCorrected4Mom" << endl;

  fBremCorrected4MomArray = dynamic_cast<TClonesArray *> (ioman->GetObject("BremCorrected4Mom"));
  if ( ! fBremCorrected4MomArray ) {
    cout << "-W- BremPidReader::Init: "
	 << "No BremCorrected4Mom array!" << endl;
    return kERROR;
  }

  f = new TFile(output_name,"RECREATE");
  t = new TTree("t","t");
  t->Branch("nch",&nChCand,"nch/I");
  t->Branch("neut",&nNeutCand,"neut/I");

  t->Branch("ch_charge",&ch_charge,"ch_charge[nch]/I");
  t->Branch("ch_mom_mc",&ch_mom_mc,"ch_mom_mc[nch]/F");
  t->Branch("ch_mom_rec",&ch_mom_rec,"ch_mom_rec[nch]/F");
  t->Branch("ch_mom_cor",&ch_mom_cor,"ch_mom_cor[nch]/F");
  t->Branch("ch_mom_mrg",&ch_mom_mrg,"ch_mom_mrg[nch]/F");
  t->Branch("ch_mom_sep",&ch_mom_sep,"ch_mom_sep[nch]/F");
  t->Branch("ch_mom_stored",&ch_mom_stored,"ch_mom_stored[nch]/F");
  t->Branch("ch_phi",&ch_phi,"ch_phi[nch]/F");
  t->Branch("ch_the",&ch_the,"ch_the[nch]/F");
  t->Branch("ch_nphot_sep",&ch_nphot_sep,"ch_nphot_sep[nch]/I");
  t->Branch("ch_nphot_mrg",&ch_nphot_mrg,"ch_nphot_mrg[nch]/I");
  t->Branch("ch_is_prim",&ch_is_prim,"ch_is_prim[nch]/I");

  t->Branch("brem_rad",&brem_rad,"brem_rad[nch]/F");
  t->Branch("brem_e",&brem_rad,"brem_e[nch]/F");
  t->Branch("nbrem_trk",&nbrem_trk,"nbrem_trk[nch]/I");

}

void BremPidReader::print_cands() {

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

void BremPidReader::Exec(Option_t* opt)
{

  // Reset output array
  //if ( ! fBremCorrected4MomArray ) Fatal("Exec", "No BremCorrected4Mom Array");
  //fBremCorrected4MomArray->Delete();
  int verb = 0;

  if (verb>0||(nEvt%1000==0))
    cout << "==== New Event " << nEvt << " ===== " << endl;

  nChCand = fChargedCandidateArray->GetEntriesFast();
  nNeutCand = fNeutralCandidateArray->GetEntriesFast();

  for (int iCand = 0; iCand<nChCand; ++iCand){

    if (verb>1)
      cout << "--------------- New Reco Track " << iCand << "------------------" << endl;

    PndPidCandidate* theChargedCand = (PndPidCandidate*) fChargedCandidateArray->At(iCand);
    fRecMomOfEle = theChargedCand->GetMomentum().Mag();
    fRecThetaOfEle = theChargedCand->GetMomentum().Theta()*TMath::RadToDeg();
    fRecPhiOfEle = theChargedCand->GetMomentum().Phi()*TMath::RadToDeg();
    fCharge = theChargedCand->GetCharge();

    ch_charge[iCand] = fCharge;
    ch_mom_rec[iCand] = fRecMomOfEle;
    ch_phi[iCand] = fRecPhiOfEle;
    ch_the[iCand] = fRecThetaOfEle;

    // get the mctruth
    Int_t mcidx = theChargedCand->GetMcIndex();
    // The following line (condition on mcidx) should not be allowed in production code!!
    // it will lead to mismatch between PndPidCand idx and PndPidBremCorr Index
    if (mcidx>fMcArray->GetEntriesFast() || mcidx<0) continue;
    if (verb>1)
      cout << "mcidx= " << mcidx << endl;
    ch_is_prim[iCand] = (mcidx<=4);

    PndMCTrack *truth = (PndMCTrack*) fMcArray->At(mcidx);
    ch_mom_mc[iCand] = truth->GetMomentum().Mag();

    TVector3 mom = theChargedCand->GetMomentum();
    double ene = theChargedCand->GetEnergy();
    TLorentzVector m4 = TLorentzVector(mom,ene);
    //double mass = m4.M();
    double mass = 0.000511; // Electron mass hypothesis (GeV)

    int nPhotSep = 0;
    //fSepPhotonE = GetSepPhotonE(theChargedCand, nPhotSep);
    fSepPhotonE = GetSepPhotonE_fromBumps(theChargedCand, nPhotSep);
    ch_nphot_sep[iCand] = nPhotSep;
    TVector3 momCorrSep = ((fSepPhotonE+fRecMomOfEle)/fRecMomOfEle) * mom;
    ch_mom_sep[iCand] = momCorrSep.Mag();

    int nPhotMrg = 0;
    fMergPhotonE = GetMergPhotonE(theChargedCand, nPhotMrg);
    ch_nphot_mrg[iCand] = nPhotMrg;
    TVector3 momCorrMrg = ((fMergPhotonE+fRecMomOfEle)/fRecMomOfEle) * mom;
    ch_mom_mrg[iCand] = momCorrMrg.Mag();

    double energy_gamma = fSepPhotonE + fMergPhotonE;
    TVector3 momCorr = ((energy_gamma+fRecMomOfEle)/fRecMomOfEle) * mom;
    double eneCorr = TMath::Hypot(mass, momCorr.Mag());
    ch_mom_cor[iCand] = momCorr.Mag();

    PndPidBremCorrected4Mom *tmp = (PndPidBremCorrected4Mom*)fBremCorrected4MomArray->At(iCand);
    if(tmp) {
      TVector3 storedCorrMom = tmp->GetMomentum();
      ch_mom_stored[iCand] = storedCorrMom.Mag();
    } else {
      ch_mom_stored[iCand] = -9999.;
    }

    if (verb>1) {
      cout << "momRec= " << mom.Mag() << " ";
      mom.Print();
      cout << "momCorr= " << momCorr.Mag() << " ";
      momCorr.Print();
      cout << "storedCorrMom= ";
      if (tmp) {
	TVector3 storedCorrMom = tmp->GetMomentum();
	cout << storedCorrMom.Mag() << " ";
	storedCorrMom.Print();
      } else {
	cout << "not found... :("  << endl;
      }
    }

    Double_t _radT=0, _e=0;
    int nBremTrk = 0;
    get_earliest_brem(nBremTrk,_radT,_e);
    brem_rad[iCand] = _radT;
    brem_e[iCand] = _e;
    nbrem_trk[iCand] = nBremTrk;

  }

  //print_cands();

  t->Fill();

  nEvt++;

}

void BremPidReader::get_earliest_brem(int &n, double &r, double &e) {
  n = 0;
  r = 1e9;
  int nMcTrack = fMcArray->GetEntriesFast();
  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack){
    PndMCTrack *McTrack = (PndMCTrack *) fMcArray->At(iMcTrack);
    if (McTrack->GetPdgCode()!=22) continue; // only interested in photons
    if (McTrack->GetMotherID()!=0) continue; // only photons emerging from primary electron
    double radiusT = McTrack->GetStartVertex().Pt();
    if (radiusT<44) n++; // counter for the number of photons emmitted directly by electron in the tracking volume
    if (radiusT<r) {
      r = radiusT;
      e = McTrack->GetMomentum().Mag();
    }
  }
}

PndPidBremCorrected4Mom* BremPidReader::AddBremCorrected4Mom(){
  TClonesArray& clref = *fBremCorrected4MomArray;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) PndPidBremCorrected4Mom();
}

double BremPidReader::GetSepPhotonE(PndPidCandidate *ChargedCand, int &nphotsep){

  fSepClustId.clear();

  Float_t PhotonTotEnergySep = 0;
  for(Int_t iNeutralCand = 0; iNeutralCand<nNeutCand; ++iNeutralCand)
    {
      Float_t PhotonEnergySep = 0;
      PndPidCandidate *PhotonCand = (PndPidCandidate *) fNeutralCandidateArray->At(iNeutralCand);
      PhotonEnergySep = PhotonCand->GetEmcCalEnergy();
      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(PhotonCand->GetEmcIndex());

      Int_t iSepClust = PhotonCand->GetEmcIndex();

      if ( PhotonCand->GetEmcIndex() == ChargedCand->GetEmcIndex() ) continue;
      //if ( PhotonEnergySep > 0.8* ChargedCand->GetEnergy() ) continue;

      const double PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
      const double PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();

      const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
      const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
      const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;
      //cout << "DphBar= " << DeltaPhiBarrel << " DphFor= " << DeltaPhiForward << endl;

      const Float_t RealDeltaPhi = fCharge<0?PhotonPhiSep-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiSep;
      const Float_t RealDeltaTheta = fCharge<0?PhotonThetaSep-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaSep;

      const Float_t ThetaCutUp = 2.;
      const Float_t ThetaCutDown = -2.;
      const Float_t PhiCutUp = (fRecThetaOfEle <= 23.)?DeltaPhiForward:DeltaPhiBarrel;
      const Float_t PhiCutDown = -1;

      const Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
      const Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
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

//double BremPidReader::GetSepPhotonE(PndPidCandidate *ChargedCand, int &nphotsep){
//
//  fSepClustId.clear();
//
//  Float_t PhotonTotEnergySep = 0;
//  for(Int_t iNeutralCand = 0; iNeutralCand<nNeutCand; ++iNeutralCand)
//    {
//      Float_t PhotonEnergySep = 0;
//      PndPidCandidate *PhotonCand = (PndPidCandidate *) fNeutralCandidateArray->At(iNeutralCand);
//      PhotonEnergySep = PhotonCand->GetEmcCalEnergy();
//      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(PhotonCand->GetEmcIndex());
//
//      Int_t iSepClust = PhotonCand->GetEmcIndex();
//
//      if ( PhotonCand->GetEmcIndex() == ChargedCand->GetEmcIndex() ) continue;
//      //if ( PhotonEnergySep > 0.8* ChargedCand->GetEnergy() ) continue;
//
//      double PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
//      double PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();
//
//      Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
//      Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
//      Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;
//      Float_t RealDeltaPhi = 0, RealDeltaTheta =0;
//
//      if (fCharge < 0){
//	RealDeltaPhi = PhotonPhiSep - fRecPhiOfEle;
//	RealDeltaTheta = PhotonThetaSep - fRecThetaOfEle;
//      } else {
//	RealDeltaPhi = fRecPhiOfEle - PhotonPhiSep;
//	RealDeltaTheta =fRecThetaOfEle - PhotonThetaSep;
//      }
//
//      Float_t ThetaCutUp = 4.;
//      Float_t ThetaCutDown = -4.;
//      Float_t PhiCutUp = (fRecThetaOfEle <= 23.)?DeltaPhiForward:DeltaPhiBarrel;
//      Float_t PhiCutDown = -1;
//
//      Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
//      Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
//      if (PhiCut && ThetaCut) {
//	fSepClustId.push_back(iSepClust);
//	PhotonTotEnergySep += PhotonEnergySep;
//	nphotsep++;
//      }
//
//    }//loop neutralcand
//
//  if(PhotonTotEnergySep > fRecMomOfEle/100.)
//    {
//      return PhotonTotEnergySep;
//    }
//  else return 0.;
//
//}

double BremPidReader::GetSepPhotonE_fromBumps(PndPidCandidate *ChargedCand, int &nphotsep){

  fSepClustId.clear();

  Float_t PhotonTotEnergySep = 0;
  const int nBump = fBumpArray->GetEntriesFast();
  for(Int_t iBump = 0; iBump<nBump; ++iBump)
    {
      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(iBump);
      const Float_t PhotonEnergySep = PhotonBump->GetEnergyCorrected();

      const Int_t iSepClust = PhotonBump->GetClusterIndex();
      if ( PhotonBump->GetClusterIndex() == ChargedCand->GetEmcIndex() ) continue;
      //if ( PhotonEnergySep > 0.8* ChargedCand->GetEnergy() ) continue;

      const double PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
      const double PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();

      const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
      const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
      const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;
      //cout << "DphBar= " << DeltaPhiBarrel << " DphFor= " << DeltaPhiForward << endl;

      const Float_t RealDeltaPhi = fCharge<0?PhotonPhiSep-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiSep;
      const Float_t RealDeltaTheta = fCharge<0?PhotonThetaSep-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaSep;

      const Float_t ThetaCutUp = 2.;
      const Float_t ThetaCutDown = -2.;
      const Float_t PhiCutUp = (fRecThetaOfEle <= 23.)?DeltaPhiForward:DeltaPhiBarrel;
      const Float_t PhiCutDown = -1;

      const Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
      const Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
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


double BremPidReader::GetMergPhotonE(PndPidCandidate *ChargedCand, int &nphotmrg){

  Double_t PhotonTotEnergyMerg = 0.0;

  // no EMcal cluster associated with track ...
  if (ChargedCand->GetEmcIndex() < 0) return 0.0;

  PndEmcBump *EleBump = (PndEmcBump *) fBumpArray->At(ChargedCand->GetEmcIndex());
  Int_t EleRefCluster = EleBump->GetClusterIndex();

  if (EleRefCluster < 0) return 0.0;
  GetEmcPhiBumpList(EleRefCluster);

  // This WORKS
  int iMax = 0;
  Float_t eMax = -1e9;
  //for (Int_t ib = fEmcPhiBumpList.size()-1; ib >= 0; --ib) {
  for (Int_t ib = 0; ib < fEmcPhiBumpList.size(); ++ib) {
    if( fEmcPhiBumpList[ib]->energy() > eMax) {
      iMax = ib;
      eMax = fEmcPhiBumpList[ib]->energy();
    }
  }
  Int_t iS = fCharge<0?0:iMax+1;
  Int_t iE = fCharge<0?iMax-1:fEmcPhiBumpList.size()-1;
  for(Int_t r = iS; r<=iE; r++) {
    PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();
    nphotmrg++;
  }

  //Double_t EnergyCut = 0.15/TMath::Sin(fRecThetaOfEle*TMath::DegToRad());
  //if (fCharge<0) {
  //  for (Int_t i_phibump = fEmcPhiBumpList.size()-1; i_phibump >= 0; --i_phibump) {
  //    if( fEmcPhiBumpList[i_phibump]->energy() > EnergyCut) {
  //	for(Int_t r = 0; r<i_phibump; r++) {
  //	  PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();
  //	  nphotmrg++;
  //	}
  //	i_phibump = -1;
  //    }
  //  }
  //} else {
  //  for (Int_t i_phibump =0; i_phibump < fEmcPhiBumpList.size(); ++i_phibump) {
  //    if( fEmcPhiBumpList[i_phibump]->energy() > EnergyCut) {
  //	for(Int_t r = i_phibump+1; r<fEmcPhiBumpList.size(); r++) {
  //	  PhotonTotEnergyMerg += fEmcPhiBumpList[r]->energy();
  //	  nphotmrg++;
  //	}
  //	i_phibump = fEmcPhiBumpList.size()+1;
  //    }
  //  }
  //}

  if(PhotonTotEnergyMerg > fRecMomOfEle/100.)  {
    return PhotonTotEnergyMerg;
  } else {
    return 0.0;
  }

}

void BremPidReader::GetEmcPhiBumpList(int iClust) {
  fEmcPhiBumpList.clear();
  int nPhiBump = fPhiBumpArray->GetEntriesFast();
  for (int ipb=0; ipb<nPhiBump; ++ipb) {
    PndEmcBump *phibump = (PndEmcBump*) fPhiBumpArray->At(ipb);
    if ( phibump->GetClusterIndex() == iClust ) {
      fEmcPhiBumpList.push_back(phibump);
    }
  }
}

void BremPidReader::FinishTask() {
  f->cd();
  t->Write();
}

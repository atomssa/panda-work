//-----------------------
// This Class's Header --
//-----------------------
#include "McBumpMatch.h"
#include "TFile.h"
#include "TTree.h"

//---------------
// C++ Headers --
//---------------
#include <vector>
#include <iostream>

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

#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>

using namespace std;

McBumpMatch::McBumpMatch():
  nEvt(0),
  nMcTrack(0),
  nChCand(0),
  nNeutCand(0),
  radMaxTracking_cm(0.0),
  fBumpArray(0),
  fClusterArray(0),
  fDigiArray(0),
  fHitArray(0),
  fChargedCandidateArray(0),
  fNeutralCandidateArray(0),
  fMcArray(0)
{
  nEvt = 0;
  radMaxTracking_cm = 42.0;
}

InitStatus McBumpMatch::Init() {

  cout <<"McBumpMatch::Init " << endl;

  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- McBumpMatch::Init: "
	 << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  fBumpArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcBump"));
  if ( ! fBumpArray ) {
    cout << "-W- McBumpMatch::Init: "
	 << "No PndEmcBump array!" << endl;
    return kERROR;
  }

  fClusterArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcCluster"));
  if ( ! fClusterArray ) {
    cout << "-W- McBumpMatch::Init: "
	 << "No PndEmcCluster array!" << endl;
    return kERROR;
  }

  fDigiArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcDigi"));
  if ( ! fDigiArray ) {
    cout << "-W- McBumpMatch::Init: "
	 << "No PndEmcDigi array!" << endl;
    return kERROR;
  }

  fHitArray = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcHit"));
  if ( ! fHitArray ) {
    cout << "-W- McBumpMatch::Init: "
	 << "No PndEmcHit array!" << endl;
    return kERROR;
  }

  fChargedCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidChargedCand"));
  if ( ! fChargedCandidateArray ) {
    cout << "-W- McBumpMatch::Init: "
	 << "No PidChargedCand array!" << endl;
    return kERROR;
  }

  fNeutralCandidateArray = dynamic_cast<TClonesArray *> (ioman->GetObject("PidNeutralCand"));
  if ( ! fNeutralCandidateArray ) {
    cout << "-W- McBumpMatch::Init: "
	 << "No PidNeutralCand array!" << endl;
    return kERROR;
  }

  fMcArray = dynamic_cast<TClonesArray *> (ioman->GetObject("MCTrack"));
  if ( ! fMcArray ) {
    cout << "-W- McBumpMatch::Init: "
	 << "No McTrack array!" << endl;
    return kERROR;
  }

  return kSUCCESS;
}

void McBumpMatch::Exec(Option_t* opt)
{
  int verb = 0;
  if (verb>0||(nEvt%1000==0))
    cout << "==== New Event " << nEvt << " ===== " << endl;
  nChCand = fChargedCandidateArray->GetEntriesFast();
  nNeutCand = fNeutralCandidateArray->GetEntriesFast();
  nMcTrack = fMcArray->GetEntriesFast();
  if (split_detector())
    print_cands();
  ++nEvt;
  return;
}

inline
bool bpr_comp_pair(pair<int, pair<int,int> > i, pair<int,pair<int,int> > j) { return (i.first>j.first); }

inline
void McBumpMatch::print_vect(vector<int> &v) {
  cout << "N= " << v.size() << ": ";
  for (int i=0; i<v.size(); ++i) cout << v[i] << ", ";
  cout << endl;
}

/**
 * @breif MCTrack-Bump association based on main-contributor ordering
 * Function that associates a list of PndMCTrack with a list of PndEmcBump in decreasing order of contribution
 * Ie. The MC track whose descendents contribute the most to a given bump will be assoicated to the bump
 * An MC track is a contributor to a bump if it is part of the MC tracks that gnenerated one of the "points"
 * that ultimately made it into the bump
 * The best matching MC track id is set to the MCTruth of the Bump.
 * If none of the MC tracks in the list
 * contribute to a given Bump, its MCTruth is set to (some-crazy-number) to communicate with caller that
 * no match has been found. In the same way,
 * @param vector<PndEmcBump*> bump_list : list of bumps to be associated
 * @param vector<int> mc_list : list of MCTracks to be associated
 */
void McBumpMatch::mc_match_bumps_main_contrib(vector<PndEmcBump*> bump_list, vector<int>& mc_list) {
  vector<vector<int> > ancestory;
  for (int imc=0; imc < mc_list.size(); ++imc) {
    vector<int> a;
    find_ancestory(mc_list[imc], a);
    ancestory.push_back(a);
  }
  vector<vector<int> > contributors;
  for (int ibump=0; ibump < bump_list.size(); ++ibump) {
    vector<int> a;
    find_contributors(bump_list[ibump], a);
    contributors.push_back(a);
  }
  vector< pair<int, pair<int,int> > > score(bump_list.size()*mc_list.size());
  for (int imc = 0; imc < mc_list.size(); ++imc) {
    for (int ibump = 0; ibump < bump_list.size(); ++ibump) {
      vector<int> intersn;
      set_intersection(ancestory[imc].begin(),ancestory[imc].end(),
		       contributors[ibump].begin(),contributors[ibump].end(),
		       back_inserter(intersn));
      score.push_back( make_pair(intersn.size(), make_pair(ibump,imc)) );
    }
  }
  vector<bool> matched_bump(bump_list.size(), false);
  vector<bool> matched_mc(mc_list.size(), false);
  sort(score.begin(),score.end(),bpr_comp_pair);
  for (int ii = 0; ii < score.size(); ++ii) {
    int ibump = score[ii].second.first;
    int imc = score[ii].second.second;
    if (matched_bump[ibump]||matched_mc[imc]) continue;
    //if (score[ii].first) break;
    //bump_list[ibump]->SetMCTruth(mc_list[imc]);
    matched_bump[ibump] = true;
    matched_mc[imc] = true;
  }
}

void McBumpMatch::find_contributors(PndEmcBump *bump, vector<int>& tree) {
  vector<int> tmp;
  int digisize = bump->DigiList().size();
  for (int id=0; id<digisize; ++id) {
    PndEmcDigi *digi = (PndEmcDigi*) fDigiArray->At(bump->DigiList()[id]);
    PndEmcHit *hit = (PndEmcHit*) fHitArray->At(digi->GetHitIndex());
    int mcsize = hit->GetMcList().size();
    for (int imc=0; imc<mcsize; ++imc) {
      tmp.push_back(hit->GetMcList()[imc]);
    }
  }
  sort(tmp.begin(),tmp.end());
}

void McBumpMatch::find_ancestory(const int& mcidx, vector<int>& tree) {
  vector<bool> added(nMcTrack,false);
  added[mcidx] = true;
  for (int ii = mcidx+1; ii < nMcTrack; ++ii) {
    PndMCTrack* mctrk = (PndMCTrack *) fMcArray->At(ii);
    int tmp_mid = mctrk->GetMotherID();
    while(!added[tmp_mid]) {
      PndMCTrack* tmp_mctrk = (PndMCTrack *) fMcArray->At(tmp_mid);
      tmp_mid = tmp_mctrk->GetMotherID();
      if (tmp_mid<mcidx) break;
    }
    if (tmp_mid>=mcidx) {
      tree.push_back(ii);
      added[ii] = true;
    }
  }
}

void McBumpMatch::
get_mc_brem_photons(const int &mcidx, vector<int>& brem_photons) {
  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack){
    PndMCTrack *McTrack = (PndMCTrack *) fMcArray->At(iMcTrack);
    if (McTrack->GetPdgCode()!=22) continue; // only interested in photons
    if (McTrack->GetMotherID()!=mcidx) continue; // only photons emerging from primary electron
    const double radiusT = McTrack->GetStartVertex().Pt();
    if (radiusT>radMaxTracking_cm) continue;
    brem_photons.push_back(iMcTrack);
  }
}

void McBumpMatch::print_mc() {
  cout << "------ nMcTrack= " << nMcTrack << " ------" << endl;
  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack){
    PndMCTrack *McTrack = (PndMCTrack *) fMcArray->At(iMcTrack);
    int motherid = McTrack->GetMotherID();
    cout << "id= " << iMcTrack << " mId= " << motherid << endl;
  }
}

bool McBumpMatch::split_detector() {
  int nbump = fBumpArray->GetEntriesFast();
  for (int ibump=0; ibump<nbump; ++ibump) {
    PndEmcBump *bumpi = (PndEmcBump*) fBumpArray->At(ibump);
    for (int jbump=ibump+1; jbump<nbump; ++jbump) {
      PndEmcBump *bumpj = (PndEmcBump*) fBumpArray->At(jbump);
      if (bumpi->GetClusterIndex()==bumpj->GetClusterIndex()) return true;
    }
  }
  return false;
}

void McBumpMatch::print_cands() {

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
    cout << "ClustId = " << ii;
    cout << " Ang=(" << clust->position().Theta()*TMath::RadToDeg();
    cout << ", " << clust->position().Phi()*TMath::RadToDeg() << ") ";
    cout << " E= " << clust->energy();
    cout << " Nbump= " << clust->NBumps();

    int mcsize = clust->GetMcSize();
    cout << " McSize= " << mcsize << " McList: (";
    for (int i=0; i<mcsize; ++i) {
      cout << clust->GetMcIndex(i) << (i<(mcsize-1)?", ": "");
    }
    cout << ")";

    int digisize = clust->DigiList().size();
    cout << " Digi= " << digisize << "(";
    for (int i=0; i<digisize; ++i) {
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

    int mcsize = bump->GetMcSize();
    cout << " McSize= " << mcsize << " McList: (";
    for (int i=0; i<mcsize; ++i) {
      cout << bump->GetMcIndex(i) << (i<(mcsize-1)?", ": "");
    }
    cout << ")";

    int digisize = bump->DigiList().size();
    cout << " NDigi= " << digisize << "(";
    for (int i=0; i<digisize; ++i) {
      cout << bump->DigiList()[i] << (i<(digisize-1)?", ": ")");
    }
    cout << ")";
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

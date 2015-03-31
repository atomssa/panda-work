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

#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>

using namespace std;

BremPidReader::BremPidReader():
  fClusterArray(0), fPhiBumpArray(0), fBumpArray(0), fChargedCandidateArray(0), fNeutralCandidateArray(0), fBremCorrected4MomArray(0),fRecMomOfEle(0), fRecThetaOfEle(0), fRecPhiOfEle(0), fCharge(0), /*fSepPhotonE(0.),*/ fMergPhotonE(0.), fPersistance(kTRUE)
{
  output_name = "bremcorr.root";
  nEvt = 0;
  radMaxTracking_cm = 42.0;
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

  t->Branch("nch",&nch,"nch/I");

  t->Branch("charge",&charge,"charge[nch]/I");
  t->Branch("mom_mc",&mom_mc,"mom_mc[nch]/F");
  t->Branch("mom_rec",&mom_rec,"mom_rec[nch]/F");
  t->Branch("mom_cor",&mom_cor,"mom_cor[nch]/F");
  t->Branch("mom_wcor",&mom_wcor,"mom_wcor[nch]/F");
  t->Branch("mom_mrg",&mom_mrg,"mom_mrg[nch]/F");
  t->Branch("mom_sep",&mom_sep,"mom_sep[nch]/F");
  t->Branch("mom_stored",&mom_stored,"mom_stored[nch]/F");
  t->Branch("phi",&phi,"phi[nch]/F");
  t->Branch("the",&the,"the[nch]/F");
  t->Branch("phi_mc",&phi_mc,"phi_mc[nch]/F");
  t->Branch("the_mc",&the_mc,"the_mc[nch]/F");
  t->Branch("pdg_mc",&pdg_mc,"pdg_mc[nch]/I");
  t->Branch("nphot_sep",&nphot_sep,"nphot_sep[nch]/I");
  t->Branch("nphot_mrg",&nphot_mrg,"nphot_mrg[nch]/I");
  t->Branch("is_prim",&is_prim,"is_prim[nch]/I");

  t->Branch("_nmcb",&_nmcb,"_nmcb[nch]/I"); // number of MC brem photons associated with this track with r<42cm
  t->Branch("imcb_s",&imcb_s,"imcb_s[nch]/I"); // Start index (inclusive) of MC brem photons associated with this tracks with r<42cm
  t->Branch("imcb_e",&imcb_e,"imcb_e[nch]/I"); // End index (inclusive) of MC brem photons associated with this tracks with r<42cm

  // Phtons originating from primary MC tracks with initial point r<42cm
  t->Branch("nmcb",&nmcb,"nmcb/I"); // number of mc brem photons with r<42cm for all tracks
  t->Branch("mcb_phi",&mcb_phi,"mcb_phi[nmcb]/F"); // True phi angle of this MC brem photon
  t->Branch("mcb_the",&mcb_the,"mcb_the[nmcb]/F"); // True phi angle of this MC brem photon
  t->Branch("mcb_ene",&mcb_ene,"mcb_ene[nmcb]/F"); // True energy of this MC brem photon
  t->Branch("mcb_rad",&mcb_rad,"mcb_rad[nmcb]/F"); // True creation radius of this MC brem photon
  t->Branch("mcb_zed",&mcb_zed,"mcb_zed[nmcb]/F"); // True creation radius of this MC brem photon
  t->Branch("mcb_match",&mcb_match,"mcb_match[nmcb]/I"); // The index of the best matching separated bump to this MC brem photon
  t->Branch("mcb_score",&mcb_score,"mcb_score[nmcb]/I"); // The number of common MC tracks with the matching separated bump
  t->Branch("mcb_match_ab",&mcb_match_ab,"mcb_match_ab[nmcb]/I"); // The index of the best matching separated bump to this MC brem photon
  t->Branch("mcb_score_ab",&mcb_score_ab,"mcb_score_ab[nmcb]/I"); // The number of common MC tracks with the matching separated bump

  // Separated bumps found within the brem cut window of a reconstructed track
  t->Branch("_nsb",&_nsb,"_nsb[nch]/I"); // number of separated bumps associated with this track
  t->Branch("isb_s",&isb_s,"isb_s[nch]/I"); // Start index (inclusive) of separate bumps associated with this track
  t->Branch("isb_e",&isb_e,"isb_e[nch]/I"); // End index (inclusive) of separate bumps associated with this track

  t->Branch("nsb",&nsb,"nsb/I"); // number of separated bumps found for all tracks
  t->Branch("sb_idx",&sb_idx,"sb_idx[nsb]/F"); // index of the separated bum found for this track in the list of all bumps
  t->Branch("sb_phi",&sb_phi,"sb_phi[nsb]/F"); // phi agnle of of this separated bump /* REDUNDENT WITH ab_phi of the right index */
  t->Branch("sb_the",&sb_the,"sb_the[nsb]/F"); // theta agnle of of this separated bump /* REDUNDENT WITH ab_the of the right index */
  t->Branch("sb_ene",&sb_ene,"sb_ene[nsb]/F"); // reconstructed energy of this separated bump  /* REDUNDENT WITH ab_ene of the right index */
  t->Branch("sb_rcalc",&sb_rcalc,"sb_rcalc[nsb]/F"); // The recalculated radius of this separated bump based on DeltaPhi and pT
  t->Branch("sb_match",&sb_match,"sb_match[nsb]/I"); // The index of the best matching MC brem photon to this separated bump
  t->Branch("sb_score",&sb_score,"sb_score[nsb]/I"); // The number of common MC tracks with the matching MC track

  // All bumps list
  t->Branch("nab",&nab,"nab/I");
  t->Branch("ab_phi",&ab_phi,"ab_phi[nab]/F"); // phi agnle of of this separated bump
  t->Branch("ab_the",&ab_the,"ab_the[nab]/F"); // theta agnle of of this separated bump
  t->Branch("ab_ene",&ab_ene,"ab_ene[nab]/F"); // reconstructed energy of this separated bump
  t->Branch("ab_isb",&ab_isb,"ab_isb[nab]/I"); // Index of the track for which this bump had been found to be a separted bump. -1 otherwise
  t->Branch("ab_ich",&ab_ich,"ab_ich[nab]/I"); // Index of the track with which this bump shares EmcIndex
  t->Branch("ab_match",&ab_match,"ab_match[nsb]/I"); // The index of the best matching MC brem photon to this bump
  t->Branch("ab_score",&ab_score,"ab_score[nsb]/I"); // The number of common MC tracks with the matching MC track

}

void BremPidReader::Exec(Option_t* opt)
{

  int verb = 0;

  if (verb>0||(nEvt%1000==0))
    cout << "==== New Event " << nEvt << " ===== " << endl;

  nChCand = fChargedCandidateArray->GetEntriesFast();
  nNeutCand = fNeutralCandidateArray->GetEntriesFast();
  nMcTrack = fMcArray->GetEntriesFast();

  nch = 0;
  nsb = 0;
  nmcb = 0;

  vector<vector<int> > ab_tree;
  fill_bump_list(ab_tree);

  for (int iCand = 0; iCand<nChCand; ++iCand){

    if (verb>1)
      cout << "--------------- New Reco Track " << iCand << "------------------" << endl;

    PndPidCandidate* theChargedCand = (PndPidCandidate*) fChargedCandidateArray->At(iCand);
    fRecMomOfEle = theChargedCand->GetMomentum().Mag();
    fRecThetaOfEle = theChargedCand->GetMomentum().Theta()*TMath::RadToDeg();
    fRecPhiOfEle = theChargedCand->GetMomentum().Phi()*TMath::RadToDeg();
    fCharge = theChargedCand->GetCharge();

    charge[nch] = fCharge;
    mom_rec[nch] = fRecMomOfEle;
    phi[nch] = fRecPhiOfEle;
    the[nch] = fRecThetaOfEle;

    // get the mctruth
    Int_t mcidx = theChargedCand->GetMcIndex();
    // The following line (condition on mcidx) should not be allowed in production code!!
    // it will lead to mismatch between PndPidCand idx and PndPidBremCorr Index
    if (mcidx>nMcTrack || mcidx<0) continue;
    PndMCTrack *truth = (PndMCTrack*) fMcArray->At(mcidx);

    int pdg = truth->GetPdgCode();
    int mid = truth->GetMotherID();
    int mpdg = -1;
    if (mid>=0){
      PndMCTrack *mother = (PndMCTrack*) fMcArray->At(mcidx);
      mpdg  = mother->GetPdgCode();
    }
    //if (pdg==11||pdg==-11) {
    //  if (mid!=-1&&mpdg!=433) continue; // keep only primary electrons. Its complicated as it is
    //} else if (pdg==211||pdg==-211) {
    //  if (mid!=-1) continue;
    //} else {
    //  continue;
    //}
    if (mid!=-1) continue;
    
    is_prim[nch] = truth->GetMotherID()==-1;
    mom_mc[nch] = truth->GetMomentum().Mag();
    phi_mc[nch] = truth->GetMomentum().Phi()*TMath::RadToDeg();
    the_mc[nch] = truth->GetMomentum().Theta()*TMath::RadToDeg();
    pdg_mc[nch] = truth->GetPdgCode();

    TVector3 mom = theChargedCand->GetMomentum();
    double ene = theChargedCand->GetEnergy();
    TLorentzVector m4 = TLorentzVector(mom,ene);
    double mass = 0.000511; // Electron mass hypothesis (GeV)

    //int nPhotSep = 0;
    //fSepPhotonE = GetSepPhotonE(theChargedCand, nPhotSep);

    //vector<int> _idx_sb;
    //vector<double> _rad_sb,_ene_sb,_phi_sb, _the_sb;
    //vector<vector<int> > sb_tree;
    //double sepPhotonE = 0, sepPhotonEWtd = 0;
    //GetSepPhotonE_fromBumps(theChargedCand, sepPhotonE, sepPhotonEWtd, _idx_sb, _rad_sb, _ene_sb, _phi_sb, _the_sb, sb_tree);
    //assert(_rad_sb.size()+nsb <= nsb_max);
    //nphot_sep[nch] = _rad_sb.size(); // Obsolete...
    //TVector3 momCorrSep = ((sepPhotonE+fRecMomOfEle)/fRecMomOfEle) * mom;
    //mom_sep[nch] = momCorrSep.Mag();
    //for (int i=0; i<_rad_sb.size(); ++i) {
    //  ab_isb[_idx_sb[nsb+i]] = nch;
    //  sb_idx[nsb+i] = _idx_sb[i];
    //  sb_phi[nsb+i] = _phi_sb[i];
    //  sb_the[nsb+i] = _the_sb[i];
    //  sb_rcalc[nsb+i] = _rad_sb[i];
    //  sb_ene[nsb+i] = _ene_sb[i];
    //}
    //isb_s[nch] = nsb;
    //nsb += _rad_sb.size();
    //isb_e[nch] = nsb;
    //_nsb[nch] = _rad_sb.size();

    vector<vector<int> > sb_tree;
    double sepPhotonE = 0, sepPhotonEWtd = 0;
    GetSepPhotonE_fromBumps(theChargedCand, sepPhotonE, sepPhotonEWtd, sb_tree);
    TVector3 momCorrSep = ((sepPhotonE+fRecMomOfEle)/fRecMomOfEle) * mom;
    mom_sep[nch] = momCorrSep.Mag();

    int nPhotMrg = 0;
    fMergPhotonE = GetMergPhotonE(theChargedCand, nPhotMrg);
    nphot_mrg[nch] = nPhotMrg;
    TVector3 momCorrMrg = ((fMergPhotonE+fRecMomOfEle)/fRecMomOfEle) * mom;
    mom_mrg[nch] = momCorrMrg.Mag();

    const double energy_gamma = sepPhotonE + fMergPhotonE;
    TVector3 momCorr = ((energy_gamma+fRecMomOfEle)/fRecMomOfEle) * mom;
    //const double eneCorr = TMath::Hypot(mass, momCorr.Mag());
    mom_cor[nch] = momCorr.Mag();

    const double energy_gamma_wtd = sepPhotonEWtd + fMergPhotonE;
    TVector3 momCorrWtd = ((energy_gamma_wtd+fRecMomOfEle)/fRecMomOfEle) * mom;
    //const double eneCorrWtd = TMath::Hypot(mass, momCorrWtd.Mag());
    mom_wcor[nch] = momCorrWtd.Mag();

    PndPidBremCorrected4Mom *tmp = (PndPidBremCorrected4Mom*)fBremCorrected4MomArray->At(iCand);
    if(tmp) {
      TVector3 storedCorrMom = tmp->GetMomentum();
      mom_stored[nch] = storedCorrMom.Mag();
    } else {
      mom_stored[nch] = -9999.;
    }

    //vector<double> _rad_mc, _zed_mc, _ene_mc, _phi_mc, _the_mc;
    //vector<vector<int> > mcb_tree;
    //get_mc_brem_photons(mcidx, _rad_mc, _zed_mc, _ene_mc, _phi_mc, _the_mc, mcb_tree);
    //assert(_rad_mc.size()+nmcb < nmcb_max);
    //for (int ii = 0; ii < _rad_mc.size(); ++ii) {
    //  mcb_phi[nmcb+ii] = _phi_mc[ii];
    //  mcb_the[nmcb+ii] = _the_mc[ii];
    //  mcb_rad[nmcb+ii] = _rad_mc[ii];
    //  mcb_zed[nmcb+ii] = _zed_mc[ii];
    //  mcb_ene[nmcb+ii] = _ene_mc[ii];
    //}
    //imcb_s[nch] = nmcb;
    //nmcb += _rad_mc.size();
    //imcb_e[nch] = nmcb;
    //_nmcb[nch] = _rad_mc.size();

    vector<vector<int> > mcb_tree;
    get_mc_brem_photons(mcidx, mcb_tree);

    //cout << "Matching for track "  << nch << endl;
    // Do matching between separate bump and mc-brem photon
    vector<int> _sb_match(_nsb[nch],-1), _mcb_match(_nmcb[nch],-1);
    vector<int> _sb_score(_nsb[nch],-1), _mcb_score(_nmcb[nch],-1);
    brem_matching(sb_tree, mcb_tree, _sb_match, _mcb_match, _sb_score, _mcb_score);
    for (int ii = imcb_s[nch]; ii < imcb_e[nch]; ++ii) {
      mcb_match[ii] = isb_s[nch] + _mcb_match.at(ii-imcb_s[nch]);
      mcb_score[ii] = _mcb_score.at(ii-imcb_s[nch]);
    }
    for (int ii = isb_s[nch]; ii < isb_e[nch]; ++ii) {
      sb_match[ii] = imcb_s[nch] + _sb_match.at(ii-isb_s[nch]);
      sb_score[ii] = _sb_score.at(ii-isb_s[nch]);
      if (sb_match[ii] > nmcb) { cout << "DISASTER!!" << endl;}
    }

    // Do matching between all bumps and mc-brem photons
    vector<int> _ab_match(nab,-1), _mcb_match_ab(_nmcb[nch],-1);
    vector<int> _ab_score(nab,-1), _mcb_score_ab(_nmcb[nch],-1);
    brem_matching(ab_tree, mcb_tree, _ab_match, _mcb_match_ab, _ab_score, _mcb_score_ab);
    for (int ii = imcb_s[nch]; ii < imcb_e[nch]; ++ii) {
      mcb_match_ab[ii] = _mcb_match_ab.at(ii-imcb_s[nch]);
      mcb_score_ab[ii] = _mcb_score_ab.at(ii-imcb_s[nch]);
    }
    for (int ii = 0; ii < nab; ++ii) {
      if (_ab_match.at(ii)<0) continue;
      ab_match[ii] = imcb_s[nch] + _ab_match.at(ii);
      ab_score[ii] = _ab_score.at(ii);
      if (ab_match[ii] > nmcb) { cout << "DISASTER!!" << endl;}
    }

    ++nch;

  }

  //print_cands();
  t->Fill();

  nEvt++;

}

inline
bool bpr_comp_pair(pair<int, pair<int,int> > i, pair<int,pair<int,int> > j) { return (i.first>j.first); }

inline
void BremPidReader::print_vect(vector<int> &v) {
  cout << "N= " << v.size() << ": ";
  for (int i=0; i<v.size(); ++i) cout << v[i] << ", ";
  cout << endl;
}
void BremPidReader::brem_matching(vector<vector<int> >& _sb_tree, vector<vector<int> >&_mcb_tree,
				  vector<int>& _sb_match, vector<int>& _mcb_match,
				  vector<int>& _sb_score, vector<int>& _mcb_score) {
  //cout << "--------------------------------------------" << endl;
  int n_sb = _sb_tree.size();
  int n_mcb = _mcb_tree.size();
  vector< pair<int, pair<int,int> > > score(n_sb*n_mcb);
  for (int isb = 0; isb < n_sb; ++isb) {
    for (int imcb = 0; imcb < n_mcb; ++imcb) {
      vector<int> intersn;
      //cout << "--------------------------------------------" << endl;
      set_intersection(_sb_tree[isb].begin(),_sb_tree[isb].end(),_mcb_tree[imcb].begin(),_mcb_tree[imcb].end(),back_inserter(intersn));
      //cout << "Common (score= "<< intersn.size()<<  "): ";
      //print_vect(intersn);
      score.push_back( make_pair(intersn.size(), make_pair(isb,imcb)) );
    }
  }
  sort(score.begin(),score.end(),bpr_comp_pair);
  //cout << " n_sb= " << n_sb << " n_mcb= " << n_mcb << endl;
  for (int ii = 0; ii < score.size(); ++ii) {
    int isb = score[ii].second.first;
    int imcb = score[ii].second.second;
    if (_sb_match[isb]!=-1||_mcb_match[imcb]!=-1) continue;
    //cout << "score= " << score[ii].first << " isb= " << isb << " imcb= " << imcb << endl;
    _sb_match[isb] = imcb;
    _mcb_match[imcb] = isb;
    _mcb_score[imcb] = _sb_score[isb] = score[ii].first;
  }
}

void BremPidReader::find_ancestory(const int& mcidx, vector<int>& tree) {

  std::vector<int> fMotherIds;
  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack)
    fMotherIds.push_back(((PndMCTrack *) fMcArray->At(iMcTrack))->GetMotherID());

  vector<int> added(nMcTrack,0);
  added[mcidx] = 1;
  for (int ii = mcidx+1; ii < nMcTrack; ++ii) {
    int tmp_mid = fMotherIds[ii];
    while(added[tmp_mid]==0) {
      tmp_mid = fMotherIds[tmp_mid];
      if (tmp_mid<mcidx) break;
    }
    if (tmp_mid>=mcidx) {
      tree.push_back(ii);
      added[ii] = 1;
    }
  }
}

void BremPidReader::
get_mc_brem_photons(const int &mcidx, vector<vector<int> >& _mcb_tree) {
  imcb_s[nch] = nmcb;
  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack){
    PndMCTrack *McTrack = (PndMCTrack *) fMcArray->At(iMcTrack);
    if (McTrack->GetPdgCode()!=22) continue; // only interested in photons
    if (McTrack->GetMotherID()!=mcidx) continue; // only photons emerging from primary electron
    const double radiusT = McTrack->GetStartVertex().Pt();
    if (radiusT>radMaxTracking_cm) continue;
    assert(nmcb < nmcb_max);
    mcb_rad[nmcb] = radiusT;
    mcb_zed[nmcb] = McTrack->GetStartVertex().Z();
    mcb_ene[nmcb] = McTrack->GetMomentum().Mag();
    mcb_phi[nmcb] = McTrack->GetMomentum().Phi()*TMath::RadToDeg();;
    mcb_the[nmcb] = McTrack->GetMomentum().Theta()*TMath::RadToDeg();;
    vector<int> tree;
    find_ancestory(iMcTrack,tree);
    _mcb_tree.push_back(tree);
    nmcb++;

  }
  imcb_e[nch] = nmcb;
  _nmcb[nch] = imcb_e[nch] - imcb_s[nch];
}

//void BremPidReader::
//get_mc_brem_photons(const int &mcidx, vector<double> &_rad, vector<double> &_zed,
//		    vector<double> &_ene, vector<double> &_phi, vector<double> &the, vector<vector<int> >& _mcb_tree) {
//  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack){
//    PndMCTrack *McTrack = (PndMCTrack *) fMcArray->At(iMcTrack);
//    if (McTrack->GetPdgCode()!=22) continue; // only interested in photons
//    if (McTrack->GetMotherID()!=mcidx) continue; // only photons emerging from primary electron
//    const double radiusT = McTrack->GetStartVertex().Pt();
//    const double zed = McTrack->GetStartVertex().Z();
//    if (radiusT>radMaxTracking_cm) continue;
//    _rad.push_back(radiusT);
//    _zed.push_back(zed);
//    // temporary BS to make sure kin is okay
//    _ene.push_back(McTrack->GetMomentum().Mag());
//    _phi.push_back(McTrack->GetMomentum().Phi());
//    _the.push_back(McTrack->GetMomentum().Theta());
//    vector<int> tree;
//    find_ancestory(iMcTrack,tree);
//    _mcb_tree.push_back(tree);
//  }
//}

PndPidBremCorrected4Mom* BremPidReader::AddBremCorrected4Mom(){
  TClonesArray& clref = *fBremCorrected4MomArray;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) PndPidBremCorrected4Mom();
}

void BremPidReader::
fill_bump_list(vector<vector<int> > &_ab_tree) {
  nab = fBumpArray->GetEntriesFast();
  for(Int_t iBump = 0; iBump<nab; ++iBump) {
    PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(iBump);
    assert(iBump<nab_max);
    int _ich = -1;
    for (int iCand = 0; iCand<nChCand; ++iCand){
      PndPidCandidate* theChargedCand = (PndPidCandidate*) fChargedCandidateArray->At(iCand);
      if ( PhotonBump->GetClusterIndex() == theChargedCand->GetEmcIndex() ) {
	_ich = iCand;
	break; // assumes a bump can be associated with only one charge candidate
      }
    }
    ab_ene[iBump] = PhotonBump->GetEnergyCorrected();
    ab_the[iBump] = PhotonBump->position().Theta()*TMath::RadToDeg();
    ab_phi[iBump] = PhotonBump->position().Phi()*TMath::RadToDeg();
    ab_isb[iBump] = -1;
    ab_ich[iBump] = _ich;

    const Int_t iSepClust = PhotonBump->GetClusterIndex();
    PndEmcCluster *PhotonCluster = (PndEmcCluster*) fClusterArray->At(iSepClust);
    vector<int> tmp;
    int digisize = PhotonCluster->DigiList().size();
    for (int id=0; id<digisize; ++id) {
      PndEmcDigi *digi = (PndEmcDigi*) fDigiArray->At(PhotonCluster->DigiList()[id]);
      PndEmcHit *hit = (PndEmcHit*) fHitArray->At(digi->GetHitIndex());
      int mcsize = hit->GetMcList().size();
      for (int imc=0; imc<mcsize; ++imc) {
	tmp.push_back(hit->GetMcList()[imc]);
      }
    }
    sort(tmp.begin(),tmp.end());
    std::vector<int>::iterator it = unique(tmp.begin(),tmp.end());
    tmp.resize(distance(tmp.begin(),it));
    _ab_tree.push_back(tmp);

  }
}

double BremPidReader::
GetSepPhotonE_fromBumps(PndPidCandidate *ChargedCand, double &esep, double &esep_wtd, vector<vector<int> > &_sb_tree) {

  Float_t PhotonTotEnergySep = 0;
  Float_t PhotonTotEnergySepWtd = 0;
  isb_s[nch] = nsb;
  const int nBump = fBumpArray->GetEntriesFast();
  for(Int_t iBump = 0; iBump<nBump; ++iBump)
    {
      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(iBump);
      const Float_t PhotonEnergySep = PhotonBump->GetEnergyCorrected();

      const Int_t iSepClust = PhotonBump->GetClusterIndex();
      if ( PhotonBump->GetClusterIndex() == ChargedCand->GetEmcIndex() ) continue;
      //if ( PhotonEnergySep > 0.8* ChargedCand->GetEnergy() ) continue;

      PndEmcCluster *PhotonCluster = (PndEmcCluster*) fClusterArray->At(iSepClust);

      const double PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
      const double PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();

      const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
      const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
      const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;

      const Float_t RealDeltaPhi = fCharge<0?PhotonPhiSep-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiSep;
      const Float_t RealDeltaTheta = fCharge<0?PhotonThetaSep-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaSep;

      const Float_t RealDeltaPhiRad = RealDeltaPhi*TMath::DegToRad();
      const Float_t rad_calc = 100*TMath::Sin(RealDeltaPhiRad/2.)*2*Pt/0.3/2.0; // B=2T

      //const Float_t wt = -rad_calc/42. + 1.;
      const Float_t wt = 1.0/(1.+TMath::Exp((rad_calc-21.)/5));
      const Float_t ThetaCutUp = 2.;
      const Float_t ThetaCutDown = -2.;
      const Float_t PhiCutUp = (fRecThetaOfEle <= 23.)?DeltaPhiForward:DeltaPhiBarrel;
      const Float_t PhiCutDown = -1;

      const Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
      const Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;

      if (PhiCut && ThetaCut) {
	assert(nmcb < nmcb_max);
	ab_isb[iBump] = nch;
	sb_idx[nsb] = iBump;
	sb_phi[nsb] = PhotonPhiSep;
	sb_the[nsb] = PhotonThetaSep;
	sb_rcalc[nsb] = rad_calc;
	sb_ene[nsb] = PhotonEnergySep;
	nsb++;

	vector<int> tmp;
	int digisize = PhotonCluster->DigiList().size();
	for (int id=0; id<digisize; ++id) {
	  PndEmcDigi *digi = (PndEmcDigi*) fDigiArray->At(PhotonCluster->DigiList()[id]);
	  PndEmcHit *hit = (PndEmcHit*) fHitArray->At(digi->GetHitIndex());
	  int mcsize = hit->GetMcList().size();
	  for (int imc=0; imc<mcsize; ++imc) {
	    tmp.push_back(hit->GetMcList()[imc]);
	  }
	}
	sort(tmp.begin(),tmp.end());
	std::vector<int>::iterator it = unique(tmp.begin(),tmp.end());
	tmp.resize(distance(tmp.begin(),it));
	_sb_tree.push_back(tmp);

	PhotonTotEnergySep += PhotonEnergySep;
	PhotonTotEnergySepWtd += wt*PhotonEnergySep;

      }

    }//loop neutralcand
  isb_e[nch] = nsb;
  nphot_sep[nch] = _nsb[nch] = isb_e[nch] - isb_s[nch];

  if (PhotonTotEnergySep < fRecMomOfEle/100.) PhotonTotEnergySep = 0;
  if (PhotonTotEnergySepWtd < fRecMomOfEle/100.) PhotonTotEnergySepWtd = 0;

  esep = PhotonTotEnergySep;
  esep_wtd = PhotonTotEnergySepWtd;

}

//double BremPidReader::
//GetSepPhotonE_fromBumps(PndPidCandidate *ChargedCand, double &esep, double &esep_wtd, vector<int>& _idx,
//			vector<double>& _rcalc, vector<double>& _ene, vector<double>& _phi,
//			vector<double>& _the, vector<vector<int> > &_sb_tree) {
//
//  Float_t PhotonTotEnergySep = 0;
//  Float_t PhotonTotEnergySepWtd = 0;
//  const int nBump = fBumpArray->GetEntriesFast();
//  for(Int_t iBump = 0; iBump<nBump; ++iBump)
//    {
//      PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(iBump);
//      const Float_t PhotonEnergySep = PhotonBump->GetEnergyCorrected();
//
//      const Int_t iSepClust = PhotonBump->GetClusterIndex();
//      if ( PhotonBump->GetClusterIndex() == ChargedCand->GetEmcIndex() ) continue;
//      //if ( PhotonEnergySep > 0.8* ChargedCand->GetEnergy() ) continue;
//
//      PndEmcCluster *PhotonCluster = (PndEmcCluster*) fClusterArray->At(iSepClust);
//
//      const double PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
//      const double PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();
//
//      const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
//      const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
//      const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;
//
//      const Float_t RealDeltaPhi = fCharge<0?PhotonPhiSep-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiSep;
//      const Float_t RealDeltaTheta = fCharge<0?PhotonThetaSep-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaSep;
//
//      const Float_t RealDeltaPhiRad = RealDeltaPhi*TMath::DegToRad();
//      const Float_t rad_calc = 100*TMath::Sin(RealDeltaPhiRad/2.)*2*Pt/0.3/2.0;
//
//      //const Float_t wt = -rad_calc/42. + 1.;
//      const Float_t wt = 1.0/(1.+TMath::Exp((rad_calc-21.)/5));
//      const Float_t ThetaCutUp = 2.;
//      const Float_t ThetaCutDown = -2.;
//      const Float_t PhiCutUp = (fRecThetaOfEle <= 23.)?DeltaPhiForward:DeltaPhiBarrel;
//      const Float_t PhiCutDown = -1;
//
//      const Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
//      const Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
//
//      if (PhiCut && ThetaCut) {
//	_idx.push_back(iBump);
//	_rcalc.push_back(rad_calc);
//	_ene.push_back(PhotonEnergySep);
//	_phi.push_back(PhotonPhiSep);
//	_the.push_back(PhotonThetaSep);
//
//	vector<int> tmp;
//	int digisize = PhotonCluster->DigiList().size();
//	for (int id=0; id<digisize; ++id) {
//	  PndEmcDigi *digi = (PndEmcDigi*) fDigiArray->At(PhotonCluster->DigiList()[id]);
//	  PndEmcHit *hit = (PndEmcHit*) fHitArray->At(digi->GetHitIndex());
//	  int mcsize = hit->GetMcList().size();
//	  for (int imc=0; imc<mcsize; ++imc) {
//	    tmp.push_back(hit->GetMcList()[imc]);
//	  }
//	}
//	_sb_tree.push_back(tmp);
//	sort(tmp.begin(),tmp.end());
//	std::vector<int>::iterator it = unique(tmp.begin(),tmp.end());
//	tmp.resize(distance(tmp.begin(),it));
//	PhotonTotEnergySep += PhotonEnergySep;
//	PhotonTotEnergySepWtd += wt*PhotonEnergySep;
//      }
//
//    }//loop neutralcand
//
//  if (PhotonTotEnergySep < fRecMomOfEle/100.) PhotonTotEnergySep = 0;
//  if (PhotonTotEnergySepWtd < fRecMomOfEle/100.) PhotonTotEnergySepWtd = 0;
//
//  esep = PhotonTotEnergySep;
//  esep_wtd = PhotonTotEnergySepWtd;
//
//}

double BremPidReader::GetMergPhotonE(PndPidCandidate *ChargedCand, int &nphotmrg){

  Double_t PhotonTotEnergyMerg = 0.0;

  // no EMcal cluster associated with track ...
  if (ChargedCand->GetEmcIndex() < 0) return 0.0;

  PndEmcBump *EleBump = (PndEmcBump *) fBumpArray->At(ChargedCand->GetEmcIndex());
  Int_t EleRefCluster = EleBump->GetClusterIndex();

  if (EleRefCluster < 0) return 0.0;

  std::vector<PndEmcBump*> EmcPhiBumpList;
  int nPhiBump = fPhiBumpArray->GetEntriesFast();
  for (int ipb=0; ipb<nPhiBump; ++ipb) {
    PndEmcBump *phibump = (PndEmcBump*) fPhiBumpArray->At(ipb);
    if ( phibump->GetClusterIndex() == EleRefCluster ) {
      EmcPhiBumpList.push_back(phibump);
    }
  }

  // This WORKS
  int iMax = 0;
  Float_t eMax = -1e9;
  //for (Int_t ib = EmcPhiBumpList.size()-1; ib >= 0; --ib) {
  for (Int_t ib = 0; ib < EmcPhiBumpList.size(); ++ib) {
    if( EmcPhiBumpList[ib]->energy() > eMax) {
      iMax = ib;
      eMax = EmcPhiBumpList[ib]->energy();
    }
  }
  Int_t iS = fCharge<0?0:iMax+1;
  Int_t iE = fCharge<0?iMax-1:EmcPhiBumpList.size()-1;
  for(Int_t r = iS; r<=iE; r++) {
    PhotonTotEnergyMerg += EmcPhiBumpList[r]->energy();
    nphotmrg++;
  }

  //Double_t EnergyCut = 0.15/TMath::Sin(fRecThetaOfEle*TMath::DegToRad());
  //if (fCharge<0) {
  //  for (Int_t i_phibump = EmcPhiBumpList.size()-1; i_phibump >= 0; --i_phibump) {
  //    if( EmcPhiBumpList[i_phibump]->energy() > EnergyCut) {
  //	for(Int_t r = 0; r<i_phibump; r++) {
  //	  PhotonTotEnergyMerg += EmcPhiBumpList[r]->energy();
  //	  nphotmrg++;
  //	}
  //	i_phibump = -1;
  //    }
  //  }
  //} else {
  //  for (Int_t i_phibump =0; i_phibump < EmcPhiBumpList.size(); ++i_phibump) {
  //    if( EmcPhiBumpList[i_phibump]->energy() > EnergyCut) {
  //	for(Int_t r = i_phibump+1; r<EmcPhiBumpList.size(); r++) {
  //	  PhotonTotEnergyMerg += EmcPhiBumpList[r]->energy();
  //	  nphotmrg++;
  //	}
  //	i_phibump = EmcPhiBumpList.size()+1;
  //    }
  //  }
  //}

  if(PhotonTotEnergyMerg > fRecMomOfEle/100.)  {
    return PhotonTotEnergyMerg;
  } else {
    return 0.0;
  }

}

void BremPidReader::FinishTask() {
  f->cd();
  t->Write();
}

void BremPidReader::print_mc() {
  cout << "------ nMcTrack= " << nMcTrack << " ------" << endl;
  for (int iMcTrack = 0; iMcTrack<nMcTrack; ++iMcTrack){
    PndMCTrack *McTrack = (PndMCTrack *) fMcArray->At(iMcTrack);
    int motherid = McTrack->GetMotherID();
    cout << "id= " << iMcTrack << " mId= " << motherid << endl;
  }
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

double BremPidReader::GetSepPhotonE(PndPidCandidate *ChargedCand, int &nphotsep){

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

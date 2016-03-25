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

BremPidReader::BremPidReader(int ibinsize):
  iBinSize(ibinsize), fClusterArray(0), fPhiBumpArray(), fBumpArray(0),
  fChargedCandidateArray(0), fNeutralCandidateArray(0), fBremCorrected4MomArray(0),
  fRecMomOfEle(0), fRecThetaOfEle(0), fRecPhiOfEle(0), fCharge(0), fPersistance(kTRUE)
{
  output_name = "bremcorr.root";
  nEvt = 0;
  radMaxTracking_cm = 42.0;
  if (iBinSize>nbs) iBinSize=0;
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

  fMCHeader = dynamic_cast<FairMCEventHeader*>(ioman->GetObject("MCEventHeader."));
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

  //const char* tca_name = (iBinSize==0?"EmcPhiBump":Form("EmcPhiBump%d",iBinSize));
  //cout << "BremPidReader::Init: using tca " << tca_name << " as phibump data source" << endl;
  //fPhiBumpArray = dynamic_cast<TClonesArray *> (ioman->GetObject(tca_name));
  //if ( ! fPhiBumpArray ) {
  //  cout << "-W- PndEmcMakePhiBump::Init: "
  //	 << "No " << tca_name << " array!" << endl;
  //  return kERROR;
  //}

  for (int ibs=0; ibs<nbs; ++ibs) {
    const char* tca_name = (ibs==0?"EmcPhiBump":Form("EmcPhiBump%d",ibs));
    cout << "BremPidReader::Init: using tca " << tca_name << " as phibump data source" << endl;
    fPhiBumpArray[ibs] = dynamic_cast<TClonesArray *> (ioman->GetObject(tca_name));
    if ( ! fPhiBumpArray[ibs] ) {
      cout << "-W- PndEmcMakePhiBump::Init: "
	   << "No " << tca_name << " array!" << endl;
      return kERROR;
    }
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
  t->Branch("mom_sep",&mom_sep,"mom_sep[nch]/F");
  t->Branch("mom_stored",&mom_stored,"mom_stored[nch]/F");
  t->Branch("mom_sep_w",&mom_sep_w,"mom_sep_w[nch]/F");
  t->Branch("mom_sep_w_bf",&mom_sep_w_bf,"mom_sep_w_bf[nch]/F");

  for (int ibs=0; ibs< nbs; ++ibs) {
    t->Branch(Form("b%d_mom_cor",ibs),&mom_cor[ibs],Form("b%d_mom_cor[nch]/F",ibs));
    t->Branch(Form("b%d_mom_wcor",ibs),&mom_wcor[ibs],Form("b%d_mom_wcor[nch]/F",ibs));
    t->Branch(Form("b%d_mom_mrg",ibs),&mom_mrg[ibs],Form("b%d_mom_mrg[nch]/F",ibs));
    t->Branch(Form("b%d_mom_mrg_w",ibs),&mom_mrg_w[ibs],Form("b%d_mom_mrg_w[nch]/F",ibs));
    t->Branch(Form("b%d_mom_mrg_w_bf",ibs),&mom_mrg_w_bf[ibs],Form("b%d_mom_mrg_w_bf[nch]/F",ibs));
    t->Branch(Form("b%d_mom_mrg_pc",ibs),&mom_mrg_pc[ibs],Form("b%d_mom_mrg_pc[nch]/F",ibs));
    t->Branch(Form("b%d_mom_mrg_w_pc",ibs),&mom_mrg_w_pc[ibs],Form("b%d_mom_mrg_w_pc[nch]/F",ibs));
    t->Branch(Form("b%d_mom_mrg_w_bf_pc",ibs),&mom_mrg_w_bf_pc[ibs],Form("b%d_mom_mrg_w_bf_pc[nch]/F",ibs));
    t->Branch(Form("b%d_mom_wbfcor",ibs),&mom_wbfcor[ibs],Form("b%d_mom_wbfcor[nch]/F",ibs));
    t->Branch(Form("b%d_mom_wbfcor_mw_bf",ibs),&mom_wbfcor_mw_bf[ibs],Form("b%d_mom_wbfcor_mw_bf[nch]/F",ibs));
    t->Branch(Form("b%d_mom_wbfcor_mw_bf_pc",ibs),&mom_wbfcor_mw_bf_pc[ibs],Form("b%d_mom_wbfcor_mw_bf_pc[nch]/F",ibs));
  }

  t->Branch("phi",&phi,"phi[nch]/F");
  t->Branch("the",&the,"the[nch]/F");
  t->Branch("phi_mc",&phi_mc,"phi_mc[nch]/F");
  t->Branch("the_mc",&the_mc,"the_mc[nch]/F");
  t->Branch("pdg_mc",&pdg_mc,"pdg_mc[nch]/I");
  t->Branch("nphot_sep",&nphot_sep,"nphot_sep[nch]/I");
  t->Branch("nphot_mrg",&nphot_mrg,"nphot_mrg[nch]/I");
  t->Branch("nphot_mrg_pc",&nphot_mrg_pc,"nphot_mrg_pc[nch]/I");
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
  t->Branch("sb_idx",&sb_idx,"sb_idx[nsb]/F"); // index of the separated bump found for all tracks in the list of all bumps
  t->Branch("sb_phi",&sb_phi,"sb_phi[nsb]/F"); // phi agnle of of this separated bump /* REDUNDENT WITH ab_phi of the right index */
  t->Branch("sb_the",&sb_the,"sb_the[nsb]/F"); // theta agnle of of this separated bump /* REDUNDENT WITH ab_the of the right index */
  t->Branch("sb_ene",&sb_ene,"sb_ene[nsb]/F"); // reconstructed energy of this separated bump  /* REDUNDENT WITH ab_ene of the right index */
  t->Branch("sb_rcalc",&sb_rcalc,"sb_rcalc[nsb]/F"); // The recalculated radius of this separated bump based on DeltaPhi and pT
  t->Branch("sb_zcalc",&sb_zcalc,"sb_zcalc[nsb]/F"); // The recalculated radius of this separated bump based on DeltaPhi and pT
  t->Branch("sb_match",&sb_match,"sb_match[nsb]/I"); // The index of the best matching MC brem photon to this separated bump
  t->Branch("sb_score",&sb_score,"sb_score[nsb]/I"); // The number of common MC tracks with the matching MC brem photon

  // Merged phibumps found in the electron's cluster
  for (int ibs=0; ibs<nbs; ++ibs) {
    t->Branch(Form("b%d__npb",ibs),&_npb[ibs],Form("b%d__npb[nch]/I",ibs)); // number of merged phi bumps associated with this track
    t->Branch(Form("b%d_ipb_s",ibs),&ipb_s[ibs],Form("b%d_ipb_s[nch]/I",ibs)); // Start index (inclusive) of merged phi bumps associated with this track
    t->Branch(Form("b%d_ipb_e",ibs),&ipb_e[ibs],Form("b%d_ipb_e[nch]/I",ibs)); // End index (exclusive) of merged phi bumps associated with this track
    t->Branch(Form("b%d_npb",ibs),&npb[ibs],Form("b%d_npb/I",ibs)); // number of merged phi bumps found for all tracks
    t->Branch(Form("b%d_pb_acc",ibs),&pb_acc[ibs],Form("b%d_pb_acc[b%d_npb]/I",ibs,ibs)); // whether this phi-bump falls within the max delta phi cut or not
    t->Branch(Form("b%d_pb_phi",ibs),&pb_phi[ibs],Form("b%d_pb_phi[b%d_npb]/F",ibs,ibs)); // phi agnle of of this merged phi bump
    t->Branch(Form("b%d_pb_the",ibs),&pb_the[ibs],Form("b%d_pb_the[b%d_npb]/F",ibs,ibs)); // theta agnle of of this merged phi bump (should be same as the cluster's phi)
    t->Branch(Form("b%d_pb_ene",ibs),&pb_ene[ibs],Form("b%d_pb_ene[b%d_npb]/F",ibs,ibs)); // energy of this phi bump (found binsong's splitting)
    t->Branch(Form("b%d_pb_rcalc",ibs),&pb_rcalc[ibs],Form("b%d_pb_rcalc[b%d_npb]/F",ibs,ibs)); // Recalculated radius of emission of this merged phi bump based on DeltaPhi and pT
    t->Branch(Form("b%d_pb_zcalc",ibs),&pb_zcalc[ibs],Form("b%d_pb_zcalc[b%d_npb]/F",ibs,ibs)); // Recalculated radius of emission of this merged phi bump based on DeltaPhi and pT
  }
  // Kind of impossible to figure out the best matching mc brem photon because digis are mangled when making phi-bumps
  //t->Branch("pb_match",&pb_match,"pb_match[npb]/I"); // The index of the best matching MC brem photon to this separated bump: no way to find this
  //t->Branch("pb_score",&pb_score,"pb_score[npb]/I"); // The number of common MC tracks with the matching MC track by def, this will be 100%

  // All bumps list
  t->Branch("nab",&nab,"nab/I");
  t->Branch("ab_phi",&ab_phi,"ab_phi[nab]/F"); // phi agnle of of this separated bump
  t->Branch("ab_the",&ab_the,"ab_the[nab]/F"); // theta agnle of of this separated bump
  t->Branch("ab_ene",&ab_ene,"ab_ene[nab]/F"); // reconstructed energy of this separated bump
  t->Branch("ab_isb",&ab_isb,"ab_isb[nab]/I"); // Index of the track for which this bump had been found to be a separted bump. -1 otherwise
  t->Branch("ab_ich",&ab_ich,"ab_ich[nab]/I"); // Index of the track with which this bump shares EmcIndex
  t->Branch("ab_match",&ab_match,"ab_match[nsb]/I"); // The index of the best matching MC brem photon to this bump
  t->Branch("ab_score",&ab_score,"ab_score[nsb]/I"); // The number of common MC tracks with the matching MC track

  return kSUCCESS;
}


void BremPidReader::Exec(Option_t* opt)
{

  int verb = 0;

  if (verb>0||(nEvt%1000==0))
    cout << "==== New Event " << nEvt << " ===== " << endl;

  //if (split_detector())
  //  print_cands();
  //++nEvt;
  //return;

  nChCand = fChargedCandidateArray->GetEntriesFast();
  nNeutCand = fNeutralCandidateArray->GetEntriesFast();
  nMcTrack = fMcArray->GetEntriesFast();

  nch = 0;
  nsb = 0;
  for (int ibs=0; ibs<nbs; ++ibs) npb[ibs] = 0;
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

    if (pdg==11||pdg==-11||pdg==13||pdg==-13) {
      if (mid!=-1&&mpdg!=433) continue; // keep only primary electrons. Its complicated as it is
    } else if (pdg==211||pdg==-211) {
      if (mid!=-1) continue;
    } else {
      continue;
    }

    is_prim[nch] = (mid==-1||(mpdg==433&&(pdg==1||pdg==-1)));
    mom_mc[nch] = truth->GetMomentum().Mag();
    phi_mc[nch] = truth->GetMomentum().Phi()*TMath::RadToDeg();
    the_mc[nch] = truth->GetMomentum().Theta()*TMath::RadToDeg();
    pdg_mc[nch] = truth->GetPdgCode();

    TVector3 mom = theChargedCand->GetMomentum();
    double ene = theChargedCand->GetEnergy();
    //TLorentzVector m4 = TLorentzVector(mom,ene);
    double mass = 0.000511; // Electron mass hypothesis (GeV)

    vector<vector<int> > sb_tree;

    double sepPhotonE = 0, sepPhotonE_wtd = 0, sepPhotonE_wtd_bf;
    GetSepPhotonE_fromBumps(theChargedCand, sepPhotonE, sepPhotonE_wtd, sepPhotonE_wtd_bf, sb_tree);
    // Different level of weighting only in separated
    mom_sep[nch] = corrected_mom(sepPhotonE);
    mom_sep_w[nch] = corrected_mom(sepPhotonE_wtd);
    mom_sep_w_bf[nch] = corrected_mom(sepPhotonE_wtd_bf);

    double mergPhotonE = 0, mergPhotonE_wtd = 0, mergPhotonE_wtd_bf = 0;
    double mergPhotonE_pc = 0, mergPhotonE_wtd_pc = 0, mergPhotonE_wtd_bf_pc = 0; // pc for phicut aka. political correctness
    //GetMergPhotonE(theChargedCand, i, mergPhotonE, mergPhotonE_wtd, mergPhotonE_wtd_bf, mergPhotonE_pc, mergPhotonE_wtd_pc, mergPhotonE_wtd_bf_pc);
    for (int ibs =0; ibs<nbs; ++ibs) {
      GetMergPhotonE(theChargedCand, ibs, mergPhotonE, mergPhotonE_wtd, mergPhotonE_wtd_bf, mergPhotonE_pc, mergPhotonE_wtd_pc, mergPhotonE_wtd_bf_pc);
      // Different level of weighting only in merged
      mom_mrg[ibs][nch] = corrected_mom(mergPhotonE);
      mom_mrg_w[ibs][nch] = corrected_mom(mergPhotonE_wtd);
      mom_mrg_w_bf[ibs][nch] = corrected_mom(mergPhotonE_wtd_bf);
      mom_mrg_pc[ibs][nch] = corrected_mom(mergPhotonE_pc);
      mom_mrg_w_pc[ibs][nch] = corrected_mom(mergPhotonE_wtd_pc);
      mom_mrg_w_bf_pc[ibs][nch] = corrected_mom(mergPhotonE_wtd_bf_pc);

      // Combine b/f weighted separated with different level in merged
      mom_cor[ibs][nch] = corrected_mom(sepPhotonE + mergPhotonE); // original correction
      mom_wcor[ibs][nch] = corrected_mom(sepPhotonE_wtd + mergPhotonE); // weight only separated
      mom_wbfcor[ibs][nch] = corrected_mom(sepPhotonE_wtd_bf + mergPhotonE); // weight separated b/f separately
      mom_wbfcor_mw_bf[ibs][nch] = corrected_mom(sepPhotonE_wtd_bf + mergPhotonE_wtd_bf); // weight both separated and merged b/f separately
      mom_wbfcor_mw_bf_pc[ibs][nch] = corrected_mom(sepPhotonE_wtd_bf + mergPhotonE_wtd_bf_pc); // weight both separated and merged b/f separately, phi cut merged
    }

    PndPidBremCorrected4Mom *tmp = (PndPidBremCorrected4Mom*)fBremCorrected4MomArray->At(iCand);
    if(tmp) {
      TVector3 storedCorrMom = tmp->GetMomentum();
      mom_stored[nch] = storedCorrMom.Mag();
    } else {
      mom_stored[nch] = -9999.;
    }

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
double BremPidReader::corrected_mom(const double& photTotEne ) {
  return ((photTotEne+fRecMomOfEle)/fRecMomOfEle) * mom_rec[nch];
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
      set_intersection(_sb_tree[isb].begin(),_sb_tree[isb].end(),
		       _mcb_tree[imcb].begin(),_mcb_tree[imcb].end(),
		       back_inserter(intersn));
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
void BremPidReader::mc_match_bumps_main_contrib(vector<PndEmcBump*> bump_list, vector<int>& mc_list) {
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

void BremPidReader::find_contributors(PndEmcBump *bump, vector<int>& tree) {
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

void BremPidReader::find_ancestory(const int& mcidx, vector<int>& tree) {
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

void BremPidReader::
get_mc_brem_photons(const int &mcidx, vector<vector<int> >& _mcb_tree) {
  imcb_s[nch] = nmcb;
  _nmcb[nch] = 0;
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
    _nmcb[nch]++;
  }
  imcb_e[nch] = nmcb;
  assert(imcb_e[nch] - imcb_s[nch] == _nmcb[nch]);
}

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

void BremPidReader::
GetSepPhotonE_fromBumps(PndPidCandidate *ChargedCand, double &esep, double &esep_wtd,
			double &esep_wtd_bf, vector<vector<int> > &_sb_tree) {

  esep = 0;
  esep_wtd = 0;
  esep_wtd_bf = 0;

  isb_s[nch] = nsb;
  _nsb[nch] = 0;

  const int iTrkEmcIdx = ChargedCand->GetEmcIndex();
  if (iTrkEmcIdx>=0) {

    Float_t PhotonTotEnergySep = 0;
    Float_t PhotonTotEnergySepWtd = 0;
    Float_t PhotonTotEnergySepWtdBf = 0;

    const int nBump = fBumpArray->GetEntriesFast();
    for(Int_t iBump = 0; iBump<nBump; ++iBump)
      {
	PndEmcBump *PhotonBump = (PndEmcBump *) fBumpArray->At(iBump);
	const Float_t PhotonEnergySep = PhotonBump->GetEnergyCorrected();

	const Int_t iSepClust = PhotonBump->GetClusterIndex();

	// Exclude bumps that have the same EmcIdx with *any* track
	const int ncand = fChargedCandidateArray->GetEntriesFast();
	bool is_asso = false;
	for (int icand=0; icand<ncand; ++icand) {
	  PndPidCandidate* tmpChargedCand = (PndPidCandidate*) fChargedCandidateArray->At(icand);
	  int tmpiTrkEmcIdx = tmpChargedCand->GetEmcIndex();
	  if ( iSepClust == tmpiTrkEmcIdx ) {
	    is_asso = true;
	    break;
	  }
	}
	if (is_asso) continue;

	PndEmcCluster *PhotonCluster = (PndEmcCluster*) fClusterArray->At(iSepClust);

	const Double_t PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
	const Double_t PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();

	const Bool_t fwd = fRecThetaOfEle <= 23.;
	const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
	const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
	const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;

	const Float_t RealDeltaPhi = fCharge<0?PhotonPhiSep-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiSep;
	const Float_t RealDeltaTheta = fCharge<0?PhotonThetaSep-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaSep;

	const Float_t rad_calc = 100*TMath::Sin(RealDeltaPhi*TMath::DegToRad()/2.)*2*Pt/0.3/2.0; // B=2T
	const Float_t zed_calc = rad_calc/TMath::Tan(TMath::DegToRad()*fRecThetaOfEle);

	const Float_t wt = 1.0/(1.+TMath::Exp((rad_calc-21.)/5.));
	const Float_t wt_bf = fwd ? 1.0/(1.+TMath::Exp((zed_calc-90.)/25.)) : 1.0/(1.+TMath::Exp((rad_calc-21.)/5.));
	const Float_t ThetaCutUp = 2.;
	const Float_t ThetaCutDown = -2.;
	const Float_t PhiCutUp = fwd ? DeltaPhiForward : DeltaPhiBarrel;
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
	  sb_zcalc[nsb] = zed_calc;
	  sb_ene[nsb] = PhotonEnergySep;
	  _nsb[nch]++;
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

	  esep += PhotonEnergySep;
	  esep_wtd += wt*PhotonEnergySep;
	  esep_wtd_bf += wt_bf*PhotonEnergySep;
	}

      }//loop neutralcand

  } // iTrkEmcIdx >= 0

  isb_e[nch] = nsb;
  assert(isb_e[nch] - isb_s[nch] == _nsb[nch]);
  nphot_sep[nch] = _nsb[nch];

  if (esep < fRecMomOfEle/100.) esep_wtd = 0;
  if (esep_wtd < fRecMomOfEle/100.) esep_wtd = 0;
  if (esep_wtd_bf < fRecMomOfEle/100.) esep_wtd_bf = 0;

}


void BremPidReader::GetMergPhotonE(PndPidCandidate *ChargedCand, int binSize,
				     double &emrg, double &emrg_wtd, double &emrg_wtd_bf,
				     double &emrg_pc, double &emrg_wtd_pc,  double &emrg_wtd_bf_pc){

  emrg = 0;
  emrg_wtd = 0;
  emrg_wtd_bf = 0;
  emrg_pc = 0;
  emrg_wtd_pc =0;
  emrg_wtd_bf_pc = 0;

  nphot_mrg[nch] = 0;
  nphot_mrg_pc[nch] = 0;

  ipb_s[binSize][nch] = npb[binSize];
  _npb[binSize][nch] = 0;

  // EMcal cluster associated with track ...
  if (ChargedCand->GetEmcIndex() >= 0) {

    PndEmcBump *EleBump = (PndEmcBump *) fBumpArray->At(ChargedCand->GetEmcIndex());
    Int_t EleRefCluster = EleBump->GetClusterIndex();

    if (EleRefCluster >= 0) {

      std::vector<PndEmcBump*> EmcPhiBumpList;
      int nPhiBump = fPhiBumpArray[binSize]->GetEntriesFast();
      for (int ipb=0; ipb<nPhiBump; ++ipb) {
	PndEmcBump *phibump = (PndEmcBump*) fPhiBumpArray[binSize]->At(ipb);
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

      Double_t PhotonTotEnergyMerg = 0.0;

      for(Int_t r = iS; r<=iE; r++) {
	const Double_t PhotonThetaMrg = EmcPhiBumpList[r]->position().Theta()*TMath::RadToDeg();
	const Double_t PhotonPhiMrg = EmcPhiBumpList[r]->position().Phi()*TMath::RadToDeg();
	const Bool_t fwd = fRecThetaOfEle <= 23.;
	const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
	const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
	const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;
	const Float_t RealDeltaPhi = fCharge<0?PhotonPhiMrg-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiMrg;
	const Float_t RealDeltaTheta = fCharge<0?PhotonThetaMrg-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaMrg;
	const Float_t rad_calc = 100*TMath::Sin(RealDeltaPhi*TMath::DegToRad()/2.)*2*Pt/0.3/2.0; // B=2T
	const Float_t zed_calc = rad_calc/TMath::Tan(TMath::DegToRad()*fRecThetaOfEle);
	const Float_t wt = 1.0/(1.+TMath::Exp((rad_calc-21.)/5.));
	const Float_t wt_bf = fwd ? 1.0/(1.+TMath::Exp((zed_calc-90.)/25.)) : 1.0/(1.+TMath::Exp((rad_calc-21.)/5.));
	const Float_t ThetaCutUp = 2.;
	const Float_t ThetaCutDown = -2.;
	const Float_t PhiCutUp = fwd ? DeltaPhiForward : DeltaPhiBarrel;
	const Float_t PhiCutDown = -1;
	const Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
	const Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;

	pb_phi[binSize][npb[binSize]] = PhotonPhiMrg;
	pb_the[binSize][npb[binSize]] = PhotonThetaMrg;
	pb_rcalc[binSize][npb[binSize]] = rad_calc;
	pb_zcalc[binSize][npb[binSize]] = zed_calc;
	pb_ene[binSize][npb[binSize]] = EmcPhiBumpList[r]->energy();
	pb_acc[binSize][npb[binSize]] = 0;

	emrg += EmcPhiBumpList[r]->energy();
	emrg_wtd = wt * EmcPhiBumpList[r]->energy();
	emrg_wtd_bf = wt_bf * EmcPhiBumpList[r]->energy();
	nphot_mrg[nch]++;
	if (PhiCut&&ThetaCut) {
	  pb_acc[binSize][npb[binSize]] = 1;
	  emrg_pc += EmcPhiBumpList[r]->energy();
	  emrg_wtd_pc = wt * EmcPhiBumpList[r]->energy();
	  emrg_wtd_bf_pc = wt_bf * EmcPhiBumpList[r]->energy();
	  nphot_mrg_pc[nch]++;
	}
	_npb[binSize][nch]++;
	npb[binSize]++;
      }

      if (emrg < fRecMomOfEle/100.) emrg = 0;
      if (emrg_wtd < fRecMomOfEle/100.) emrg_wtd = 0;
      if (emrg_wtd_bf < fRecMomOfEle/100.) emrg_wtd_bf = 0;
      if (emrg_pc < fRecMomOfEle/100.) emrg_pc = 0;
      if (emrg_wtd_pc < fRecMomOfEle/100.) emrg_wtd_pc = 0;
      if (emrg_wtd_bf_pc < fRecMomOfEle/100.) emrg_wtd_bf_pc = 0;

    } // EleRefCluster >= 0

  } // ChargedCand->GetEmcIndex() >= 0

  ipb_e[binSize][nch] = npb[binSize];
  assert(ipb_e[binSize][nch] - ipb_s[binSize][nch] == _npb[binSize][nch]);

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

bool BremPidReader::split_detector() {
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

    int mcsize = bump->GetMcSize();
    cout << " McSize= " << mcsize << " McList: (";
    for (int i=0; i<mcsize; ++i) {
      cout << bump->GetMcIndex(i) << (i<(mcsize-1)?", ": "");
    }
    cout << ")";

    int digisize = bump->DigiList().size();
    cout << " NDigi= " << digisize << "(";
    for (int i=0; i<digisize; ++i) {
      //cout << ((PndEmcDigi*)fDigiArray->At(clust->DigiList()[i]))->GetTrackId() << (i<(digisize-1)?", ": "");
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

//double BremPidReader::GetSepPhotonE(PndPidCandidate *ChargedCand, int &nphotsep){
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
//      const double PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
//      const double PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();
//
//      const Float_t Pt = fRecMomOfEle*TMath::Sin(fRecThetaOfEle/TMath::RadToDeg());
//      const Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
//      const Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(fRecThetaOfEle/57.3)*57.3;
//      //cout << "DphBar= " << DeltaPhiBarrel << " DphFor= " << DeltaPhiForward << endl;
//
//      const Float_t RealDeltaPhi = fCharge<0?PhotonPhiSep-fRecPhiOfEle:fRecPhiOfEle-PhotonPhiSep;
//      const Float_t RealDeltaTheta = fCharge<0?PhotonThetaSep-fRecThetaOfEle:fRecThetaOfEle-PhotonThetaSep;
//
//      const Float_t ThetaCutUp = 2.;
//      const Float_t ThetaCutDown = -2.;
//      const Float_t PhiCutUp = (fRecThetaOfEle <= 23.)?DeltaPhiForward:DeltaPhiBarrel;
//      const Float_t PhiCutDown = -1;
//
//      const Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
//      const Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
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

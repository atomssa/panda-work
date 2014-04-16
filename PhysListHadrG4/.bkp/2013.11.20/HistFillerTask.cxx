// The header file
#include "HistFillerTask.h"

// C++ headers
#include <string>
#include <iostream>

// FAIR headers
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

// ROOT headers
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

// other headers
#include "PndEmcBump.h"
#include "PndEmcCluster.h"
#include "PndEmcBump.h"
#include "PndEmcRecoHit.h"
#include "PndEmcMapper.h"
#include "PndMCTrack.h"

using std::cout;
using std::endl;

HistFillerTask::HistFillerTask() :
  FairTask("Simple Histo Filler Task"),nevt(0) { 
}

HistFillerTask::~HistFillerTask() { }


InitStatus HistFillerTask::Init() 
{		

  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- HistFillerTask::Init: " << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  track_array = dynamic_cast<TClonesArray *> (ioman->GetObject("MCTrack"));
  if ( ! track_array ) {
    cout << "-W- HistFillerTask::Init: " << "No PndMCTrack array!" << endl;
    return kERROR;
  }

  clust_array = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcCluster"));
  if ( ! clust_array ) {
    cout << "-W- HistFillerTask::Init: " << "No EmcCluster array!" << endl;
    return kERROR;
  }

  bump_array = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcBump"));
  if ( ! bump_array ) {
    cout << "-W- HistFillerTask::Init: " << "No PndEmcBump array!" << endl;
    return kERROR;
  }

  reco_hit_array = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcRecoHit"));
  if ( ! bump_array ) {
    cout << "-W- HistFillerTask::Init: " << "No EmcRecoHit array!" << endl;
    return kERROR;
  }

  initialize_hists();
  
  return kSUCCESS;
}

void HistFillerTask::initialize_hists() {

  h_energy_b = new TH1F("h_energy_b","Reconstructed Energy",200,0,1.1);
  h_energy_corr_b = new TH1F("h_energy_corr_b","Reconstructed Corrected Energy",200,0,1.1);

  h_energy = new TH1F("h_energy","Reconstructed Energy (bef. Acc cut)",200,0,1.1);
  h_energy_corr = new TH1F("h_energy_corr","Reconstructed Corrected Energy (bef. Acc cut)",200,0,1.1);

  h_phi_the_b = new TH2F("h_phi_the_b","Track #phi vs #theta before Acc. cut",200,-180,180,200,0,180);
  h_phi_the = new TH2F("h_phi_the","Track #phi vs #theta before Acc. cut",200,-180,180,200,0,180);
  
  h_dphi = new TH1F("d_dphi","",200,-25,25);
  h_dthe = new TH1F("d_dthe","",200,-10,10);
  
}

void HistFillerTask::Exec(Option_t* opt)
{

  if (nevt%1000==0)
    cout << "===== HistFillerTask::Exec -- Event " << nevt << " ====="<< endl;
  PndEmcMapper::Init(1);

  if ( ! track_array ) Fatal("Exec", "No Track Array");
  if ( ! clust_array ) Fatal("Exec", "No Cluster Array");
  if ( ! bump_array ) Fatal("Exec", "No Bump Array");  
  if ( ! reco_hit_array ) Fatal("Exec", "No RecoHit Array");  

  int debug = 0;
  
  const int ntrk = track_array->GetEntriesFast();
  //if (debug>1) cout << "ntrk= " << ntrk << endl;
  const PndMCTrack* track = (PndMCTrack*) track_array->At(0);
  const TVector3 trk_mom = track->GetMomentum();
  const TLorentzVector trk_mom4 = track->Get4Momentum();
  const double trk_phi = trk_mom.Phi();
  const double trk_the = trk_mom.Theta();
  const double trk_the_deg = trk_the * 180. / TMath::Pi();
  const double trk_phi_deg = trk_phi * 180. / TMath::Pi();
  const double trk_mom_mag = trk_mom.Mag();
  const TVector3 trk_vtx = track->GetStartVertex();
  const double trk_vtxr = TMath::Hypot(trk_vtx(0),trk_vtx(1));
  const double trk_vtxz = trk_vtx(2);
  const double trk_energy = trk_mom4.E();
    
  const int pdg = track->GetPdgCode();
  const int motherid = track->GetMotherID();
  const int grammaid = track->GetSecondMotherID();
      
  if (debug>1) cout << "Track: pdg= " << pdg << " pID= " << motherid
		    << " gpID= " << grammaid << " vR= " << trk_vtxr << " vZ= " << trk_vtxz
		    << " phi= " << trk_phi << " the= " << trk_the << " |mom|= " << trk_mom_mag << " E= " << trk_energy << endl;

  h_phi_the_b->Fill(trk_phi_deg,trk_the_deg);
  
  // acceptance cut for beam tubes
  const bool inacc = !((trk_phi_deg>-100&&trk_phi_deg<-90) || (trk_phi_deg>90&&trk_phi_deg<100) );

  if (inacc)
    h_phi_the->Fill(trk_phi_deg,trk_the_deg);
  
  //-------- clusters -----------//
  const int nclust = clust_array->GetEntriesFast();
  if (debug>2) cout << "nclust= " << nclust << endl;
  double del_the_min = 1e10;
  PndEmcCluster* clust;
  double energy_tot_clust = 0.0;
  if (debug>3) cout << "Nclust= " << nclust;
  for (int iclust=0; iclust<nclust; ++iclust) {
    PndEmcCluster *tmp_clust = (PndEmcCluster*)clust_array->At(iclust);
    energy_tot_clust += tmp_clust->energy();
    if (debug>3) cout << " [E_" << iclust << "= " << tmp_clust->energy() << ", Th_" << iclust << "= " << tmp_clust->where().Theta() << "]   ";
    if ( del_the_min > fabs(tmp_clust->where().Theta()-trk_the) ) {
      clust = tmp_clust;
      del_the_min = fabs(tmp_clust->where().Theta()-trk_the);
    }
  }
  if (debug>3) cout << endl;
  if (nclust > 0 ) { 
    const TVector3 clust_pos = clust->where();
    const double clust_energy = clust->energy();
    const double clust_phi = clust_pos.Phi();
    const double clust_the = clust_pos.Theta();
    if (debug>2) cout << "Clust: E_clust = " << clust_energy << " phi= " << clust_phi << " the= " << clust_the << " E_tot= " << energy_tot_clust <<endl;
  }
  
  //---------- bump ---------//
  const int nbump = bump_array->GetEntriesFast();
  if (debug>2) cout << "nbump= " << nbump << endl;
  del_the_min = 1e10;
  PndEmcBump* bump;
  double energy_tot_bump = 0.0;
  if (debug>3) cout << "Nbump= " << nbump;
  for (int ibump=0; ibump<nbump; ++ibump) {
    PndEmcBump *tmp_bump = (PndEmcBump*)bump_array->At(ibump);
    energy_tot_bump += tmp_bump->energy();
    if (debug>3) cout << " [E_" << ibump << "= " << tmp_bump->energy() << ", Th_" << ibump << "= " << tmp_bump->where().Theta() << "]   ";
    if ( del_the_min > fabs(tmp_bump->where().Theta()-trk_the) ) {
      bump = tmp_bump;
      del_the_min = fabs(tmp_bump->where().Theta()-trk_the);
    }
  }
  if (debug>3) cout << endl;
  if (nbump>0) {
    const TVector3 bump_pos = bump->where();
    const double bump_energy = bump->energy();
    const double bump_phi = bump_pos.Phi();
    const double bump_the = bump_pos.Theta();
    if (debug>2) cout << "Bump: E_bump = " << bump_energy << " phi= " << bump_phi << " the= " << bump_the << " E_tot= " << energy_tot_bump <<endl;
  }


  
  //---------- reco_hit ---------//
  const int nreco_hit = reco_hit_array->GetEntriesFast();
  if (debug>2) cout << "nreco_hit= " << nreco_hit << endl;
  del_the_min = 1e10;
  PndEmcRecoHit* reco_hit;
  double energy_tot_reco_hit = 0.0;
  double energy_tot_reco_hit_corr = 0.0;
  if (debug>3) cout << "Nreco_hit= " << nreco_hit;
  for (int ireco_hit=0; ireco_hit<nreco_hit; ++ireco_hit) {
    PndEmcRecoHit *tmp_reco_hit = (PndEmcRecoHit*)reco_hit_array->At(ireco_hit);
    energy_tot_reco_hit += tmp_reco_hit->GetEnergy();
    energy_tot_reco_hit_corr += tmp_reco_hit->GetEnergyCorrected();
    if (debug>3) cout << " [E_" << ireco_hit << "= " << tmp_reco_hit->GetEnergy() << ", Th_" << ireco_hit << "= " << tmp_reco_hit->GetPosition().Theta() << "]   ";
    if ( del_the_min > fabs(tmp_reco_hit->GetPosition().Theta()-trk_the) ) {
      reco_hit = tmp_reco_hit;
      del_the_min = fabs(tmp_reco_hit->GetPosition().Theta()-trk_the);
    }
  }
  if (debug>3) cout << endl;
  if (nreco_hit>0) {
    const TVector3 reco_hit_pos = reco_hit->GetPosition();
    const double reco_hit_energy = reco_hit->GetEnergy();
    const double reco_hit_energy_corr = reco_hit->GetEnergyCorrected();
    const double reco_hit_phi = reco_hit_pos.Phi();
    const double reco_hit_the = reco_hit_pos.Theta();
    if (debug>2) cout << "Reco_Hit: E_reco_hit = " << reco_hit_energy << " E_reco_hit_corr = " << reco_hit_energy_corr
		      << " phi= " << reco_hit_phi << " the= " << reco_hit_the << " E_tot= " << energy_tot_reco_hit
		      << " E_tot_corr= " << energy_tot_reco_hit_corr <<endl;
    
    const double reco_hit_phi_deg = reco_hit_phi*180./TMath::Pi();
    const double reco_hit_the_deg = reco_hit_the*180./TMath::Pi();
    if (inacc) {
      h_energy->Fill(reco_hit_energy/trk_energy);
      h_energy_corr->Fill(reco_hit_energy_corr/trk_energy);
      h_dphi->Fill(trk_phi_deg-reco_hit_phi_deg);
      h_dthe->Fill(trk_the_deg-reco_hit_the_deg);
    }
    h_energy_b->Fill(reco_hit_energy/trk_energy);
    h_energy_corr_b->Fill(reco_hit_energy_corr/trk_energy);
    
  }

  

  
  ++nevt;

}


void HistFillerTask::Finish()
{	

  h_energy->Write();
  h_energy_corr->Write();

  h_energy_b->Write();
  h_energy_corr_b->Write();

  h_phi_the->Write();
  h_phi_the_b->Write();

  h_dphi->Write();
  h_dthe->Write();
  
}

ClassImp(HistFillerTask)

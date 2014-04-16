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

  // Reconstructed cluster energy distribution before acceptance cut
  h_energy_b = new TH1F("h_energy_b","Reconstructed Energy",200,0,1.1);
  h_energy_corr_b = new TH1F("h_energy_corr_b","Reconstructed Corrected Energy",200,0,1.1);

  // Reconstructed cluster energy distribution after acceptance cut
  h_energy = new TH1F("h_energy","Reconstructed Energy (bef. Acc cut)",200,0,1.1);
  h_energy_corr = new TH1F("h_energy_corr","Reconstructed Corrected Energy (bef. Acc cut)",200,0,1.1);

  // angluar position of cluster before & after acceptance cuts
  h_phi_the_b = new TH2F("h_phi_the_b","Track #phi vs #theta before Acc. cut",200,-180,180,200,0,180);
  h_phi_the = new TH2F("h_phi_the","Track #phi vs #theta before Acc. cut",200,-180,180,200,0,180);

  // Residual of angular position between cluster and track initial theta angle before & after acceptance cuts
  h_dphi = new TH1F("d_dphi","",200,-25,25);
  h_dthe = new TH1F("d_dthe","",200,-10,10);

  
  // Cluster multiplicity
  h_nclust = new TH1F("h_nclust","Number of clusters",20,0,20);
  h_nbump = new TH1F("h_nbump","Number of bumps",20,0,20);
  h_nreco_hit = new TH1F("h_nreco_hit","Number of reco hits",20,0,20);

  // Angluar distance between closest cluster to track and other clusters
  h_extra_clust_dphi_dthe = new TH2F("h_extra_clust_dphi_dthe","Distance of extra clusturs from main cluster",200,-180,180,200,0,180);
  h_extra_bump_dphi_dthe = new TH2F("h_extra_bump_dphi_dthe","Distance of extra bumps from main cluster",200,-180,180,200,0,180);
  h_extra_reco_hit_dphi_dthe = new TH2F("h_extra_reco_hit_dphi_dthe","Distance of extra bumps from main cluster",200,-180,180,200,0,180);

  // Energy distribution of extra clusters
  h_extra_clust_energy = new TH1F("h_extra_clust_energy","Reconstructed energy of extra clusters",200,0,1.1);
  h_extra_bump_energy = new TH1F("h_extra_bump_energy","Reconstructed energy of extra bumps",200,0,1.1);
  h_extra_reco_hit_energy = new TH1F("h_extra_reco_hit_energy","Reconstructed energy of extra reco hits",200,0,1.1);

  // Energy distribution of extra clusters
  h_all_clust_energy = new TH1F("h_all_clust_energy","Reconstructed energy of all clusters",200,0,1.1);
  h_all_bump_energy = new TH1F("h_all_bump_energy","Reconstructed energy of all bumps",200,0,1.1);
  h_all_reco_hit_energy = new TH1F("h_all_reco_hit_energy","Reconstructed energy of all reco hits",200,0,1.1);

  // Total energy (sum) of all clusters
  h_clust_tot_energy = new TH1F("h_clust_tot_energy","Sum total of reconstructed energy of all clusters",200,0,1.1);
  h_bump_tot_energy = new TH1F("h_bump_tot_energy","Sum total of reconstructed energy of all bumps",200,0,1.1);
  h_reco_hit_tot_energy = new TH1F("h_reco_hit_tot_energy","Sum total of reconstructed energy of all reco hits",200,0,1.1);
  h_reco_hit_tot_energy_corr = new TH1F("h_reco_hit_tot_energy_corr","Sum total of reconstructed corrected energy of all reco hits",200,0,1.1);
  
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
  h_nclust->Fill(nclust);
  if (debug>2) cout << "nclust= " << nclust << endl;
  double del_the_min = 1e10;
  int i_main_clust=0;
  PndEmcCluster* main_clust;
  double energy_tot_clust = 0.0;
  if (debug>3) cout << "Nclust= " << nclust;
  for (int iclust=0; iclust<nclust; ++iclust) {
    PndEmcCluster *tmp_clust = (PndEmcCluster*)clust_array->At(iclust);
    energy_tot_clust += tmp_clust->energy();
    h_all_clust_energy->Fill(tmp_clust->energy());
    if (debug>3) cout << " [E_" << iclust << "= " << tmp_clust->energy() << ", Th_" << iclust << "= " << tmp_clust->where().Theta() << "]   ";
    if ( del_the_min > fabs(tmp_clust->where().Theta()-trk_the) ) {
      main_clust = tmp_clust;
      i_main_clust = iclust;
      del_the_min = fabs(tmp_clust->where().Theta()-trk_the);
    }
  }
  h_clust_tot_energy->Fill(energy_tot_clust);    
  if (debug>3) cout << endl;
  if (nclust > 0 ) { 
    const TVector3 clust_pos = main_clust->where();
    const double clust_energy = main_clust->energy();
    const double clust_phi = clust_pos.Phi();
    const double clust_the = clust_pos.Theta();
    for (int iclust=0; iclust<nclust; ++iclust) {
      if (iclust==i_main_clust) continue;
      PndEmcCluster *tmp_clust = (PndEmcCluster*)clust_array->At(iclust);
      const double _dphi = tmp_clust->where().Phi() - main_clust->where().Phi();
      const double _dthe = tmp_clust->where().Theta() - main_clust->where().Theta();
      h_extra_clust_dphi_dthe->Fill(_dphi,_dthe);
      h_extra_clust_energy->Fill(clust_energy);
    }
    if (debug>2) cout << "Clust: E_clust = " << clust_energy << " phi= " << clust_phi << " the= " << clust_the << " E_tot= " << energy_tot_clust <<endl;
  }





  //-------- bumps -----------//
  const int nbump = bump_array->GetEntriesFast();
  h_nbump->Fill(nbump);
  if (debug>2) cout << "nbump= " << nbump << endl;
  double del_the_min_bump = 1e10;
  int i_main_bump=0;
  PndEmcBump* main_bump;
  double energy_tot_bump = 0.0;
  if (debug>3) cout << "Nbump= " << nbump;
  for (int ibump=0; ibump<nbump; ++ibump) {
    PndEmcBump *tmp_bump = (PndEmcBump*)bump_array->At(ibump);
    energy_tot_bump += tmp_bump->energy();
    h_all_bump_energy->Fill(tmp_bump->energy());
    if (debug>3) cout << " [E_" << ibump << "= " << tmp_bump->energy() << ", Th_" << ibump << "= " << tmp_bump->where().Theta() << "]   ";
    if ( del_the_min_bump > fabs(tmp_bump->where().Theta()-trk_the) ) {
      main_bump = tmp_bump;
      i_main_bump = ibump;
      del_the_min_bump = fabs(tmp_bump->where().Theta()-trk_the);
    }
  }
  h_bump_tot_energy->Fill(energy_tot_bump);    
  if (debug>3) cout << endl;
  if (nbump > 0 ) { 
    const TVector3 bump_pos = main_bump->where();
    const double bump_energy = main_bump->energy();
    const double bump_phi = bump_pos.Phi();
    const double bump_the = bump_pos.Theta();
    for (int ibump=0; ibump<nbump; ++ibump) {
      if (ibump==i_main_bump) continue;
      PndEmcBump *tmp_bump = (PndEmcBump*)bump_array->At(ibump);
      const double _dphi = tmp_bump->where().Phi() - main_bump->where().Phi();
      const double _dthe = tmp_bump->where().Theta() - main_bump->where().Theta();
      h_extra_bump_dphi_dthe->Fill(_dphi,_dthe);
      h_extra_bump_energy->Fill(bump_energy);
    }
    if (debug>2) cout << "Bump: E_bump = " << bump_energy << " phi= " << bump_phi << " the= " << bump_the << " E_tot= " << energy_tot_bump <<endl;
  }

  





  ////---------- bump ---------//
  //const int nbump = bump_array->GetEntriesFast();
  //h_nbump->Fill(nbump);
  //if (debug>2) cout << "nbump= " << nbump << endl;
  //del_the_min = 1e10;
  //PndEmcBump* bump;
  //double energy_tot_bump = 0.0;
  //if (debug>3) cout << "Nbump= " << nbump;
  //for (int ibump=0; ibump<nbump; ++ibump) {
  //  PndEmcBump *tmp_bump = (PndEmcBump*)bump_array->At(ibump);
  //  energy_tot_bump += tmp_bump->energy();
  //  if (debug>3) cout << " [E_" << ibump << "= " << tmp_bump->energy() << ", Th_" << ibump << "= " << tmp_bump->where().Theta() << "]   ";
  //  if ( del_the_min > fabs(tmp_bump->where().Theta()-trk_the) ) {
  //    bump = tmp_bump;
  //    del_the_min = fabs(tmp_bump->where().Theta()-trk_the);
  //  }
  //}
  //if (debug>3) cout << endl;
  //if (nbump>0) {
  //  const TVector3 bump_pos = bump->where();
  //  const double bump_energy = bump->energy();
  //  const double bump_phi = bump_pos.Phi();
  //  const double bump_the = bump_pos.Theta();
  //  if (debug>2) cout << "Bump: E_bump = " << bump_energy << " phi= " << bump_phi << " the= " << bump_the << " E_tot= " << energy_tot_bump <<reco;
  //}











  //-------- endl_hit -----------//
  const int nreco_hit = reco_hit_array->GetEntriesFast();
  h_nreco_hit->Fill(nreco_hit);
  if (debug>2) cout << "nreco_hit= " << nreco_hit << endl;
  double del_the_reco_hit = 1e10;
  int i_main_reco_hit=0;
  PndEmcRecoHit* main_reco_hit;
  double energy_tot_reco_hit = 0.0;
  double energy_tot_reco_hit_corr = 0.0;
  if (debug>3) cout << "Nreco_hit= " << nreco_hit;
  for (int ireco_hit=0; ireco_hit<nreco_hit; ++ireco_hit) {
    PndEmcRecoHit *tmp_reco_hit = (PndEmcRecoHit*)reco_hit_array->At(ireco_hit);
    energy_tot_reco_hit += tmp_reco_hit->GetEnergy();
    h_all_reco_hit_energy->Fill(tmp_reco_hit->GetEnergyCorrected());
    if (debug>3) cout << " [E_" << ireco_hit << "= " << tmp_reco_hit->GetEnergy() << ", Th_" << ireco_hit << "= " << tmp_reco_hit->GetPosition().Theta() << "]   ";
    if ( del_the_reco_hit > fabs(tmp_reco_hit->GetPosition().Theta()-trk_the) ) {
      main_reco_hit = tmp_reco_hit;
      i_main_reco_hit = ireco_hit;
      del_the_reco_hit = fabs(tmp_reco_hit->GetPosition().Theta()-trk_the);
    }
  }
  h_reco_hit_tot_energy->Fill(energy_tot_reco_hit);
  h_reco_hit_tot_energy_corr->Fill(energy_tot_reco_hit_corr);    
  if (debug>3) cout << endl;
  if (nreco_hit > 0 ) { 
    const TVector3 reco_hit_pos = main_reco_hit->GetPosition();
    const double reco_hit_energy = main_reco_hit->GetEnergy();
    const double reco_hit_energy_corr = main_reco_hit->GetEnergyCorrected();
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
    
    for (int ireco_hit=0; ireco_hit<nreco_hit; ++ireco_hit) {
      if (ireco_hit==i_main_reco_hit) continue;
      PndEmcRecoHit *tmp_reco_hit = (PndEmcRecoHit*)reco_hit_array->At(ireco_hit);
      const double _dphi = tmp_reco_hit->GetPosition().Phi() - main_reco_hit->GetPosition().Phi();
      const double _dthe = tmp_reco_hit->GetPosition().Theta() - main_reco_hit->GetPosition().Theta();
      h_extra_reco_hit_dphi_dthe->Fill(_dphi,_dthe);
      h_extra_reco_hit_energy->Fill(reco_hit_energy);
    }
    if (debug>2) cout << "Reco_Hit: E_reco_hit = " << reco_hit_energy << " phi= " << reco_hit_phi << " the= " << reco_hit_the << " E_tot= " << energy_tot_reco_hit <<endl;
  }

  
  
  ////---------- reco_hit ---------//
  //const int nreco_hit = reco_hit_array->GetEntriesFast();
  //h_nreco_hit->Fill(nreco_hit);
  //if (debug>2) cout << "nreco_hit= " << nreco_hit << endl;
  //del_the_min = 1e10;
  //PndEmcRecoHit* reco_hit;
  //double energy_tot_reco_hit = 0.0;
  //double energy_tot_reco_hit_corr = 0.0;
  //if (debug>3) cout << "Nreco_hit= " << nreco_hit;
  //for (int ireco_hit=0; ireco_hit<nreco_hit; ++ireco_hit) {
  //  PndEmcRecoHit *tmp_reco_hit = (PndEmcRecoHit*)reco_hit_array->At(ireco_hit);
  //  energy_tot_reco_hit += tmp_reco_hit->GetEnergy();
  //  energy_tot_reco_hit_corr += tmp_reco_hit->GetEnergyCorrected();
  //  if (debug>3) cout << " [E_" << ireco_hit << "= " << tmp_reco_hit->GetEnergy() << ", Th_" << ireco_hit << "= " << tmp_reco_hit->GetPosition().Theta() << "]   ";
  //  if ( del_the_min > fabs(tmp_reco_hit->GetPosition().Theta()-trk_the) ) {
  //    reco_hit = tmp_reco_hit;
  //    del_the_min = fabs(tmp_reco_hit->GetPosition().Theta()-trk_the);
  //  }
  //}
  //if (debug>3) cout << endl;
  //if (nreco_hit>0) {
  //  const TVector3 reco_hit_pos = reco_hit->GetPosition();
  //  const double reco_hit_energy = reco_hit->GetEnergy();
  //  const double reco_hit_energy_corr = reco_hit->GetEnergyCorrected();
  //  const double reco_hit_phi = reco_hit_pos.Phi();
  //  const double reco_hit_the = reco_hit_pos.Theta();
  //  if (debug>2) cout << "Reco_Hit: E_reco_hit = " << reco_hit_energy << " E_reco_hit_corr = " << reco_hit_energy_corr
  //		      << " phi= " << reco_hit_phi << " the= " << reco_hit_the << " E_tot= " << energy_tot_reco_hit
  //		      << " E_tot_corr= " << energy_tot_reco_hit_corr <<endl;
  //  
  //  const double reco_hit_phi_deg = reco_hit_phi*180./TMath::Pi();
  //  const double reco_hit_the_deg = reco_hit_the*180./TMath::Pi();
  //  if (inacc) {
  //    h_energy->Fill(reco_hit_energy/trk_energy);
  //    h_energy_corr->Fill(reco_hit_energy_corr/trk_energy);
  //    h_dphi->Fill(trk_phi_deg-reco_hit_phi_deg);
  //    h_dthe->Fill(trk_the_deg-reco_hit_the_deg);
  //  }
  //  h_energy_b->Fill(reco_hit_energy/trk_energy);
  //  h_energy_corr_b->Fill(reco_hit_energy_corr/trk_energy);
  //}

  

  
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

// The header file
#include "HistFillerTask.h"

// C++ headers
#include <string>
#include <iostream>
#include <algorithm>

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
#include "PndEmcHit.h"
#include "PndEmcDigi.h"
#include "PndEmcBump.h"
#include "PndEmcCluster.h"
#include "PndEmcBump.h"
#include "PndEmcRecoHit.h"
#include "PndEmcMapper.h"
#include "PndMCTrack.h"

using std::cout;
using std::endl;

HistFillerTask::HistFillerTask(double rmin, double rmax) :
  FairTask("Simple Histo Filler Task"),nevt(0),ratio_min(rmin),ratio_max(rmax) {
  debug=0;
  initialize_remainder();
  initialize_child();

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

  digi_array = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcDigi"));
  if ( ! digi_array ) {
    cout << "-W- HistFillerTask::Init: " << "No PndEmcDigi array!" << endl;
    return kERROR;
  }

  /*
  hit_array = dynamic_cast<TClonesArray *> (ioman->GetObject("EmcHit"));
  if ( ! hit_array ) {
    cout << "-W- HistFillerTask::Init: " << "No PndEmcHit array!" << endl;
    return kERROR;
  }
  */
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

  // classified energy distributions
  for (int i=0 ;i<ncl; ++i) {
    h_energy_class[i] = new TH1F(Form("h_energy_class_%d",i),Form("h_energy_class_%d",i),200,0.,1.1);
    h_energy_class_nopi[i] = new TH1F(Form("h_energy_class_nopi_%d",i),Form("h_energy_class_nopi_%d",i),200,0.,1.1);
    h_eratio_vs_r[i] = new TH2F(Form("h_eratio_vs_r_class_%d",i),Form("h_raito_vs_r_class_%d",i),200,0.,1.1,200,50,80);
    h_eratio_vs_e_child_pip[i] = new TH2F(Form("h_eratio_vs_e_child_pip_class_%d",i),Form("h_raito_vs_e_child_pip_class_%d",i),200,0.,1.1,200,0.0,0.6);
    h_eratio_vs_mom_child_pip[i] = new TH2F(Form("h_eratio_vs_mom_child_pip_class_%d",i),Form("h_raito_vs_mom_child_pip_class_%d",i),200,0.,1.1,200,0.0,0.6);    
  }

  for (int ing=0; ing<ngb; ++ing) {
    for (int inp=0; inp<npb; ++inp) {
      for (int inn=0; inn<nnb; ++inn) {
	const char *name = Form("h_%dg_%dp_%dn",ing,inp,inn);
	const char *title = Form("h_%dg_%dp_%dn",ing,inp,inn);
	h_energy_class_4d[ing][inp][inn] = new TH1F(name,title,200,0,1.1);
      }
    }
  }

  h_c_ng = new TH1F("h_c_ng","h_c_ng",12,0,12);
  h_c_nep = new TH1F("h_c_nep","h_c_nep",4,0,4);
  h_c_nem = new TH1F("h_c_nem","h_c_nem",8,0,8);
  h_c_np = new TH1F("h_c_np","h_c_np",9,0,9);
  h_c_nn = new TH1F("h_c_nn","h_c_nn",25,0,25);
  h_c_npip = new TH1F("h_c_npip","h_c_npip",5,0,5);
  h_c_npim = new TH1F("h_c_npim","h_c_npim",5,0,5);
  h_c_npi0 = new TH1F("h_c_npi0","h_c_npi0",5,0,5);
  h_c_nmum = new TH1F("h_c_nmum","h_c_nmum",4,0,4);
  h_c_nmup = new TH1F("h_c_nmup","h_c_nmup",4,0,4);
  h_c_nnumu = new TH1F("h_c_nnumu","h_c_nnumu",4,0,4);
  h_c_nanumu = new TH1F("h_c_nanumu","h_c_nanumu",4,0,4);
  
  h_r_ng = new TH1F("h_r_ng","h_r_ng",70,0,70);
  h_r_nep = new TH1F("h_r_nep","h_r_nep",16,0,16);
  h_r_nem = new TH1F("h_r_nem","h_r_nem",20,0,20);
  h_r_np = new TH1F("h_r_np","h_r_np",7,0,7);
  h_r_nn = new TH1F("h_r_nn","h_r_nn",30,0,30);
  h_r_npip = new TH1F("h_r_npip","h_r_npip",6,0,6);
  h_r_npim = new TH1F("h_r_npim","h_r_npim",5,0,5);
  h_r_npi0 = new TH1F("h_r_npi0","h_r_npi0",5,0,5);
  h_r_nmum = new TH1F("h_r_nmum","h_r_nmum",4,0,4);
  h_r_nmup = new TH1F("h_r_nmup","h_r_nmup",4,0,4);
  h_r_nnumu = new TH1F("h_r_nnumu","h_r_nnumu",4,0,4);
  h_r_nanumu = new TH1F("h_r_nanumu","h_r_nanumu",4,0,4);

}

void HistFillerTask::Exec(Option_t* opt)
{

  if (nevt%1000==0) cout << "===== HistFillerTask::Exec -- Event " << nevt << " ====="<< endl;
  

  PndEmcMapper::Init(1);

  if ( ! track_array ) Fatal("Exec", "No Track Array");
  if ( ! clust_array ) Fatal("Exec", "No Cluster Array");
  if ( ! bump_array ) Fatal("Exec", "No Bump Array");  
  if ( ! reco_hit_array ) Fatal("Exec", "No RecoHit Array");  

  const int ntrk = track_array->GetEntriesFast();
  if (debug>1)
    cout << "ntrk= " << ntrk << endl;
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
      
  if (debug>1)
    cout << "Track: pdg= " << pdg << " pID= " << motherid
		    << " gpID= " << grammaid << " vR= " << trk_vtxr << " vZ= " << trk_vtxz
		    << " phi= " << trk_phi << " the= " << trk_the << " |mom|= " << trk_mom_mag << " E= " << trk_energy << endl;

  h_phi_the_b->Fill(trk_phi_deg,trk_the_deg);
  
  // acceptance cut for beam tubes
  const bool inacc = !((trk_phi_deg>-100&&trk_phi_deg<-90) || (trk_phi_deg>90&&trk_phi_deg<100) );

  if (inacc) h_phi_the->Fill(trk_phi_deg,trk_the_deg);


  
  
  
  
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

  std::vector<int> mclist;
  std::vector<int> digilist;
  if (nclust > 0 ) { 
    const TVector3 clust_pos = main_clust->where();
    const double clust_energy = main_clust->energy();
    const double clust_phi = clust_pos.Phi();
    const double clust_the = clust_pos.Theta();
    mclist = main_clust->GetMcList();
    digilist = main_clust->DigiList();
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

  //int iii=0;
  //cout << "MCList.size()= " << mclist.size() << endl;
  //for (std::vector<int>::iterator it=mclist.begin(); it!=mclist.end(); ++it) {
  //  cout << "McList" << iii << " : " << *it << endl;
  //  ++iii;
  //}  


  
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

  double _reco_hit_energy = 0;
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

    int iclass = classify(reco_hit_energy/trk_energy);
    
    h_energy_b->Fill(reco_hit_energy/trk_energy);
    h_energy_corr_b->Fill(reco_hit_energy_corr/trk_energy);

    _reco_hit_energy = reco_hit_energy;
    
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

  for (int i=0; i<ncl; ++i) {
    h_energy_class[i]->Write();
    h_energy_class_nopi[i]->Write();    
    h_eratio_vs_r[i]->Write();
    h_eratio_vs_e_child_pip[i]->Write();
    h_eratio_vs_mom_child_pip[i]->Write();
  }
  
  for (int ing=0; ing<ngb; ++ing) {
    for (int inp=0; inp<npb; ++inp) {
      for (int inn=0; inn<nnb; ++inn) {
	h_energy_class_4d[ing][inp][inn]->Write();
      }
    }
  }

  h_r_ng->Write();
  h_r_nep->Write();
  h_r_nem->Write();
  h_r_np->Write();
  h_r_nn->Write();
  h_r_npip->Write();
  h_r_npim->Write();
  h_r_npi0->Write();
  h_r_nmum->Write();
  h_r_nmup->Write();
  h_r_nnumu->Write();
  h_r_nanumu->Write();
  h_c_ng->Write();
  h_c_nep->Write();
  h_c_nem->Write();
  h_c_np->Write();
  h_c_nn->Write();
  h_c_npip->Write();
  h_c_npim->Write();
  h_c_npi0->Write();
  h_c_nmum->Write();
  h_c_nmup->Write();
  h_c_nnumu->Write();
  h_c_nanumu->Write();

  

}


void HistFillerTask::fill_count_hists() {
  h_r_ng->Fill(r_ng);
  h_r_nep->Fill(r_nep);
  h_r_nem->Fill(r_nem);
  h_r_np->Fill(r_np);
  h_r_nn->Fill(r_nn);
  h_r_npip->Fill(r_npip);
  h_r_npim->Fill(r_npim);
  h_r_npi0->Fill(r_npi0);
  h_r_nmum->Fill(r_nmum);
  h_r_nmup->Fill(r_nmup);
  h_r_nnumu->Fill(r_nnumu);
  h_r_nanumu->Fill(r_nanumu);
  h_c_ng->Fill(c_ng);
  h_c_nep->Fill(c_nep);
  h_c_nem->Fill(c_nem);
  h_c_np->Fill(c_np);
  h_c_nn->Fill(c_nn);
  h_c_npip->Fill(c_npip);
  h_c_npim->Fill(c_npim);
  h_c_npi0->Fill(c_npi0);
  h_c_nmum->Fill(c_nmum);
  h_c_nmup->Fill(c_nmup);
  h_c_nnumu->Fill(c_nnumu);
  h_c_nanumu->Fill(c_nanumu);
}

const int HistFillerTask::classify(double ratio) {

  const int _ntrk = track_array->GetEntriesFast();
  if (debug>3)
    cout << "ntrk= " << _ntrk << endl;
  //cout << "========================================================================" << endl;
  initialize_remainder();
  initialize_child();

  int iclass = 0;
  int isubclass = 0;  

  double radius = 0.0;
  double radius_e[4] = {0.0};
  double e_ch_pi = 0.0;
  double mom_ch_pi = 0.0;
  bool fill4d = true;
  
  std::vector<int> child_ids;
  std::vector<int> remainder_ids;
  
  for (int itrk=0; itrk<_ntrk; ++itrk) {

    const PndMCTrack* _track = (PndMCTrack*) track_array->At(itrk);
    const TVector3 _trk_mom = _track->GetMomentum();
    const TLorentzVector _trk_mom4 = _track->Get4Momentum();
    const double _trk_mom_mag = _trk_mom.Mag();
    const TVector3 _trk_vtx = _track->GetStartVertex();
    const double _trk_vtxr = TMath::Hypot(_trk_vtx(0),_trk_vtx(1));
    const double _trk_vtxz = _trk_vtx(2);
    const double _trk_energy = _trk_mom4.E();
    
    const int _pdg = _track->GetPdgCode();
    const int _motherid = _track->GetMotherID();
    const int _grammaid = _track->GetSecondMotherID();

    if (_motherid==-1) {
      // skip
    } else if (_motherid==0) {
      child_ids.push_back( _pdg);
      if (abs(_pdg)!=11) {
	radius = _trk_vtxr;
      }
      if (_pdg==11) {
	radius_e[c_nem] = _trk_vtxr;
      }
      if (_pdg == 211) {
	mom_ch_pi = _trk_mom_mag;
	e_ch_pi = _trk_energy;
      }
      increment_child(_pdg);
    } else {
      remainder_ids.push_back( _pdg);      
      increment_remainder(_pdg);
    }
  }

  // pi0+EMShower
  if ( c_npi0 > 0 ) {
    fill4d = false;
    // case where there is atleast 1pi0 in children list 
    if ( (r_nep+r_nem)>=5 || r_ng >= 5)  {
      // associated with an EM shower (dominant in 0.9+)
      iclass = 0;
    } else {
      // not associated with an EM shower
      // this cas is very rare, so it will be treated same case as previous
      iclass = 0;
    }
  }
  // case of no child with no other remainder particles (dominant at mip)
  else if ( child_ids.size()==0 && remainder_ids.size()==0) {
    iclass = 1;
    fill4d = false;
  }
  else if ( child_ids.size()==1 && c_npip >= 1) {
    fill4d = false;
    // only one pi+ in the children list 
    if (remainder_ids.size()==0) {
      // and nothing in remainder list
      iclass = 2;
    } else {
      // something in the remainder list
      if ( r_npip > 0 ) {
	// atleast one pi+ in the remainder list
	iclass = 3;
      } else {
	// no pi+ in the remainder list
	iclass = 4;
      }
    }
  }

  else if ( c_nmup>0||c_nmum>0||c_nnumu>0||c_nanumu>0||r_nmup>0||r_nmum>0||r_nnumu>0||r_nanumu>0) {
    // a muon/neutrino muon has been created
    iclass = 5;
    fill4d = false;
  }

  // this is the case of pair absorption pi+nn->pn
  else {
    iclass = 6;
    fill4d = true;
    if ( c_nem==0 ) {

      //if ( child_ids.size()==2 && c_np == 1 && c_nn ==1 ) {
      //
      //	if (remainder_ids.size()==0) {
      //	  isubclass = 0;
      //	} else {
      //	  isubclass = 2;
      //	}
      //
      //} else {
      //}
      
      if ( c_np == 0 ) {
	isubclass = 0;	  
      } else {
	isubclass = 1;
      }

    } else {

      bool _elec_diff_rad = false;
      //for (int i=0; i<4; ++i) {
      //	if (radius_e[i]!=0 && radius != radius_e[i] ) {
      //	  _elec_diff_rad = true;
      //	  break;
      //	}
      //}
      if ( c_nem>0 ) {
	isubclass = 2;
      } else {
	isubclass = 3;
      }

    } 

  }

  if (ratio_min<ratio && ratio<ratio_max) {
    //if (iclass==6) {    
    //if (isubclass == 1 )
      cout << "===========================================" << endl;
      //else
      //cout << "!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=-" << endl;
    cout << "ratio = " << ratio << endl;
    shout_child();
    shout_remainder();
    //print_event();
  }
  
  h_energy_class[iclass]->Fill(ratio);
  h_eratio_vs_r[iclass]->Fill(ratio,radius);
  h_eratio_vs_e_child_pip[iclass]->Fill(ratio,e_ch_pi);
  h_eratio_vs_mom_child_pip[iclass]->Fill(ratio,mom_ch_pi);

  //int c_ne = c_nep+c_nem;
  int c_ne = c_ng;

  if (fill4d) {
    h_energy_class_4d[c_ng>=ngb?ngb-1:c_ng][c_np>=npb?npb-1:c_np][c_nn>=nnb?nnb-1:c_nn]->Fill(ratio);
    fill_count_hists();
    h_energy_class_nopi[isubclass]->Fill(ratio);
  }
    
  return iclass;
  
}


const TString HistFillerTask::tpdg(int ipdg) {

  if (ipdg == 22) {
    return "g";
  } else if (ipdg == 11) {
    return "e-";
  } else if (ipdg == -11) {
    return "e+";
  } else if (ipdg == 2212) {
    return "p";
  } else if (ipdg == 2112) {
    return "n";
  } else if (ipdg == 211) {
    return "pi+";
  } else if (ipdg == -211) {
    return "pi-";
  } else if (ipdg == 13) {
    return "mu-";
  } else if (ipdg == -13) {
    return "mu+";
  } else if (ipdg == 14) {
    return "numu";
  } else if (ipdg == -14) {
    return "anumu";
  } else if (ipdg == 111) {
    return "pi0";
  } else {
    return Form("?%d",ipdg);
  }

}

void HistFillerTask::shout_child() {
  cout << "Children: ";
  if (c_ng>0) cout << c_ng << "g ";
  if (c_nep>0) cout << c_nep << "e+ ";
  if (c_nem>0) cout << c_nem << "e- ";  
  if (c_np>0) cout << c_np << "p ";
  if (c_nn>0) cout << c_nn << "n ";
  if (c_npi0>0) cout << c_npi0 << "pi0 ";
  if (c_npip>0) cout << c_npip << "pip ";
  if (c_npim>0) cout << c_npim << "pim ";
  if (c_nmup>0) cout << c_nmup << "mup ";
  if (c_nmum>0) cout << c_nmum << "mum ";
  if (c_nnumu>0) cout << c_nnumu << "numu ";
  if (c_nanumu>0) cout << c_nanumu << "anumu ";        
  cout << endl; 
}

void HistFillerTask::shout_remainder() {
  cout << "Remainder: ";
  if (r_ng>0) cout << r_ng << "g ";
  if (r_nep>0) cout << r_nep << "e+ ";
  if (r_nem>0) cout << r_nem << "e- ";  
  if (r_np>0) cout << r_np << "p ";
  if (r_nn>0) cout << r_nn << "n ";
  if (r_npi0>0) cout << r_npi0 << "pi0 ";
  if (r_npip>0) cout << r_npip << "pip ";
  if (r_npim>0) cout << r_npim << "pim ";
  if (r_nmup>0) cout << r_nmup << "mup ";
  if (r_nmum>0) cout << r_nmum << "mum ";
  if (r_nnumu>0) cout << r_nnumu << "numu ";
  if (r_nanumu>0) cout << r_nanumu << "anumu ";        
  cout << endl; 
}

void HistFillerTask::initialize_remainder() {
  r_ng = 0;
  r_nep = 0;
  r_nem = 0;  
  r_np = 0;
  r_nn = 0;
  r_npip = 0;
  r_npim = 0;
  r_npi0 = 0;
  r_nmum = 0;
  r_nmup = 0;
  r_nnumu = 0;
  r_nanumu = 0;    
}

void HistFillerTask::initialize_child() {
  c_ng = 0;
  c_nep = 0;
  c_nem = 0;  
  c_np = 0;
  c_nn = 0;
  c_npip = 0;
  c_npim = 0;
  c_npi0 = 0;
  c_nmum = 0;
  c_nmup = 0;
  c_nnumu = 0;
  c_nanumu = 0;    
}


void HistFillerTask::increment_remainder(int _pdg) {
  if (_pdg == 22) {
    ++r_ng;
  } else if (_pdg == 11) {
    ++r_nem;
  } else if (_pdg == -11) {
    ++r_nep;
  } else if (_pdg == 2212) {
    ++r_np;
  } else if (_pdg == 2112) {
    ++r_nn;
  } else if (_pdg == 111) {
    ++r_npi0;
  } else if (_pdg == 211) {
    ++r_npip;
  } else if (_pdg == -211) {
    ++r_npim;
  } else if (_pdg == 13) {
    ++r_nmum;
  } else if (_pdg == -13) {
    ++r_nmup;
  } else if (_pdg == 14) {
    ++r_nnumu;
  } else if (_pdg == -14) {
    ++r_nanumu;
  } else {
    //cout << "HistFillerTask::increment_child -- Unkown pdgid " << _pdg << endl;
  }
}

void HistFillerTask::increment_child(int _pdg) {
  if (_pdg == 22) {
    ++c_ng;
  } else if (_pdg == 11) {
    ++c_nem;
  } else if (_pdg == -11) {
    ++c_nep;
  } else if (_pdg == 2212) {
    ++c_np;
  } else if (_pdg == 2112) {
    ++c_nn;
  } else if (_pdg == 111) {
    ++c_npi0;
  } else if (_pdg == 211) {
    ++c_npip;
  } else if (_pdg == -211) {
    ++c_npim;
  } else if (_pdg == 13) {
    ++c_nmum;
  } else if (_pdg == -13) {
    ++c_nmup;
  } else if (_pdg == 14) {
    ++c_nnumu;
  } else if (_pdg == -14) {
    ++c_nanumu;
  } else {
    //cout << "HistFillerTask::increment_child -- Unkown pdgid " << _pdg << endl;
  }
}

void HistFillerTask::print_event() {

  int ntrk = track_array->GetEntriesFast();
  for (int itrk=0; itrk<ntrk; ++itrk) {
    const PndMCTrack* _track = (PndMCTrack*) track_array->At(itrk);
    const TVector3 _trk_mom = _track->GetMomentum();
    const TLorentzVector _trk_mom4 = _track->Get4Momentum();
    const double _trk_mom_mag = _trk_mom.Mag();
    const TVector3 _trk_vtx = _track->GetStartVertex();
    const double _trk_vtxr = TMath::Hypot(_trk_vtx(0),_trk_vtx(1));
    const double _trk_vtxz = _trk_vtx(2);
    const double _trk_energy = _trk_mom4.E();
    const int _pdg = _track->GetPdgCode();
    const int _motherid = _track->GetMotherID();

    //if (std::find(mclist.begin(),mclist.end(),itrk) != mclist.end() ) cout << " -> ";
    cout << "Track " << itrk << ": pdg= " << tpdg(_pdg) << " pID= " << _motherid
	 << " vR= " << _trk_vtxr
      //<< " vZ= " << _trk_vtxz
	 << " |mom|= " << _trk_mom_mag
	 << " E= " << _trk_energy
      //<< " M= " << sqrt(_trk_energy*_trk_energy-_trk_mom_mag*_trk_mom_mag)
	 << endl;
  }
}


//------ run over the digis that contributed to this cluster ---- //
//std::vector<int> mc_from_digi;
//std::vector<int> mc_from_hit;  
//std::vector<Int_t>::iterator digipos;
//double etot_digi = 0.0;
//for (digipos=digilist.begin();digipos!=digilist.end();++digipos){
//  PndEmcDigi *digi = (PndEmcDigi *) digi_array->At(*digipos);
//  etot_digi+= digi->GetEnergy();
//  TVector3 pos = digi->where();
//  double _rad = TMath::Hypot(pos(0),pos(1));
//  //cout << " digid " << *digipos << " track= " << digi->GetTrackId() << " r= " << _rad <<  " E= " << digi->GetEnergy() << " HitId= " << digi->GetHitIndex() <<endl;
//  mc_from_digi.push_back(digi->GetTrackId());
//
//  /*
//	PndEmcHit *hit = (PndEmcHit *) hit_array->At(digi->GetHitIndex());
//	mc_from_hit = hit->GetMcList();
//	for (int ii=0; ii<mc_from_hit.size(); ++ii) {
//	cout << "     * hit " << ii << " : track= " <<  hit[ii] << endl;
//	}
//  */
//}
//cout << "Ratio = " << ratio << endl;
//cout << "E_tot(digi) = " << etot_digi << " E_reco= " << _reco_hit_energy << endl;


  //else if ( c_npip==1 ) {
  //
  //  if (child_ids.size() == 1) {
  //
  //    if (remainder_ids.size()==0) {
  //	iclass = 3;
  //    }
  //
  //    else if (remainder_ids.size()<=4) {
  //	iclass = 4;
  //    }
  //
  //    else {
  //	iclass = 5;
  //    }
  //  }
  //
  //  else if (child_ids.size() == 2) {
  //
  //    if ( c_np == 1 ) {
  //	iclass = 6;
  //    }
  //
  //    else if ( c_nn == 1) {
  //	iclass = 7;
  //    }
  //
  //    else if ( c_ng == 1) {
  //	iclass = 8;
  //    }
  //
  //    else {
  //	iclass = 9;
  //    }
  //
  //  }
  //
  //  else {
  //    if ( c_nem == 0 && c_nep == 0 ) {
  //	iclass = 10;
  //    } else {
  //	iclass = 11;
  //    }
  //  }
  //}
  //
  //else if ( c_npip>=2) {
  //
  //  if (child_ids.size() == 2) {
  //
  //    if (remainder_ids.size()==0) {
  //	iclass = 12;
  //    }
  //
  //    if (remainder_ids.size()<=5) {
  //	iclass = 12;
  //    }      
  //
  //    else {
  //	iclass = 12;
  //    }
  //  }
  //
  //  else {
  //    iclass = 12;
  //  }
  //
  //}
  
  //else if ( c_npip==1 && child_ids.size()==1 && remainder_ids.size()==0 ) {
  //  iclass = 3;
  //}
  //
  //else if ( c_npip==1 && child_ids.size()==1 ) {
  //  iclass = 4;
  //}  
  //
  //else if ( c_npip==1 && (r_ng>2||r_nn>2||r_np>2) ) {
  //  iclass = 5;
  //}

  //else if ( c_npip==1 && r_np>2 ) {
  //  iclass = 6;
  //}

  //else if ( child_ids.size()>3) {
  //  iclass = 7;
  //}
  //
  //else if ( remainder_ids.size()>3) {
  //  iclass = 8;
  //}
  
  //else if ( remainder_ids.size()==0 && child_ids.size()>1 && c_npip==1 ) {
  //  // there is a pi+ in the child list associated with 1p and/or 1n and/or 1gamma and no remainder - dominant in 0.7 - 0.75)
  //  if ((c_np == 0 && c_nn <= 1)||(c_nn == 0 && c_np <= 1)) {
  //    iclass = 3;
  //  }
  //
  //  else if ((c_np == 1 && c_nn <= 1)||(c_nn == 1 && c_np <= 1)) {
  //    iclass = 4;
  //  }
  //
  //  // there is a pi+ in the child list associated with 1p and/or 1n and/or 1gamma and no remainder - dominant in 0.7 - 0.75)
  //  else if ((c_np == 2 && c_nn <= 1)||(c_nn == 2 && c_np <= 1)) {
  //    iclass = 5;
  //  }
  //
  //
  //  // there is a pi+ in the child list associated with 1p and/or 1n and/or 1gamma and no remainder - dominant in 0.7 - 0.75)
  //  else if ((c_np == 3 && c_nn <= 1)||(c_nn == 3 && c_np <= 1)) {
  //    iclass = 6;
  //  }
  //  
  //  else if (c_np>3 || c_nn>3) {
  //    iclass = 7;
  //  }
  //  
  //}
  
  //else {
  //  //if (ratio>0.55) cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  //  //else cout << "----------------------------------------------------------------------------------------------------" << endl;
  //  //cout << "Ratio= " << ratio << endl;
  //  //shout_child();
  //  //shout_remainder();
  //  //print_event();
  //  //iclass = 8;
  //
  //  if (c_npip==0) {
  //    iclass = 9;
  //  } else {
  //    iclass = 10;
  //  }
  //}

  //cout << "r= " << ratio << " rad= " << radius << " e= " << e_ch_pi << " p= " << mom_ch_pi << endl;



ClassImp(HistFillerTask)

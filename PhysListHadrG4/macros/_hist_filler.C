//#include "TFile.h"
//#include "TTree.h"
//#include "TClonesArray.h"
//#include "TVector3.h"
//
////#include "../../emc/EmcReco/PndEmcCluster.h"
////#include "../../emc/EmcReco/PndMCTrack.h"
////#include "../../emc/EmcTools/PndEmcMapper.h"
////
//#include <iostream>

{

  string data_file_name = "output/data_QGSP_BERT_EMV_MOM_1.5_pip.root";
  
  int debug = 5;
  
  TFile *data_file = TFile::Open(data_file_name.c_str());
  TTree *t = (TTree*) data_file->Get("cbmsim");

  TClonesArray *clust_array = new TClonesArray("PndEmcCluster");
  t->SetBranchAddress("EmcCluster",&clust_array);

  TClonesArray *track_array = new TClonesArray("PndMCTrack");
  t->SetBranchAddress("MCTrack",&track_array);

  TClonesArray *bump_array = new TClonesArray("PndEmcBump");
  t->SetBranchAddress("EmcBump",&bump_array);

  //TClonesArray *hit_array = new TClonesArray("PndEmcHit");
  //t->SetBranchAddress("EmcHit",&hit_array);

  TClonesArray *reco_hit_array = new TClonesArray("PndEmcRecoHit");
  t->SetBranchAddress("EmcRecoHit",&reco_hit_array);
  
  PndEmcMapper::Init(1);

  int nevt = t->GetEntriesFast();
  for (int ievt=0; ievt<nevt; ++ievt) {
    t->GetEntry(ievt);

    if (ievt>100) break;
    if (debug>1) cout << "=============== Event " << ievt << "  ====================" << endl;

    const int ntrk = track_array->GetEntriesFast();
    //if (debug>1) cout << "ntrk= " << ntrk << endl;
    const PndMCTrack* track = (PndMCTrack*) track_array->At(0);
    const TVector3 trk_mom = track->GetMomentum();
    const double trk_phi = trk_mom.Phi();
    const double trk_the = trk_mom.Theta();
    const double trk_mom_mag = trk_mom.Mag();
    const TVector3 trk_vtx = track->GetStartVertex();
    const double trk_vtxr = TMath::Hypot(trk_vtx(0),trk_vtx(1));
    const double trk_vtxz = trk_vtx(2);

    const int pdg = track->GetPdgCode();
    const int motherid = track->GetMotherID();
    const int grammaid = track->GetSecondMotherID();
      
    if (debug>1) cout << "Track: pdg= " << pdg << " pID= " << motherid
		      << " gpID= " << grammaid << " vR= " << trk_vtxr << " vZ= " << trk_vtxz
		      << " phi= " << trk_phi << " the= " << trk_the << " |mom|= " << trk_mom_mag << endl;



    //-------- clusters -----------//
    const int nclust = clust_array->GetEntriesFast();
    //if (debug>1) cout << "nclust= " << nclust << endl;
    double del_the_min = 1e10;
    PndEmcCluster* clust;
    double energy_tot_clust = 0.0;
    cout << "Nclust= " << nclust;
    for (int iclust=0; iclust<nclust; ++iclust) {
      PndEmcCluster *tmp_clust = (PndEmcCluster*)clust_array->At(iclust);
      energy_tot_clust += tmp_clust->energy();
      cout << " [E_" << iclust << "= " << tmp_clust->energy() << ", Th_" << iclust << "= " << tmp_clust->where().Theta() << "]   ";
      if ( del_the_min > fabs(tmp_clust->where().Theta()-trk_the) ) {
	clust = tmp_clust;
	del_the_min = fabs(tmp_clust->where().Theta()-trk_the);
      }
    }
    cout << endl;
    const TVector3 clust_pos = clust->where();
    const double clust_energy = clust->energy();
    const double clust_phi = clust_pos.Phi();
    const double clust_the = clust_pos.Theta();
    cout << "Clust: E_clust = " << clust_energy << " phi= " << clust_phi << " the= " << clust_the << " E_tot= " << energy_tot_clust <<endl;



    //---------- bump ---------//
    const int nbump = bump_array->GetEntriesFast();
    //if (debug>1) cout << "nbump= " << nbump << endl;
    double del_the_min = 1e10;
    PndEmcBump* bump;
    double energy_tot_bump = 0.0;
    cout << "Nbump= " << nbump;
    for (int ibump=0; ibump<nbump; ++ibump) {
      PndEmcBump *tmp_bump = (PndEmcBump*)bump_array->At(ibump);
      energy_tot_bump += tmp_bump->energy();
      cout << " [E_" << ibump << "= " << tmp_bump->energy() << ", Th_" << ibump << "= " << tmp_bump->where().Theta() << "]   ";
      if ( del_the_min > fabs(tmp_bump->where().Theta()-trk_the) ) {
	bump = tmp_bump;
	del_the_min = fabs(tmp_bump->where().Theta()-trk_the);
      }
    }
    cout << endl;
    const TVector3 bump_pos = bump->where();
    const double bump_energy = bump->energy();
    const double bump_phi = bump_pos.Phi();
    const double bump_the = bump_pos.Theta();
    cout << "Bump: E_bump = " << bump_energy << " phi= " << bump_phi << " the= " << bump_the << " E_tot= " << energy_tot_bump <<endl;


    ////---------- hit ---------//
    //const int nhit = hit_array->GetEntriesFast();
    ////if (debug>1) cout << "nhit= " << nhit << endl;
    //double del_the_min = 1e10;
    //PndEmcHit* hit;
    //double energy_tot_hit = 0.0;
    //cout << "Nhit= " << nhit;
    //for (int ihit=0; ihit<nhit; ++ihit) {
    //  PndEmcHit *tmp_hit = (PndEmcHit*)hit_array->At(ihit);
    //  energy_tot_hit += tmp_hit->GetEnergy();
    //  cout << " [E_" << ihit << "= " << tmp_hit->GetEnergy() << ", Th_" << ihit << "= " << tmp_hit->GetTheta() << "]   ";
    //  if ( del_the_min > fabs(tmp_hit->GetTheta()-trk_the) ) {
    //	hit = tmp_hit;
    //	del_the_min = fabs(tmp_hit->GetTheta()-trk_the);
    //  }
    //}
    //cout << endl;
    //const double hit_energy = hit->GetEnergy();
    //const double hit_phi = hit.GetPhi();
    //const double hit_the = hit.GetTheta();
    //cout << "Hit: E_hit = " << hit_energy << " phi= " << hit_phi << " the= " << hit_the << " E_tot= " << energy_tot_hit <<endl;


    //---------- reco_hit ---------//
    const int nreco_hit = reco_hit_array->GetEntriesFast();
    //if (debug>1) cout << "nreco_hit= " << nreco_hit << endl;
    double del_the_min = 1e10;
    PndEmcRecoHit* reco_hit;
    double energy_tot_reco_hit = 0.0;
    double energy_tot_reco_hit_corr = 0.0;
    cout << "Nreco_hit= " << nreco_hit;
    for (int ireco_hit=0; ireco_hit<nreco_hit; ++ireco_hit) {
      PndEmcRecoHit *tmp_reco_hit = (PndEmcRecoHit*)reco_hit_array->At(ireco_hit);
      energy_tot_reco_hit += tmp_reco_hit->GetEnergy();
      energy_tot_reco_hit_corr += tmp_reco_hit->GetEnergyCorrected();
      cout << " [E_" << ireco_hit << "= " << tmp_reco_hit->GetEnergy() << ", Th_" << ireco_hit << "= " << tmp_reco_hit->GetPosition().Theta() << "]   ";
      if ( del_the_min > fabs(tmp_reco_hit->GetPosition().Theta()-trk_the) ) {
	reco_hit = tmp_reco_hit;
	del_the_min = fabs(tmp_reco_hit->GetPosition().Theta()-trk_the);
      }
    }
    cout << endl;
    const TVector3 reco_hit_pos = reco_hit->GetPosition();
    const double reco_hit_energy = reco_hit->GetEnergy();
    const double reco_hit_energy_corr = reco_hit->GetEnergyCorrected();
    const double reco_hit_phi = reco_hit_pos.Phi();
    const double reco_hit_the = reco_hit_pos.Theta();
    cout << "Reco_Hit: E_reco_hit = " << reco_hit_energy << " E_reco_hit_corr = " << reco_hit_energy_corr  << " phi= " << reco_hit_phi << " the= " << reco_hit_the << " E_tot= " << energy_tot_reco_hit << " E_tot_corr= " << energy_tot_reco_hit_corr <<endl;

    
    

    




    
  }

}

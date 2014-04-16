#ifndef HistFillerTask_H
#define HistFillerTask_H 1

#include "FairTask.h"
#include "TLorentzVector.h"

class TTree;
class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class HistFillerTask : public FairTask
{

public:
	
  // ** Default constructor   
  HistFillerTask();
	
  // ** Destructor 
  ~HistFillerTask();	
	
  // ** Virtual method Init 
  virtual InitStatus Init();
	
  // ** Virtual method Exec 
  virtual void Exec(Option_t* opt);
	
  virtual void Finish();

protected:

  void initialize_hists();
  
  int nevt;
  
  TTree *cbmsim;
  TClonesArray *track_array;
  TClonesArray *clust_array;
  TClonesArray *bump_array;
  TClonesArray *reco_hit_array;

  TH1F *h_energy_b;
  TH1F *h_energy_corr_b;

  TH1F *h_energy;
  TH1F *h_energy_corr;

  TH2F *h_phi_the_b;
  TH2F *h_phi_the;

  TH1F *h_dphi;
  TH1F *h_dthe;

  TH1F *h_nclust;
  TH1F *h_nbump;
  TH1F *h_nreco_hit;

  TH2F* h_extra_clust_dphi_dthe;
  TH2F* h_extra_bump_dphi_dthe;
  TH2F* h_extra_reco_hit_dphi_dthe;

  TH1F* h_extra_clust_energy;
  TH1F* h_extra_bump_energy;
  TH1F* h_extra_reco_hit_energy;

  TH1F* h_all_clust_energy;
  TH1F* h_all_bump_energy;
  TH1F* h_all_reco_hit_energy;

  TH1F* h_clust_tot_energy;
  TH1F* h_bump_tot_energy;
  TH1F* h_reco_hit_tot_energy;
  TH1F* h_reco_hit_tot_energy_corr;
  
  ClassDef(HistFillerTask,1);
  
};

#endif /* HistFillerTask_H */

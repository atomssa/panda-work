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
  HistFillerTask(double,double);
	
  // ** Destructor 
  ~HistFillerTask();	
	
  // ** Virtual method Init 
  virtual InitStatus Init();
	
  // ** Virtual method Exec 
  virtual void Exec(Option_t* opt);
	
  virtual void Finish();

protected:

  void initialize_hists();
  const TString tpdg(int ipdg);
  const int  classify(double);
  void initialize_remainder();
  void initialize_child();
  void increment_remainder(int);
  void increment_child(int);  
  void shout_remainder();
  void shout_child();
  void print_event();
  void fill_count_hists();

  int nevt;

  double ratio_min;
  double ratio_max;
  
  TTree *cbmsim;
  TClonesArray *track_array;
  TClonesArray *clust_array;
  TClonesArray *bump_array;
  TClonesArray *reco_hit_array;
  TClonesArray *digi_array;
  TClonesArray *hit_array;    

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

  static const int ncl = 13;
  TH1F* h_energy_class[ncl];
  
  TH1F* h_energy_class_nopi[ncl];

  static const int ngb = 2;
  static const int npb = 2;
  static const int nnb = 12;
  
  TH1F* h_energy_class_4d[ngb][npb][nnb];
  
  TH2F* h_eratio_vs_r[ncl];
  TH2F* h_eratio_vs_e_child_pip[ncl];
  TH2F* h_eratio_vs_mom_child_pip[ncl];  
  
  int r_ng;
  int r_nep;
  int r_nem;
  int r_np;
  int r_nn;
  int r_npip;
  int r_npim;
  int r_npi0;
  int r_nmum;
  int r_nmup;
  int r_nnumu;
  int r_nanumu;

  int c_ng;
  int c_nep;
  int c_nem;
  int c_np;
  int c_nn;
  int c_npip;
  int c_npim;
  int c_npi0;
  int c_nmum;
  int c_nmup;
  int c_nnumu;
  int c_nanumu;
  
  TH1F* h_r_ng;
  TH1F* h_r_nep;
  TH1F* h_r_nem;
  TH1F* h_r_np;
  TH1F* h_r_nn;
  TH1F* h_r_npip;
  TH1F* h_r_npim;
  TH1F* h_r_npi0;
  TH1F* h_r_nmum;
  TH1F* h_r_nmup;
  TH1F* h_r_nnumu;
  TH1F* h_r_nanumu;

  TH1F* h_c_ng;
  TH1F* h_c_nep;
  TH1F* h_c_nem;
  TH1F* h_c_np;
  TH1F* h_c_nn;
  TH1F* h_c_npip;
  TH1F* h_c_npim;
  TH1F* h_c_npi0;
  TH1F* h_c_nmum;
  TH1F* h_c_nmup;
  TH1F* h_c_nnumu;
  TH1F* h_c_nanumu;
  
  int debug;

  ClassDef(HistFillerTask,1);
  
};

#endif /* HistFillerTask_H */

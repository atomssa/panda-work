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
  
  ClassDef(HistFillerTask,1);
  
};

#endif /* HistFillerTask_H */

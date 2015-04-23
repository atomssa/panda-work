#ifndef EffHists_H
#define EffHists_H

#include "FairTask.h"
#include "PndAnalysis.h"
#include "PndPidProbability.h"

#include <vector>
#include <string>

class TH1F;
class TH2F;
class TEfficiency;
class TVector3;
class TFile;
class TString;
class TLorentzVector;
class PndPidProbability;
class FairRootManager;

class EffHists : public FairTask{

  typedef Double_t (PndPidProbability::*prob_func)(PndPidProbability*) const;

 public:
  EffHists(int);
  ~EffHists();

  virtual InitStatus Init();

  virtual void Exec(Option_t* opt);

  virtual void FinishTask();
  virtual void FinishEvent();

  void set_verbosity(int v) {verb = v;}
  void set_prob_cut(int, double);

 private:
  void init_hists();
  void write_hists();
  TClonesArray* init_tca(TString);
  InitStatus init_tcas();

  double get_comb_prob(prob_func);

 private:
  int verb;
  int nevt;
  int m_sp;

  PndAnalysis *fAna;
  FairRootManager* m_ioman;

  TClonesArray* m_cand_array;
  TClonesArray* m_drc_array;
  TClonesArray* m_disc_array;
  TClonesArray* m_mvd_array;
  TClonesArray* m_stt_array;
  TClonesArray* m_emcb_array;
  PndPidProbability *m_prob_emcb;
  PndPidProbability *m_prob_stt;
  PndPidProbability *m_prob_mvd;
  PndPidProbability *m_prob_drc;
  PndPidProbability *m_prob_disc;

 public:
  enum{iposit=0, imuonp, ipionp, ikaonp, iproton, ielec, imuonm, ipionm, ikaonm, iantiproton, nsp_max};
  static const TString s_spc[nsp_max];
  static const TString s_spc_tex[nsp_max];

  enum{iel=0,imu,ipi,ik,iprot,npid_max};
  static const TString s_pid[npid_max]; // = {"e_id", "mu_id", "pi_id", "k_id", "prot_id"};

  enum{iemc = 0, istt, imvd, idirc, idisc, ndet};
  static const TString s_det[ndet];

 private:
  double prob_cut[npid_max]; // = {0.5, 0.5, 0.5, 0.5, 0.5}
  double det_var_max[ndet]; // = {"emc", "stt", "mvd", "dirc", "disc"};
  TH2F* eff_den[npid_max], *eff_num[npid_max];
  TEfficiency *eff2d[npid_max], *eff1d_the[npid_max], *eff1d_mom[npid_max];
  double mom_max;
  double the_max;

  ClassDef(EffHists,1);

};

#endif /*EffHists_H*/

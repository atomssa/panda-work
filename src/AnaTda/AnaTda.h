#ifndef AnaTda_H
#define AnaTda_H 1

#include <vector>

#include "FairTask.h"
#include "TLorentzVector.h"
#include "RhoCandList.h"

class TTree;
class TClonesArray;
class TObjectArray;
class TProfile;
class TProfile2D;
class TH2F;
class TH1F;
class TFile;
class TLorentzVector;
class TVector3;

class PndAnalysis;
class RhoCandList;
class RhoCandidate;
class RhoMassParticleSelector;

class AnaTda : public FairTask {
  public:
  // ** Default constructor
  AnaTda(const int&);

  // ** Destructor
  ~AnaTda();

  // ** Virtual method Init
  virtual InitStatus Init();

  // ** Virtual method Exec
  virtual void Exec(Option_t* opt);

  virtual void Finish();

  protected:
  PndAnalysis *fAnalysis;             // *** the PndAnalysis object

  int nevt;

  bool fBremCorr;

  int SelectTruePid(PndAnalysis *, RhoCandList &);
  void calc_kin(RhoCandidate*, RhoCandidate *, double &, double &, double &);
  void def_hists();
  void initial_state();
  void set_selectors();

  void def_tutorial_hists();
  void def_manual_kin_fit_hists();
  void def_pair_hists();
  void def_single_hists();
  void def_kin_fit_hists(const int& type, const int&);
  void def_kin_fit_hists();

  void print_rho_cand_list(RhoCandList&, const char*);
  void cleanup_rho_cand_lists();
  void get_singles_lists();
  void fill_single_dists();
  void truth_match(RhoCandList& org, RhoCandList& dest, const int &pdg);
  void truth_match_singles();
  void fill_single_dists_tr();
  void make_pair_lists();
  void fill_pair_mass(RhoCandList& org, TH1F* dest);
  void fill_pair_dists();
  void fill_pi0s();
  void pdgm_nearest_pi0s();
  void jpsi_mass_selection();
  void pi0_kinematic_selection();
  void jpsi_truth_match();
  void kin_fit_full_sys(RhoCandList& org, const int&, const int&);
  void kin_fit_epem_pi0_btb();
  void kin_fit_epem_pi0_cts();
  void kin_fit_pippim_pi0_btb();
  void kin_fit_pippim_pi0_cts();
  void kin_fit_epem_pi0_nearest();
  void kin_fit_pippim_pi0_nearest();

  RhoCandList pip, pim, ep, em, g1, g2;
  RhoCandList pip_tr, pim_tr, ep_tr, em_tr, g1_tr, g2_tr;
  RhoCandList epem, pippim, gg;
  RhoCandList epem_tr, pippim_tr, gg_tr;
  RhoCandList pi0, pi0_true, pi0nearest, pi0_btb, pi0_cts;

  RhoCandList epem_mcut, pippim_mcut;
  RhoCandList jpsi, jpsi_true;

  RhoCandList epem_mcut_pi0_btb, epem_mcut_pi0_cts;
  RhoCandList pippim_mcut_pi0_btb, pippim_mcut_pi0_cts;
  RhoCandList epem_pi0nearest;
  RhoCandList pippim_pi0nearest;

  double m0_jpsi;
  double m0_pi0;
  RhoMassParticleSelector *jpsiMassSel;
  RhoMassParticleSelector *pi0MassSel;
  TLorentzVector ini;
  TVector3 boost_to_cm;
  double sqrt_s;

  // *** create some histograms
  TH1F *hjpsim_all;
  TH1F *hpi0m_all;
  TH1F *hjpsipi0m_all;

  TH1F *hjpsim_ftm;
  TH1F *hpi0m_ftm;
  TH1F *hjpsipi0m_ftm;

  TH1F *hjpsim_true;
  TH1F *hpi0m_true;
  TH1F *hjpsipi0m_true;

  TH1F *hjpsim_nearest;
  TH1F *hpi0m_nearest;
  TH1F *hjpsipi0m_nearest;

  TH1F *hjpsim_sel;
  TH1F *hpi0m_sel;
  TH1F *hjpsipi0m_sel;

  TH1F *hjpsim_nm;
  TH1F *hpi0m_nm;
  TH1F *hjpsipi0m_nm;

  TH1F *hjpsim_diff;
  TH1F *hpi0m_diff;
  TH1F *hjpsipi0m_diff;

  TH1F *h_m_epem;
  TH1F *h_m_epem_tr;
  TH1F *h_m_epem_mcut;
  TH1F *h_m_pippim;
  TH1F *h_m_pippim_tr;
  TH1F *h_m_gg;
  TH1F *h_m_gg_tr;

  TH1F *h_m_pi0n;
  TH1F *h_m_epem_pi0n;
  TH1F *h_m_pippim_pi0n;

  TH1F *h_4c_chi2[3][2];
  TH1F *h_4c_prob[3][2];
  TH1F *h_4c_m[3][2];
  TH2F *h_4c_prob_vs_m[3][2];

  TH1F* h_num_g;
  TH1F* h_num_epm;
  TH1F* h_num_pipm;
  TH1F* h_e_g;
  TH2F* h_mom_the_epm;
  TH2F* h_mom_the_pipm;
  TH1F* h_num_g_tr;
  TH1F* h_num_epm_tr;
  TH1F* h_num_pipm_tr;
  TH1F* h_e_g_tr;
  TH2F* h_mom_the_epm_tr;
  TH2F* h_mom_the_pipm_tr;

  TH1F *h_oa_gg_epair;
  TH1F *h_mass_gg_epair;
  TH2F *h_oa_gg_epair_vs_mass_gg;
  TH2F *h_mass_gg_epair_vs_mass_gg;
  TH2F *h_oa_vs_mass_gg_epair;
  TH2F *h_oa_vs_mass_gg_epair_btb;  // most back-to-back
  TH2F *h_oa_vs_mass_gg_epair_cts;  // closest-to-s
  TH1F *h_m_gg_btb;  // most back-to-back
  TH1F *h_m_gg_cts;  // closest-to-s

  ClassDef(AnaTda, 1);
};

#endif /* AnaTda_H */

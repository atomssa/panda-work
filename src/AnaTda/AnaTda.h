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

  std::vector<double> pi0oacut;
  std::vector<double> pi0ecut_min;
  std::vector<double> pi0ecut_max;
  double pi0mcut_min;
  double pi0mcut_max;
  double jpsi_mcut_min;
  double jpsi_mcut_max;

  void calc_kin(RhoCandidate*, RhoCandidate *, double &, double &, double &);
  void def_hists();
  void initial_state();
  void set_selectors();

  void def_tutorial_hists();
  void def_manual_kin_fit_hists(const int&);
  void def_pair_hists();
  void def_single_hists();
  void def_kin_fit_hists(const int& type, const int&);
  void def_kin_fit_hists();
  void def_gamma_from_pi0_hists();
  void def_elecs_from_jpsi_hists();
  void write_tut_hists();
  void write_manual_kin_fit_hists(const int&);
  void write_kin_fit_hists();

  double oa(RhoCandidate*, RhoCandidate*);
  double mass(RhoCandidate*, RhoCandidate*);

  void fill_single_e(RhoCandList&, TH1*);
  void fill_single_e(RhoCandList&, RhoCandList&, TH1*);
  void fill_single_mom(RhoCandList&, TH1*);
  void fill_single_mom(RhoCandList&, RhoCandList&, TH1*);
  void fill_single_mom_the(RhoCandList&, TH2F*);
  void fill_single_mom_the(RhoCandList&, RhoCandList&, TH2F*);

  void print_mc_list();
  void print_rho_cand_list(RhoCandList&, const char*);
  void cleanup_rho_cand_lists();

  void get_singles_lists();
  void fill_single_dists();

  void truth_match(RhoCandList& org, RhoCandList& dest, const int &pdg);
  void truth_match_singles();
  void fill_single_dists_tr();

  // Pull out candidates
  void make_pair_lists();
  void fill_pair_mass(RhoCandList& org, TH1F* dest);
  void fill_pair_oa(RhoCandList& org, TH1F* dest);
  void fill_pair_dists();

  void primary_match(const int&, RhoCandidate*, RhoCandidate*, RhoCandidate*, RhoCandidate*);
  void primary_match(const int&, RhoCandidate*, RhoCandidate*);
  bool primary_match_pair(const int&, RhoCandidate*, double&);
  bool primary_match_single(const int&, RhoCandidate*, double&);

  bool primary_match_pi0(RhoCandidate* c, double&);
  void pi0_truth_match();
  void fill_gamma_from_pi0s();
  void pi0_analysis_cut(RhoCandList&, const int&);
  void pi0_analysis_cut();
  void pi0_analysis_cut(RhoCandList&, RhoCandList&, const int&, bool);
  void fill_pi0_analysis_hists(RhoCandList&, const int&);
  void fill_pi0_analysis_hists(RhoCandList&, const int&, const int&);

  void jpsi_truth_match();
  bool primary_match_jpsi(RhoCandidate*, double&);
  void fill_elecs_from_jpsi();
  void jpsi_analysis_cut();
  void jpsi_analysis_cut(RhoCandList&, RhoCandList&);
  void fill_jpsi_analysis_hists(RhoCandList&, const int&);

  // has to be reworked
  void kin_fit_full_sys(RhoCandList& org, const int&, const int&);
  void kin_fit_all();
  void kin_fit_epem_pi0_btb();
  void kin_fit_epem_pi0_cts();
  void kin_fit_pippim_pi0_btb();
  void kin_fit_pippim_pi0_cts();
  void kin_fit_pi0_nearest(RhoCandList& org, const int&);
  void kin_fit_pi0_nearest_all();
  void kin_fit_epem_pi0_nearest();
  void kin_fit_pippim_pi0_nearest();

  void pi0_kinematic_selection(RhoCandList&, RhoCandList&, const int&);
  void pdgm_nearest_pi0s();


  RhoCandList mcList;
  RhoCandList pip, pim, ep, em, g1, g2;
  RhoCandList pip_tr, pim_tr, ep_tr, em_tr, g1_tr, g2_tr;
  RhoCandList epem, pippim, gg;
  RhoCandList epem_tr, pippim_tr, gg_tr;

  RhoCandList pi0, pi0_true; // TODO: refactor pi0_true -> pi0_pm
  RhoCandList pi0_ana, pi0_pm_ana;

  static const int npi0ana = 16;
  //std::vector<RhoCandList> pi0_ana_, pi0_pm_ana_;
  RhoCandList pi0_ana_[npi0ana], pi0_pm_ana_[npi0ana];
  RhoCandList pi0nearest, pi0_btb, pi0_cts;

  RhoCandList epem_mcut, pippim_mcut;
  RhoCandList jpsi, jpsi_true, jpsi_ana, jpsi_pm_ana; // TODO: refactor jpsi_true -> jpsi_pm

  RhoCandList all_ana, all_pm_ana;

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

  TH1F *h_m_pippim;
  TH1F *h_m_pippim_tr;

  TH1F *h_4c_chi2[3][2];
  TH1F *h_4c_prob[3][2];
  TH1F *h_4c_m[3][2];
  TH2F *h_4c_prob_vs_m[3][2];

  TH1F* h_num_g;
  TH1F* h_num_epm;
  TH1F* h_num_pipm;
  TH2F* h_mom_the_pipm;
  TH1F* h_num_g_tr;
  TH1F* h_num_epm_tr;
  TH1F* h_num_pipm_tr;
  TH2F* h_mom_the_pipm_tr;

  static const int nhist = 6;
  //TH1F *h_dth_gg_epair;
  //TH1F *h_mass_gg_epair;
  TH2F *h_dth_gg_epair_vs_mass_gg[nhist];
  TH2F *h_mass_gg_epair_vs_mass_gg[nhist];
  TH2F *h_dth_vs_mass_gg_epair[nhist];
  TH2F *h_dth_vs_mass_gg_epair_btb[nhist];  // most back-to-back
  TH2F *h_dth_vs_mass_gg_epair_cts[nhist];  // closest-to-s
  TH1F *h_m_gg_btb[nhist];  // most back-to-back
  TH1F *h_m_gg_cts[nhist];  // closest-to-s
  TH1F *h_m_epem_btb[nhist];
  TH1F *h_m_epem_cts[nhist];

  // static const int npi0ana = 3;
  // static const int nhist = 6;
  enum {all=0, tr=1, pm=2, pm_mc=3, ana=4, pm_ana=5};

  TH1F *h_m_gg_pi0ana[npi0ana][nhist];
  TH1F *h_oa_gg[nhist];
  TH1F *h_m_gg[nhist];
  TH1F *h_e_g[nhist];

  TH1F *h_m_epem[nhist];
  TH1F *h_mom_epm[nhist];
  TH2F *h_mom_the_epm[nhist];

  ClassDef(AnaTda, 1);
};

#endif /* AnaTda_H */

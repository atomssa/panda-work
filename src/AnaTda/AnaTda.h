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
  AnaTda(const int&, const int&, const int&);

  // ** Destructor
  ~AnaTda();

  // ** Virtual method Init
  virtual InitStatus Init();

  // ** Virtual method Exec
  virtual void Exec(Option_t* opt);

  virtual void Finish();

  void set_verbosity(const int &v) {verb = v;}

  protected:
  PndAnalysis *fAnalysis;             // *** the PndAnalysis object

  int nevt;
  bool bg_mc;
  int iplab;
  double plab[3];
  double p_antip;

  int verb;

  bool fBremCorr;

  int pdg_jpsi;
  int pdg_pi0;
  int pdg_pip;
  int pdg_pim;
  int pdg_em;
  int pdg_ep;
  double m0_jpsi;
  double m0_pi0;

  double up[3][3]; // pi0 OA vs Eavg cut upper limit parameters
  double lw[3][3]; // pi0 OA vs Eavg cut lower limit parameters
  double pi0mcut_min;
  double pi0mcut_max;
  double jpsi_mcut_min;
  double jpsi_mcut_max;
  double etot_min;
  double etot_max;
  double dth_min;
  double dth_max;

  RhoCandList mcList;
  RhoCandList pip, pim, ep, em, g;
  RhoCandList pip_tr, pim_tr, ep_tr, em_tr, g_tr;
  RhoCandList epem, pippim, gg;
  RhoCandList epem_tr, pippim_tr, gg_tr;
  RhoCandList pi0, pi0_true; // TODO: refactor pi0_true -> pi0_pm
  RhoCandList pi0_ana, pi0_pm_ana;
  static const int npi0ana = 2;
  std::vector<RhoCandList> pi0_ana_, pi0_pm_ana_;
  RhoCandList pi0nearest, pi0_btb, pi0_cts;
  RhoCandList epem_mcut, pippim_mcut;
  RhoCandList jpsi, jpsi_true, jpsi_ana, jpsi_pm_ana; // TODO: refactor jpsi_true -> jpsi_pm
  RhoCandList jpsi_mconst;
  RhoCandList pi0jpsi_ana, pi0jpsi_pm_ana;
  RhoCandList epem_mcut_pi0_btb, epem_mcut_pi0_cts;
  RhoCandList pippim_mcut_pi0_btb, pippim_mcut_pi0_cts;
  RhoCandList epem_pi0nearest;
  RhoCandList pippim_pi0nearest;

  TLorentzVector ini;
  TVector3 boost_to_cm;
  double sqrt_s;

  RhoMassParticleSelector *jpsiMassSel;
  RhoMassParticleSelector *pi0MassSel;

  void calc_kin(RhoCandidate*, RhoCandidate *, double &, double &, double &, double &, double &, double &);
  void calc_kin_from_daughters(RhoCandidate*, RhoCandidate *, RhoCandidate*, RhoCandidate *, double &, double &, double &, double &);
  void def_hists();
  void initial_state();
  void set_selectors();

  void def_eff_hists();
  void def_resid_hists();
  void def_tutorial_hists();
  void def_manual_kin_fit_hists(const int&);
  void def_pair_hists();
  void def_single_hists();
  void def_full_sys_hists();
  void def_kin_fit_hists(const int& type, const int&);
  void def_kin_fit_hists();
  void def_gamma_from_pi0_hists();
  void def_elecs_from_jpsi_hists();
  void write_tut_hists();
  void write_manual_kin_fit_hists(const int&);
  void write_kin_fit_hists();
  void write_full_sys_hists();
  void write_eff_hists();

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

  void truth_match(RhoCandList&, RhoCandList& dest, const int &);
  void truth_match_singles();
  void fill_single_dists_tr();

  // Pull out candidates
  void make_pair_lists();
  void fill_pair_mass(RhoCandList&, TH1F*);
  void fill_pair_oa(RhoCandList&, TH1F*);
  void fill_pair_oa(RhoCandList&, const int &);
  void fill_pair_oa_mc(RhoCandList&, const int &);
  void fill_pair_oa(RhoCandidate*, RhoCandidate*, const int &);
  void fill_pair_dists();

  void truth_match_residuals();

  double dist_pi0_pair_match(RhoCandidate*);
  double dist_jpsi_pair_match(RhoCandidate*);
  double dist_photon_match(RhoCandidate*, RhoCandidate*); // This needs special treatment due to lack of matches in BG MC
  double dist_chpi_match(RhoCandidate*, RhoCandidate*);

  double primary_match_pair(const int&, RhoCandidate*, int&, int&, int&);

  void primary_match(const int&, RhoCandidate*, RhoCandidate*); // TODO - Absorb in primary_match_single
  bool primary_match_single(const int&, RhoCandidate*, double&);

  void find_primary_gg();
  void fill_gamma_from_pi0s();
  bool oa_vs_avg_cut(const double&, const double &);
  void pi0_analysis_cut();
  void pi0_analysis_cut(RhoCandList&, RhoCandList&, bool);
  void fill_pi0_analysis_hists(RhoCandList&, const int&);
  void fill_pi0_analysis_hists(RhoCandList&, const int&, const int&);

  void find_primary_epem();
  void find_primary_pippim();
  void fill_elecs_from_jpsi();
  void jpsi_analysis_cut();
  void jpsi_analysis_cut(RhoCandList&, RhoCandList&);
  void fill_jpsi_analysis_hists(RhoCandList&, const int&);

  void pi0jpsi_kin_fit(RhoCandList& org_gg, RhoCandList& org_epem); // second try at kinematic fitting
  void pi0jpsi_efficiency(RhoCandList& org_gg, RhoCandList& org_epem, const int& tt); // second try at kinematic fitting

  bool passes_kin_cut(const double &, const double &);

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

  void pi0jpsi_kinematics(RhoCandList&, RhoCandList&, const int&);

  void pi0_kinematic_selection(RhoCandList&, RhoCandList&, const int&);
  void pi0jpsi_true_kinematics(RhoCandList&, RhoCandList&);
  void pdgm_nearest_pi0s();

  // *** create some histograms
  TH2F *h_resid_phth[4];
  TH2F *h_resid_pip_phth[4];
  TH1F *h_resid_pip_mom[4];
  TH2F *h_resid_pim_phth[4];
  TH1F *h_resid_pim_mom[4];
  TH2F *h_resid_ep_phth[4];
  TH1F *h_resid_ep_mom[4];
  TH2F *h_resid_em_phth[4];
  TH1F *h_resid_em_mom[4];
  TH2F *h_resid_pip_elec_hyp_phth[4];
  TH1F *h_resid_pip_elec_hyp_mom[4];
  TH2F *h_resid_pim_elec_hyp_phth[4];
  TH1F *h_resid_pim_elec_hyp_mom[4];

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

  static const int nhist = 7;
  TH2F *h_dth_gg_epair_vs_mass_gg[nhist];
  TH2F *h_mass_gg_epair_vs_mass_gg[nhist];
  TH2F *h_dth_vs_mass_gg_epair[nhist];
  TH2F *h_dth_vs_dph_gg_epair[nhist];
  TH2F *h_dth_vs_mass_gg_epair_btb[nhist];  // most back-to-back
  TH2F *h_dth_vs_mass_gg_epair_cts[nhist];  // closest-to-s
  TH1F *h_m_gg_btb[nhist];  // most back-to-back
  TH1F *h_m_gg_cts[nhist];  // closest-to-s
  TH1F *h_m_epem_btb[nhist];
  TH1F *h_m_epem_cts[nhist];

  TH1F *h_m_epem_bef_mass_fit;
  TH1F *h_m_epem_aft_mass_fit;

  // static const int npi0ana = 3;
  // static const int nhist = 6;
  enum {all=0, tr=1, pm=2, pm_mc=3, ana=4, pm_ana=5, mconst=6};

  TH1F *h_m_gg_pi0ana[npi0ana][nhist];
  TH1F *h_oa_gg[nhist];
  TH2F *h_oa_gg_vs_min_e_g[nhist];
  TH2F *h_oa_gg_vs_avg_e_g[nhist];
  TH2F *h_oa_gg_vs_asym_e_g[nhist];

  TH1F *h_m_gg[nhist];
  TH1F *h_e_g[nhist];

  TH1F *h_m_epem[nhist];
  TH1F *h_mom_epm[nhist];
  TH2F *h_mom_the_epm[nhist];

  enum {eff_ref=0, eff_pi0sel=1, eff_jpsisel=2, eff_kin=3, eff_excl=4, eff_const=5};
  static const int neffhist = 6;
  TH1F *h_eff_thpi0[neffhist];
  TH1F *h_eff_thep[neffhist];

  // Number of pairs after each analysis cut step.
  // This will measure how effective the cuts are at
  // reducing combinatoric multiplicity
  TH1F *h_num_gg[neffhist];
  TH1F *h_num_epem[neffhist];

  ClassDef(AnaTda, 1);
};

#endif /* AnaTda_H */

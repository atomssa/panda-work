#ifndef AnaTdav2_H
#define AnaTdav2_H

#include "FairTask.h"
#include "PndAnalysis.h"

#include "RhoCandidate.h"

#include <vector>
#include <string>

class RhoCandList;
class TH1F;
class TH2F;
class TF1;
class TEfficiency;
class TVector3;
class TFile;
class TLorentzVector;
class PndPidProbability;
class TClonesArray;
class FairRootManager;

class AnaTdav2 : public FairTask{

  typedef Double_t (PndPidProbability::*prob_func)(PndPidProbability*) const;

 public:

  AnaTdav2(const int&, const int&, const int&, const int&);
  ~AnaTdav2();

  virtual InitStatus Init();

  virtual void Exec(Option_t* opt);

  virtual void FinishTask();
  virtual void FinishEvent() {}

  void set_verbosity(int v) {verb = v;}

  void set_eff_file_name(std::string a) {eff_file_name= a;}
  void set_eff_hist_name(std::string a, bool rad) {eff_hist_name= a; eff_hist_rad= rad; }

  void set_pi_eff_file_name(std::string a) {pi_eff_file_name= a;}
  void set_pi_eff_hist_name(std::string a, bool rad) {pi_eff_hist_name= a; pi_eff_hist_rad= rad; }

  void set_eid_prob_min(double pmin) {eid_prob_min = pmin;}

 private:
  void fill_lists();
  void cleanup_lists() { for (int ii = 0; ii < rcl.size(); ++ii) rcl[ii].Cleanup(); }
  TClonesArray* init_tca(TString);
  void init_tcas();
  void init_hists();
  void beam_cond();

  void nocut_ref();
  void ep_uniq();
  void ep_all();
  void pi0_sel();
  void ep_pi0_asso();
  void ep_pi0_asso_all();
  void kin_excl();
  void kin_excl_all();
  void kin_fit();
  void kin_fit_4c();
  void fill_bins_excl();
  void fill_bins_kinc();
  void fill_bins_kinc_bg();
  void fill_bins_ngcut();
  void fill_bins(RhoCandList&, RhoCandList&);
  void write_hists();
  void print_mc_list();
  bool bayes_pid(RhoCandidate*);
  bool check_eid(RhoCandidate*);
  void eid_filter(RhoCandList&, RhoCandList&);
  double get_comb_prob(prob_func func);
  double eff_weight(const TVector3 &mom);

  //double _pi_eff_func(double *x, double *p);

  TH2F* smooth_hist2d(TH2F* , int);
  TEfficiency* smooth_eff2d(TEfficiency *,int);

  TH1F* smooth_hist1d(TH1F*);
  TEfficiency* smooth_eff1d(TEfficiency *);
  TEfficiency* rebin2d(TEfficiency *, int);

  double dist_chpi_match(RhoCandidate*, RhoCandidate*);
  double dist_photon_match(RhoCandidate*, RhoCandidate*);
  void charged_pion_filter(RhoCandList&, RhoCandList&, RhoCandList&, RhoCandList&, RhoCandList&, RhoCandList&);
  void mctruth_match_jpsi(RhoCandList&, RhoCandList&);
  void mctruth_match_pi0(RhoCandList&);
  int mct_uid_e, mct_uid_p;
  int mct_uid_g1, mct_uid_g2;
  int mc_elec, mc_posit;
  int mc_g1, mc_g2;

 private:

  int verb;
  int nevt;
  bool brem_corr;

  int mc_type;
  bool eid_param_method;

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

  int iplab;
  double plab[3];
  double mom_antip;
  TVector3 boost_to_cm;
  TVector3 boost_to_lab;
  TLorentzVector p4pbar, p4targ, p4sys;
  double event_t, event_u;
  double event_pi0costh_cm;
  double event_pi0theta_cm;
  double event_pi0theta_lab;
  double event_epthe_jpsi;
  double event_epcosth_jpsi;
  double m_pip_wt, m_pim_wt, m_evt_wt;
  double m_epcth_wt0, m_epcth_wt1;
  double tmin[3]; //={-0.443789, -1.0, -1.0};
  double tmid[3]; //={-0.443789, -1.0, -1.0};
  double tmax[3]; //={0.616486, 0.457248, 0.31538};
  double nevt_sim[7][3]; // number of simulated events for x-sect normalization
  double nevt_xsect[7][3]; // number of simulated events for x-sect normalization

  // Efficiency parametrizations
  TFile *eff_file;
  std::string eff_file_name;
  std::string eff_hist_name;
  TH2F* heff_epm;
  bool eff_hist_rad;

  // pion misid parametrization
  TFile *pi_eff_file;
  std::string pi_eff_file_name;
  std::string pi_eff_hist_name;
  bool pi_eff_hist_rad;
  TEfficiency* pi_eff;
  TF1* pi_eff_func;

  double eid_prob_min;

  static const int nstep = 9;
  TH1F* hnevt;
  TH1F* hwt;
  TH1F* hng;
  TH1F* hng20mev;
  TH1F* hnch;
  TH1F* hmep_mct[nstep], *hmep_non_mct[nstep], *hthe_ep_mct[nstep], *hthe_ep_mct_fwd[nstep], *hthe_ep_mct_bwd[nstep];
  TH1F* hmgg_mct[nstep], *hmgg_non_mct[nstep], *hthe_gg_mct[nstep], *hthe_gg_mct_fwd[nstep], *hthe_gg_mct_bwd[nstep], *hoa_gg_mct[nstep];
  TH1F* hmep[nstep], *hmep_valid[nstep];
  TH1F* hmgg[nstep], *hmgg_valid[nstep];
  TH1F* hnep[nstep];
  TH1F* hngg[nstep];
  TH1F* hnpi0jpsi[nstep];
  TH1F* hpi0cost_cm, *hpi0th;
  TH1F* hpi0cost_cm_mcut, *hpi0th_mcut;

  //static const int nbinth = 12;
  std::vector<double> tu_binning;
  std::vector<double> pi0th_binning;
  std::vector<double> pi0cost_cm_binning;
  std::vector<double> tu_binning_2d;
  std::vector<double> costh_binning_2d;

  std::vector<TH1F*> hmep_pi0cost_cm;
  std::vector<TH1F*> hmep_pi0th;
  std::vector<TH1F*> hmept;
  std::vector<TH1F*> hmepu;

  std::vector<TH1F*> hmeptcth,hmeptcth0,hmeptcth1;
  std::vector<TH1F*> f_hmeptcth,f_hmeptcth0,f_hmeptcth1;

  std::vector<TH1F*> hmepucth,hmepucth0,hmepucth1;
  std::vector<TH1F*> f_hmepucth,f_hmepucth0,f_hmepucth1;

  TH1F* hmmiss, *hmmiss2;
  TH1F* hmmiss_jpsi, *hmmiss2_jpsi;
  TH1F* hmtot;
  TH2F* hcmoa;
  TH1F* htrecgg, *hurecgg, *htrecep, *hurecep;
  TH1F* htrecgg_mcut, *hurecgg_mcut, *htrecep_mcut, *hurecep_mcut;
  TH1F* httrumc, *hutrumc;
  TH1F* httrumc_vb, *hutrumc_vb;
  TH1F* htrupi0thcm, *htrupi0thlab, *htrupi0costhcm;
  TH1F* htrupi0thcm_mcut, *htrupi0thlab_mcut, *htrupi0costhcm_mcut;
  TH1F* htrupi0thcm_tcut, *htrupi0thlab_tcut, *htrupi0costhcm_tcut;
  TH1F* htrupi0thcm_tcut_mcut, *htrupi0thlab_tcut_mcut, *htrupi0costhcm_tcut_mcut;
  TH2F* htrupi0thcm_vs_m, *htrupi0thlab_vs_m, *htrupi0costhcm_vs_m;
  TH2F* htrupi0thcm_mcut_vs_m, *htrupi0thlab_mcut_vs_m, *htrupi0costhcm_mcut_vs_m;
  TH1F* htresgg, *huresgg, *htresep, *huresep;

  TH2F* hepcosth_jpsi_vs_epthlab_mc_all, *hepcosth_jpsi_vs_emthlab_mc_all;
  TH1F* hepcosth_jpsi_mc_all, *hepcosth_jpsi_mc_all_wt0, *hepcosth_jpsi_mc_all_wt1;
  TH1F* hepcosth_jpsi_mc[4];
  TH2F* hepcosth_jpsi_vs_epthlab_mc[4], *hepcosth_jpsi_vs_emthlab_mc[4];
  TH1F* hepcosth_jpsi_rec_all;
  TH2F* hepcosth_jpsi_vs_epthlab_rec_all, *hepcosth_jpsi_vs_emthlab_rec_all;
  TH1F* hepcosth_jpsi_rec[4];
  TH2F* hepcosth_jpsi_vs_epthlab_rec[4], *hepcosth_jpsi_vs_emthlab_rec[4];
  TH1F* f_hepcosth_jpsi_rec_all;
  TH2F* f_hepcosth_jpsi_vs_epthlab_rec_all, *f_hepcosth_jpsi_vs_emthlab_rec_all;
  TH1F* f_hepcosth_jpsi_rec[4];
  TH2F* f_hepcosth_jpsi_vs_epthlab_rec[4], *f_hepcosth_jpsi_vs_emthlab_rec[4];

  TH1F *hmep_mconst;
  TH1F *hmtot_mconst;
  TH2F *hcmoa_mconst;
  TH1F *hmtot_mconst_cut;
  TH2F *hcmoa_mconst_cut;
  TH1F *hpi0jpsi_chi24c;
  TH1F *hpi0jpsi_chi24c_c;
  TH1F *hpi0jpsi_prob4c;
  TH1F *hpi0jpsi_pull4c;
  TH1F *hpi0jpsi_chi2diff4c;
  TH2F *hpi0jpsi_chi24c_vs_mtot_r;
  TH2F *hpi0jpsi_chi24c_vs_cm_dth_r;
  TH2F *hpi0jpsi_chi24c_vs_cm_dph_r;
  TH2F *hpi0jpsi_chi24c_vs_mtot_f;
  TH2F *hpi0jpsi_chi24c_vs_cm_dth_f;
  TH2F *hpi0jpsi_chi24c_vs_cm_dph_f;

  TH1F *hpi0pi0jpsi_chi24c;
  TH1F *hpi0pi0jpsi_chi24c_c;
  TH2F *hpi0vs2pi0_chi24c;
  TH2F *hpi0vs2pi0_chi24c_c;

  TH1F *hmom_pull_ep_r;
  TH1F *hmom_pull_ep_f;
  TH1F *hmom_pull_em_r;
  TH1F *hmom_pull_em_f;
  TH1F *hpx_pull_ep_r;
  TH1F *hpx_pull_ep_f;
  TH1F *hpx_pull_em_r;
  TH1F *hpx_pull_em_f;
  TH1F *hpy_pull_ep_r;
  TH1F *hpy_pull_ep_f;
  TH1F *hpy_pull_em_r;
  TH1F *hpy_pull_em_f;
  TH1F *hpz_pull_ep_r;
  TH1F *hpz_pull_ep_f;
  TH1F *hpz_pull_em_r;
  TH1F *hpz_pull_em_f;

  //static const int nrcl = 7;
  enum {e=0, p, g, /* Single Track Lists, no PID: e=elec p=posit, g=gamma*/
	pip, pim, /* List of pion tracks, loaded with pion hypothesis hopefully ordered the same way as electrons*/
	ie, ip, /* Single Track Lists with PID ie=elec, ip=posit */
	gg, gg_sel,/* GG pairs gg=all gg_sel=after selection cuts */
	ep, /* all e-p pairs */
	iep, /* pid'ed e-p pairs, no other cond*/
	iep_uniq, /* pid'ed e-p pairs, require uniquness, (no other charged track) */
	iep_all, /* pid'ed e-p pairs, not requireing uniquness, */
	iep_asso, /* pid'ed e-p pairs, require associated gg pair satisfying selection cuts */
	iep_asso_all, /* pid'ed e-p pairs, require associated gg pair satisfying selection cuts */
	iep_excl, /* pid'ed e-p pairs, require exclusivity of pi0-e-p */
	gg_excl, /* gg pairs, require exclusivity of pi0-e-p */
	iep_kinc, /* iep_excl + cut on singal hypothesis kinematic fit chi2 */
	gg_kinc, /* gg_excl + cut on singal hypothesis kinematic fit chi2 */
	iep_kinc_bg, /* iep_kinc +  accept only if signal hypothesis fits better than bg hypothesis */
	gg_kinc_bg, /* gg_kinc + accept only if signal hypothesis fits better than bg hypothesis */
	iep_ngcut, /* iep_kinc_bg + reject events with >4 20MeV neutral candidates */
	gg_ngcut, /* gg_kinc_bg +  reject events with >4 20MeV neutral candidates */
	nrcl /*number of entries */
  };
  std::vector<RhoCandList> rcl;
  RhoCandList mcList; // this one special..

  // Utilities
  bool check_mct_jpsi(RhoCandidate*);
  bool check_mct_pi0(RhoCandidate*);
  double oa(RhoCandidate*);
  double the_bwd(RhoCandidate*);
  double the_fwd(RhoCandidate*);
  bool oa_vs_avg_cut(const double&, const double &);
  double oa(RhoCandidate* c1, RhoCandidate* c2) { return c1->P3().Angle(c2->P3()); }
  double m(RhoCandidate* c1, RhoCandidate* c2) { return (c1->P4() + c2->P4()).M(); }
  double dth(RhoCandList* c1, RhoCandList* c2);
  double t_gg(RhoCandidate *_gg) { return (_gg->P4()-p4pbar).M2(); }
  double t_ep(RhoCandidate *_ep) { return (_ep->P4()-p4targ).M2(); }
  double u_gg(RhoCandidate *_gg) { return (_gg->P4()-p4targ).M2(); }
  double u_ep(RhoCandidate *_ep) { return (_ep->P4()-p4pbar).M2(); }
  double err_mom_sq(RhoCandidate*);
  double mom_pull_r(RhoCandidate*, RhoCandidate*);
  double mom_pull_f(RhoCandidate*, RhoCandidate*);
  double px_pull_r(RhoCandidate*, RhoCandidate*);
  double px_pull_f(RhoCandidate*, RhoCandidate*);
  double py_pull_r(RhoCandidate*, RhoCandidate*);
  double py_pull_f(RhoCandidate*, RhoCandidate*);
  double pz_pull_r(RhoCandidate*, RhoCandidate*);
  double pz_pull_f(RhoCandidate*, RhoCandidate*);

  TLorentzVector boost_transf(const TLorentzVector&, const TVector3&);
  double cost_b(const TLorentzVector&, const TVector3&);
  double the_b(const TLorentzVector&, const TVector3&);
  TLorentzVector get_p4ep(RhoCandidate*);
  TLorentzVector get_p4em(RhoCandidate*);
  int comb_bins(int,int,int);
  int find_bin(double val, const std::vector<double>& binning);
  //int t_bin(RhoCandidate*);
  //int u_bin(RhoCandidate*);
  //int pi0th_bin(RhoCandidate*);
  //int pi0cost_cm_bin(RhoCandidate*);
  double pi0cost_cm(RhoCandidate*);
  double pi0theta_cm(RhoCandidate*);
  //void fill_ang_dist(RhoCandidate*);
  double dph_cm(TLorentzVector, TLorentzVector);
  double dth_cm(TLorentzVector, TLorentzVector);
  void dth_dph_cm(TLorentzVector, TLorentzVector, double&, double&);
  void dth_dph_cm(RhoCandidate*, RhoCandidate *, double &, double &);
  void fill_mtot(RhoCandList&, RhoCandList&, TH1F*);
  void fill_mmiss(RhoCandList&, RhoCandList&, TH1F*, TH1F*);
  void fill_mmiss_jpsi(RhoCandList&, TH1F*, TH1F*);
  void fill_dth_dph_cm(RhoCandList&, RhoCandList&, TH2F*);
  void fill_pair_mass(RhoCandList&, TH1F*);
  void fill_pair_mass(RhoCandList&, TH1F*, double);
  void fill_count_hists(int, int, int);
  void fill_mctruth(RhoCandList&, RhoCandList&, int);
  void print_indices();
  void print_binning(const std::vector<double>&, const char*);
  bool calc_true_tu();
  void calc_evt_wt();

  bool is_dpm();
  bool is_evt_gen();

  ClassDef(AnaTdav2,1);

  // Analysis cut setters
 private:
  // Analysis cut flags and values
  bool apply_pi0evsoa_cut;
  double up[3][3]; // pi0 OA vs Eavg cut upper limit parameters
  double lw[3][3]; // pi0 OA vs Eavg cut lower limit parameters

  bool apply_pi0m_cut;
  double pi0m_cut_min;
  double pi0m_cut_max;

  bool apply_mtot_cut;
  double mtot_cut_min;
  double mtot_cut_max;

  bool apply_dth_dph_cut;
  double dth_sigma;
  double dph_sigma;
  double dth_dph_cm_cut_max;

  bool require_exclusivity;

  double jpsi_m_3sig_min;
  double jpsi_m_3sig_max;

  double chi2_cut[3];

  int ng20mev;

 public:
  void do_apply_pi0evsoa_cut(bool a) {apply_pi0evsoa_cut = a;}
  void do_apply_pi0m_cut(bool a) {apply_pi0m_cut = a;}
  void do_apply_mtot_cut(bool a) {apply_mtot_cut = a;}
  void do_apply_dth_dph_cut(bool a) {apply_dth_dph_cut = a;}
  void do_require_exclusivity(bool a) {require_exclusivity = a; }

  void set_pi0m_cut(double min=0.11, double max=0.16) { pi0m_cut_min = min; pi0m_cut_max = max; }
  void set_dph_sigma(double s=0.4 /* Sigmas */) { dph_sigma = s; }
  void set_dth_sigma(double s=0.4 /* Sigmas */) { dth_sigma = s; }
  void set_dth_dph_cut(double c=3.0 /* Sigmas */) { dth_dph_cm_cut_max = c; }
  void set_mtot_cut(double min=3.3, double max=3.7) { mtot_cut_min = min; mtot_cut_max = max; }

};

#endif /*AnaTdav2_H*/

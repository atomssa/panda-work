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
class TVector3;
class TFile;
class TLorentzVector;

class AnaTdav2 : public FairTask{

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

 private:
  void fill_lists();
  void cleanup_lists() { for (int ii = 0; ii < rcl.size(); ++ii) rcl[ii].Cleanup(); }
  void init_hists();
  void beam_cond();

  void nocut_ref();
  void ep_uniq();
  void pi0_sel();
  void ep_pi0_asso();
  void kin_excl();
  void kin_fit();
  void fill_bins();
  void write_hists();
  void print_mc_list();
  bool bayes_pid(RhoCandidate*);
  bool check_eid(RhoCandidate*);
  void eid_filter(RhoCandList&, RhoCandList&);

  double dist_chpi_match(RhoCandidate*, RhoCandidate*);
  void charged_pion_filter(RhoCandList&, RhoCandList&, RhoCandList&, RhoCandList&, RhoCandList&, RhoCandList&);

 private:

  int verb;
  int nevt;
  bool brem_corr;

  bool bg_mc;
  bool eid_param_method;

  PndAnalysis *fAna;

  int iplab;
  double plab[3];
  double mom_antip;
  TVector3 boost_to_cm;
  TVector3 boost_to_lab;
  TLorentzVector p4pbar, p4targ;
  double event_t, event_u;
  double tmin[3]; //={-0.443789, -1.0, -1.0};
  double tmax[3]; //={0.616486, 0.457248, 0.31538};

  // Efficiency parametrizations
  TFile *eff_file;
  std::string eff_file_name;
  std::string eff_hist_name;
  bool eff_hist_rad;
  TH2F* heff_epm;

  static const int nstep = 6;
  TH1F* hmep[nstep];
  TH1F* hnep[nstep];
  TH1F* hngg[nstep];
  TH1F* hpi0cost_cm;
  TH1F* hpi0th;
  static const int nbinth = 12;
  double tu_binning[nbinth+1];
  double pi0th_binning[nbinth+1];
  double pi0cost_cm_binning[nbinth+1];
  TH1F* hmep_pi0cost_cm[nbinth];
  TH1F* hmep_pi0th[nbinth];
  TH1F* hmept[nbinth];
  TH1F* hmepu[nbinth];
  TH1F* hmtot;
  TH2F* hcmoa;
  TH1F* htrecgg, *hurecgg, *htrecep, *hurecep;
  TH1F* httrumc, *hutrumc;
  TH1F* htrupi0thcm, *htrupi0thlab, *htrupi0costhcm;
  TH1F* htresgg, *huresgg, *htresep, *huresep;

  TH1F *hmep_mconst;
  TH1F *hmtot_mconst;
  TH2F *hcmoa_mconst;
  TH1F *hmtot_mconst_cut;
  TH2F *hcmoa_mconst_cut;

  //static const int nrcl = 7;
  enum {e=0, p, g, /* Single Track Lists, no PID: e=elec p=posit, g=gamma*/
	pip, pim, /* List of pion tracks, loaded with pion hypothesis hopefully ordered the same way as electrons*/
	ie, ip, /* Single Track Lists with PID ie=elec, ip=posit */
	gg, gg_sel,/* GG pairs gg=all gg_sel=after selection cuts */
	ep, /* all e-p pairs */
	iep, /* pid'ed e-p pairs, no other cond*/
	iep_uniq, /* pid'ed e-p pairs, require uniquness, (no other charged track) */
	iep_asso, /* pid'ed e-p pairs, require associated gg pair satisfying selection cuts */
	iep_excl, /* pid'ed e-p pairs, require exclusivity of pi0-e-p */
	gg_excl, /* gg pairs, require exclusivity of pi0-e-p */
	nrcl /*number of entries */
  };
  std::vector<RhoCandList> rcl;
  RhoCandList mcList; // this one special..

  // Utilities
  bool oa_vs_avg_cut(const double&, const double &);
  double oa(RhoCandidate* c1, RhoCandidate* c2) { return c1->P3().Angle(c2->P3()); }
  double m(RhoCandidate* c1, RhoCandidate* c2) { return (c1->P4() + c2->P4()).M(); }
  double dth(RhoCandList* c1, RhoCandList* c2);
  double t_gg(RhoCandidate *_gg) { return (_gg->P4()-p4pbar).M2(); }
  double t_ep(RhoCandidate *_ep) { return (_ep->P4()-p4targ).M2(); }
  double u_gg(RhoCandidate *_gg) { return (_gg->P4()-p4targ).M2(); }
  double u_ep(RhoCandidate *_ep) { return (_ep->P4()-p4pbar).M2(); }

  int find_bin(double val, double *binning);
  //int t_bin(RhoCandidate*);
  //int u_bin(RhoCandidate*);
  //int pi0th_bin(RhoCandidate*);
  //int pi0cost_cm_bin(RhoCandidate*);
  double pi0cost_cm(RhoCandidate*);
  double pi0theta_cm(RhoCandidate*);
  //void fill_ang_dist(RhoCandidate*);
  void dth_dph_cm(RhoCandidate*, RhoCandidate *, double &, double &);
  void fill_mtot(RhoCandList&, RhoCandList&, TH1F*);
  void fill_dth_dph_cm(RhoCandList&, RhoCandList&, TH2F*);
  void fill_pair_mass(RhoCandList&, TH1F*);
  void fill_count_hists(int, int, int);
  void print_indices();
  void print_binning(double *, const char*);
  bool calc_true_tu();

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

#ifndef AnaTda_H
#define AnaTda_H 1

#include "FairTask.h"
#include "TLorentzVector.h"

#include <vector>

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

class AnaTda : public FairTask
{

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

  TH1F *h_4c_chi2_epempi0;
  TH1F *h_4c_prob_epempi0;
  TH1F *h_4c_m_epem;
  TH2F *h_4c_prob_vs_m_epem;
  TH1F *h_4c_chi2_pippimpi0;
  TH1F *h_4c_prob_pippimpi0;
  TH1F *h_4c_m_pippim;
  TH2F *h_4c_prob_vs_m_pippim;

  TH1F *h_4c_chi2_btb_epempi0;
  TH1F *h_4c_prob_btb_epempi0;
  TH1F *h_4c_m_btb_epem;
  TH2F *h_4c_prob_vs_m_btb_epem;
  TH1F *h_4c_chi2_btb_pippimpi0;
  TH1F *h_4c_prob_btb_pippimpi0;
  TH1F *h_4c_m_btb_pippim;
  TH2F *h_4c_prob_vs_m_btb_pippim;

  TH1F *h_4c_chi2_cts_epempi0;
  TH1F *h_4c_prob_cts_epempi0;
  TH1F *h_4c_m_cts_epem;
  TH2F *h_4c_prob_vs_m_cts_epem;
  TH1F *h_4c_chi2_cts_pippimpi0;
  TH1F *h_4c_prob_cts_pippimpi0;
  TH1F *h_4c_m_cts_pippim;
  TH2F *h_4c_prob_vs_m_cts_pippim;

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
  TH2F *h_oa_vs_mass_gg_epair_btb; // most back-to-back
  TH2F *h_oa_vs_mass_gg_epair_cts; // closest-to-s
  TH1F *h_m_gg_btb; // most back-to-back
  TH1F *h_m_gg_cts; // closest-to-s

  ClassDef(AnaTda,1);

};

#endif /* AnaTda_H */

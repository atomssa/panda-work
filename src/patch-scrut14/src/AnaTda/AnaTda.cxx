// The header file
#include "AnaTda.h"

// C++ headers
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

// FAIR headers
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairRun.h"
#include "FairRadLenPoint.h"

#include "PndAnalysis.h"
#include "RhoCandList.h"
#include "RhoMassParticleSelector.h"
#include "RhoCandidate.h"
#include "PndKinFitter.h"
// ROOT headers
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"

#include "filler.h"

// other headers
// #include "PndEmcBump.h"
// #include "PndEmcCluster.h"
// #include "PndEmcBump.h"
// #include "PndEmcRecoHit.h"
// #include "PndEmcMapper.h"
#include "PndMCTrack.h"

using std::cout;
using std::endl;
using std::vector;

AnaTda::~AnaTda() { }

AnaTda::AnaTda(const int& _iplab, const int& itype, const int& brem)
  :FairTask("Radiation Length Profiler"),
   nevt(0),
   bg_mc(itype==0),
   verb(0),
   iplab((_iplab>=0&&_iplab<3)?_iplab:0),
   plab{5.513,8.,12.},
   p_antip(plab[iplab]),
   fBremCorr(brem!=0),
   pdg_jpsi(443),
   pdg_pi0(111),
   pdg_pip(211),
   pdg_pim(-pdg_pip),
   pdg_em(11),
   pdg_ep(-pdg_em),
   m0_jpsi(0.),
   m0_pi0(0.),
   lw{{0.11,-0.05,+0.00},{0.12,-0.03,-0.03},{0.15,-0.06,-0.07}},
   up{{0.14,0.21,0.07},{0.12,0.15,0.15},{0.11,0.10,0.20}},
   pi0mcut_min(0.11),
   pi0mcut_max(0.16),
   jpsi_mcut_min(2.5), // (2.96)
   jpsi_mcut_max(3.7),  // (3.22)
   etot_min(3.4),
   etot_max(3.6),
   dth_min(3.0),
   dth_max(3.3),
   mcList(),
   pip(), pim(), ep(), em(), g(),
   pip_tr(), pim_tr(), ep_tr(), em_tr(), g_tr(),
   epem(), pippim(), gg(),
   epem_tr(), pippim_tr(), gg_tr(),
   pi0(), pi0_true(), // TODO: refactor pi0_true -> pi0_pm
   pi0_ana(), pi0_pm_ana(),
   pi0_ana_(npi0ana),pi0_pm_ana_(npi0ana),
   pi0nearest(), pi0_btb(), pi0_cts(),
   epem_mcut(), pippim_mcut(),
   jpsi(), jpsi_true(), jpsi_ana(), jpsi_pm_ana(), // TODO: refactor jpsi_true -> jpsi_pm
   jpsi_mconst(),
   pi0jpsi_ana(), pi0jpsi_pm_ana(),
   epem_mcut_pi0_btb(), epem_mcut_pi0_cts(),
   pippim_mcut_pi0_btb(), pippim_mcut_pi0_cts(),
   epem_pi0nearest(),
   pippim_pi0nearest() {
     assert(iplab>=0&&iplab<=2);
   }

InitStatus AnaTda::Init() {
  fAnalysis = new PndAnalysis(),
  def_hists();
  set_selectors();
  initial_state();
  cout << "AnaTda::Init done" << endl;
  return kSUCCESS;
}

void AnaTda::def_hists() {
  def_eff_hists();
  def_resid_hists();
  def_gamma_from_pi0_hists();
  def_elecs_from_jpsi_hists();
  def_tutorial_hists();
  def_pair_hists();
  def_single_hists();
  def_full_sys_hists();
  def_kin_fit_hists();
  for (int it=4; it<nhist; ++it) def_manual_kin_fit_hists(it);
}

void AnaTda::set_selectors() {
  // *** Mass selector for the jpsi cands
  m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi
  m0_pi0 = TDatabasePDG::Instance()->GetParticle("pi0")->Mass();   // Get nominal PDG mass of the J/psi
  cout << "AnaTda::Init() m0_jpsi= " << m0_jpsi << " m0_pi0= " << m0_pi0 << endl;
  jpsiMassSel = new RhoMassParticleSelector("jpsi",m0_jpsi,0.137); // Not exactly corresponds to binsong thesis (2.96 - 3.23)
  pi0MassSel = new RhoMassParticleSelector("pi0",m0_pi0,0.06);
}

void AnaTda::initial_state(){
  // *** the lorentz vector of the initial jpsi-pi0 system
  const double mprot= 0.938;
  const double E_antip = TMath::Hypot(mprot, p_antip);
  const double beta_cm = p_antip/(E_antip + mprot);
  ini = TLorentzVector(0, 0, p_antip, mprot+E_antip);
  cout << "betac_cm = " << beta_cm << endl;
  boost_to_cm.SetZ(-beta_cm);
  boost_to_cm.Print();
  sqrt_s = TMath::Sqrt(2*mprot*(mprot+E_antip));
  cout << "E_cm = " << sqrt_s << endl;
}

void AnaTda::Exec(Option_t* opt) {

  if (verb) {
    cout << endl << endl;
    cout << "===== AnaTda::Exec -- Event " << nevt << " ====="<< endl;
    ++nevt;
  } else {
    if (++nevt%100 == 0)
      cout << "===== AnaTda::Exec -- Event " << nevt << " ====="<< endl;
  }

  // clean up for next pass
  cleanup_rho_cand_lists();

  fAnalysis->GetEvent();

  get_singles_lists();

  if (verb)
    print_mc_list();

  fill_single_dists();
  truth_match_singles();
  fill_single_dists_tr();

  truth_match_residuals();

  make_pair_lists();
  fill_pair_dists();

  find_primary_gg();
  pi0_analysis_cut(); // this preps pi0_ana

  find_primary_epem();
  jpsi_analysis_cut();

  // Only do this part if truth matching is not ambiguous
  if (pi0_true.GetLength()==1 && jpsi_true.GetLength()==1) {
    pi0jpsi_efficiency(pi0_true, jpsi_true, eff_ref);
    pi0jpsi_efficiency(pi0_pm_ana, jpsi_true, eff_pi0sel);
    pi0jpsi_efficiency(pi0_pm_ana, jpsi_pm_ana, eff_jpsisel);
    pi0jpsi_efficiency(pi0_pm_ana, jpsi_pm_ana, eff_kin);
    pi0jpsi_efficiency(pi0_pm_ana, jpsi_pm_ana, eff_excl);

    pi0jpsi_kinematics(pi0_ana,jpsi_ana,ana);
    pi0jpsi_kinematics(pi0_pm_ana,jpsi_pm_ana,pm_ana);
    pi0jpsi_true_kinematics(pi0_true,jpsi_true);
    pi0jpsi_kin_fit(pi0_pm_ana, jpsi_pm_ana);
    pi0jpsi_efficiency(pi0_pm_ana, jpsi_pm_ana, eff_const);
  }

  //pi0jpsi_ana.Combine(pi0_ana,jpsi_ana);
  //pi0jpsi_pm_ana.Combine(pi0_pm_ana,jpsi_pm_ana);
  //pi0_kinematic_selection(pi0,jpsi,rec);
  //pdgm_nearest_pi0s();
  //pi0_kinematic_selection();
  //kin_fit_all();
  //kin_fit_pi0_nearest_all();

}

void AnaTda::def_pair_hists() {
  h_m_pippim = new TH1F("h_m_pippim",
    "Invariant mass of #pi^{+}#pi^{-};M_{#pi^{+}#pi^{-}}[GeV/c^{2}]", 100, 0, 4.5);
  h_m_pippim_tr = new TH1F("h_m_pippim_tr",
    "Invariant mass of #pi^{+}#pi^{-} (Truth Match);M_{#pi^{+}#pi^{-}}[GeV/c^{2}]", 100, 0, 4.5);
}

void AnaTda::def_resid_hists() {
  double m_ph[4] = {0.1,0.3,1.0,2.1*TMath::Pi()};
  double m_th[4] = {0.04,0.1,0.4,2.1*TMath::Pi()};
  double m_mom[4] = {0.1,0.5,5,20.};
  for (int i=0; i<4; ++i) {
    // photons
    h_resid_phth[i] = new
      TH2F(Form("h_resid_phth%d", i),
	   "Residual of neutral cands. from true photons;d#phi[rad];d#theta[rad]",
	   200, -m_ph[i], m_ph[i], 200, -m_th[i], m_th[i]);

    // Electrons
    h_resid_ep_phth[i] = new
      TH2F(Form("h_resid_ep_phth%d", i),
	   "Residual of ep cands. from true electrons;d#phi[rad];d#theta[rad]",
	   200, -m_ph[i], m_ph[i], 200, -m_th[i], m_th[i]);
    h_resid_ep_mom[i] = new
      TH1F(Form("h_resid_ep_mom%d", i),
	   "Residual of ep cands. from true electrons;dp[GeV/c]",
	   200, -m_mom[i], m_mom[i]);
    h_resid_em_phth[i] = new
      TH2F(Form("h_resid_em_phth%d", i),
	   "Residual of em cands. from true electrons;d#phi[rad];d#theta[rad]",
	   200, -m_ph[i], m_ph[i], 200, -m_th[i], m_th[i]);
    h_resid_em_mom[i] = new
      TH1F(Form("h_resid_em_mom%d", i),
	   "Residual of em cands. from true electrons;dp[GeV/c]",
	   200, -m_mom[i], m_mom[i]);

    // pi+ (pion hypothesis)
    h_resid_pip_phth[i] = new
      TH2F(Form("h_resid_pip_phth%d", i),
	   "Residual of pip cands. from true pions (pi hyp.);d#phi[rad];d#theta[rad]",
	   200, -m_ph[i], m_ph[i], 200, -m_th[i], m_th[i]);
    h_resid_pip_mom[i] = new
      TH1F(Form("h_resid_pip_mom%d", i),
	   "Residual of pip cands. from true pions (pi hyp.);dp[GeV/c]",
	   200, -m_mom[i], m_mom[i]);
    h_resid_pim_phth[i] = new
      TH2F(Form("h_resid_pim_phth%d", i),
	   "Residual of pim cands. from true pions (pi hyp.);d#phi[rad];d#theta[rad]",
	   200, -m_ph[i], m_ph[i], 200, -m_th[i], m_th[i]);
    h_resid_pim_mom[i] = new
      TH1F(Form("h_resid_pim_mom%d", i),
	   "Residual of pim cands. from true pions (pi hyp.);dp[GeV/c]",
	   200, -m_mom[i], m_mom[i]);

    // pi- (electron hypothesis)
    h_resid_pip_elec_hyp_phth[i] = new
      TH2F(Form("h_resid_pip_elec_hyp_phth%d", i),
	   "Residual of pip cands. from true pions (elec. hyp.);d#phi[rad];d#theta[rad]",
	   200, -m_ph[i], m_ph[i], 200, -m_th[i], m_th[i]);
    h_resid_pip_elec_hyp_mom[i] = new
      TH1F(Form("h_resid_pip_elec_hyp_mom%d", i),
	   "Residual of pip cands. from true pions (elec. hyp.);dp[GeV/c]",
	   200, -m_mom[i], m_mom[i]);
    h_resid_pim_elec_hyp_phth[i] = new
      TH2F(Form("h_resid_pim_elec_hyp_phth%d", i),
	   "Residual of pim cands. from true pions (elec. hyp.);d#phi[rad];d#theta[rad]",
	   200, -m_ph[i], m_ph[i], 200, -m_th[i], m_th[i]);
    h_resid_pim_elec_hyp_mom[i] = new
      TH1F(Form("h_resid_pim_elec_hyp_mom%d", i),
	   "Residual of pim cands. from true pions (elec. hyp.);dp[GeV/c]",
	   200, -m_mom[i], m_mom[i]);

  }
}

void AnaTda::def_gamma_from_pi0_hists() {
  const char *n[nhist] = {"all_rec", "true_rec", "truepi0_rec", "truepi0_mc", "ana", "pm_ana", "mconst"};
  const char *t[nhist] = {"all (Reco)", "truth match (Reco)", "from true #pi^{0} (Reco)",
			  "from true #pi^{0} (MC)", "after analysis cuts",
			  "from true #pi^{0} after ana. cut", "from true #pi^{0} after mass const"};

  const double mgg_max = 0.25;
  for (int it = 0; it < nhist; ++it) {
    h_m_gg[it] =
      new TH1F(Form("h_m_gg_%s", n[it]),
	       Form("Mass of #gamma-#gamma pairs %s;M_{#gamma#gamma}[GeV/c^{2}]", t[it]),
	       100, 0, mgg_max);
    h_oa_gg[it] =
      new TH1F(Form("h_oa_gg_%s", n[it]),
	       Form("Opening angle of #gamma-#gamma pairs %s;OA[rad]", t[it]),
	       200, 0, TMath::Pi()/2.);
    double emax = 0.8;
    if (iplab==1) emax = 1.9;
    if (iplab==2) emax = 3.8;
    h_oa_gg_vs_min_e_g[it] =
      new TH2F(Form("h_oa_gg_min_e_g_%s", n[it]),
	       Form("Min(E_{#gamma 1}, E_{#gamma 2}) vs. OA of #gamma-#gamma pairs %s;OA[rad];Min(E_{#gamma 1}, E_{#gamma 2})", t[it]),
	       200, 0, TMath::Pi()/2., 200, 0, emax);
    h_oa_gg_vs_avg_e_g[it] =
      new TH2F(Form("h_oa_gg_avg_e_g_%s", n[it]),
	       Form("Avg(E_{#gamma 1}, E_{#gamma 2}) vs. OA of #gamma-#gamma pairs %s;OA[rad];Avg(E_{#gamma 1}, E_{#gamma 2})", t[it]),
	       200, 0, TMath::Pi(), 200, 0, emax);
    h_oa_gg_vs_asym_e_g[it] =
      new TH2F(Form("h_oa_gg_asym_e_g_%s", n[it]),
	       Form("Asym(E_{#gamma 1}, E_{#gamma 2}) vs. OA of #gamma-#gamma pairs %s;OA[rad];Asym(E_{#gamma 1}, E_{#gamma 2})", t[it]),
	       200, 0, TMath::Pi(), 200, 0, 1.1);
    for (int ia = 0; ia < npi0ana; ++ia) {
      if (it < 4) continue;
      h_m_gg_pi0ana[ia][it] = new TH1F(Form("h_m_gg_%s_%d", n[it], ia),
        Form("Mass of #gamma-#gamma pairs %s (Ana. Cut %d);M_{#gamma#gamma}[GeV/c^{2}]", t[it], ia),
        100, 0, mgg_max);
    }
    h_e_g[it] = new TH1F(Form("h_e_g_%s", n[it]),
      Form("Energy of #gamma %s;E[GeV]", t[it]),
      100, 0, 1.6);
  }
}

void AnaTda::def_elecs_from_jpsi_hists() {
  const char *n[nhist] = {"all_rec", "true_rec", "truejpsi_rec", "truejpsi_mc", "ana", "pm_ana", "mconst"};
  const char *t[nhist] = {"all (Reco)", "truth match (Reco)", (bg_mc?"from true #pi^{+}#pi^{-} (Reco)":"from true J/#psi (Reco)"),
			  (bg_mc?"from true #pi^{+}#pi^{-} (MC)":"from true J/#psi (MC)"), "after analysis cuts",
			  (bg_mc?"from true #pi^{+}#pi^{-} after ana. cut":"from true J/#psi after ana. cut"),
			  (bg_mc?"from true #pi^{+}#pi^{-} after mass constraint":"from true J/#psi after mass constraint"),
  };

  for (int it = 0; it < nhist; ++it) {
    h_m_epem[it] = new TH1F(Form("h_m_epem_%s", n[it]),
      Form("Mass of e^{+}e^{-} pairs %s;M_{e^{+}e^{-}}[GeV/c^{2}]", t[it]),
      100, 0, 5.0);
    h_mom_epm[it] = new TH1F(Form("h_mom_epm_%s", n[it]),
      Form("Momentum of e^{+} and e^{-} %s;p[GeV/c]", t[it]),
      100, 0, 6.0);
    h_mom_the_epm[it] = new TH2F(Form("h_mom_the_epm_%s", n[it]),
      Form("#theta vs. momentum of e^{+} and e^{-} %s;p[GeV/c];#theta[rad]", t[it]),
      100, 0, 6.0, 100, 0, TMath::Pi());
  }
}

void AnaTda::def_single_hists() {
  h_num_g = new TH1F("h_num_g",
    "Number of #gamma per event;N_{#gamma}", 25, 0, 25);
  h_num_epm = new TH1F("h_num_epm",
    "Number of electron tracks per event;N_{e^{#pm}}", 10, 0, 10);
  h_num_pipm = new TH1F("h_num_pipm",
    "Number of charged pion tracks per event;N_{#pi^{#pm}}", 10, 0, 10);
  h_mom_the_pipm = new TH2F("h_mom_the_pipm",
    "#theta vs. momentum of charge pion tracks;p[GeV/c];#theta[rad]", 100, 0, 5, 100, 0, TMath::Pi());

  h_num_g_tr = new TH1F("h_num_g_tr",
    "Number of #gamma per event (Truth Match);N_{#gamma}", 25, 0, 25);
  h_num_epm_tr = new TH1F("h_num_epm_tr",
    "Number of electron tracks per event (Truth Match);N_{e^{#pm}}", 10, 0, 10);
  h_num_pipm_tr = new TH1F("h_num_pipm_tr",
    "Number of charged pion tracks per event (Truth Match);N_{#pi^{#pm}}", 10, 0, 10);
  h_mom_the_pipm_tr = new TH2F("h_mom_the_pipm_tr",
    "#theta vs. momentum of charge pion tracks (Truth Match);p[GeV/c];#theta[rad]", 100, 0, 5, 100, 0, TMath::Pi());
}

void AnaTda::cleanup_rho_cand_lists() {

  mcList.Cleanup();

  ep.Cleanup();
  em.Cleanup();
  pip.Cleanup();
  pim.Cleanup();
  g.Cleanup();

  pip_tr.Cleanup();
  pim_tr.Cleanup();
  ep_tr.Cleanup();
  em_tr.Cleanup();
  g_tr.Cleanup();

  epem.Cleanup();
  pippim.Cleanup();
  gg.Cleanup();
  epem_tr.Cleanup();
  pippim_tr.Cleanup();
  gg_tr.Cleanup();

  pi0.Cleanup();
  pi0_true.Cleanup();
  pi0_ana.Cleanup();
  pi0_pm_ana.Cleanup();
  for (int i = 0; i < npi0ana; ++i) {
    pi0_ana_[i].Cleanup();
    pi0_pm_ana_[i].Cleanup();
  }

  jpsi.Cleanup();
  jpsi_true.Cleanup();
  jpsi_ana.Cleanup();
  jpsi_pm_ana.Cleanup();

  jpsi_mconst.Cleanup();


  pi0jpsi_ana.Cleanup();
  pi0jpsi_pm_ana.Cleanup();

  pi0nearest.Cleanup();
  pi0_btb.Cleanup();
  pi0_cts.Cleanup();

  epem_mcut.Cleanup();
  pippim_mcut.Cleanup();

  epem_mcut_pi0_btb.Cleanup();
  epem_mcut_pi0_cts.Cleanup();
  pippim_mcut_pi0_btb.Cleanup();
  pippim_mcut_pi0_cts.Cleanup();
  epem_pi0nearest.Cleanup();
  pippim_pi0nearest.Cleanup();

}

// utilities
void AnaTda::print_rho_cand_list(RhoCandList& l, const char* name) {
  cout << "========================== " << name << " ========================== " << endl;
  for (int i = 0; i < l.GetLength(); ++i) {
    cout << "i= " << i << " ";
    l[i]->PrintOn(cout);
    cout << endl;
  }
}

// Print MC Truth list with mother-daughter relations
void AnaTda::print_mc_list() {
  cout << "mcList.Length() = " << mcList.GetLength() << endl;
  int pi0id = -1;
  for (int j=0;j<mcList.GetLength();++j) {
    RhoCandidate *mcmother = mcList[j]->TheMother();
    int muid = -1;
    if (mcmother) muid = mcmother->GetTrackNumber();
    if (muid==-1 and mcList[j]->PdgCode()==111) pi0id=j;
    if ( not (muid==-1 or (pi0id!=-1&&muid==pi0id)) ) continue;
    cout << "Track "<< mcList[j]->GetTrackNumber()<<" (PDG:"<<mcList[j]->PdgCode() <<") has mother "<<muid;
    if (mcList[j]->NDaughters()>0) cout <<" and daughter(s) ";
    for (int k=0;k<mcList[j]->NDaughters();++k) cout <<mcList[j]->Daughter(k)->GetTrackNumber()<<"  ";
    cout<<endl;
  }
  cout <<endl;
}

void AnaTda::get_singles_lists() {
  // *** Select with no PID info ('All'); type and mass are set
  fAnalysis->FillList(mcList, "McTruth");
  fAnalysis->FillList(ep, (fBremCorr?"BremElectronAllPlus":"ElectronAllPlus"));
  fAnalysis->FillList(em, (fBremCorr?"BremElectronAllMinus":"ElectronAllMinus"));
  fAnalysis->FillList(g, "Neutral");
  fAnalysis->FillList(pip, (fBremCorr?"BremPionAllPlus":"PionAllPlus"));
  fAnalysis->FillList(pim, (fBremCorr?"BremPionAllMinus":"PionAllMinus"));
}

void AnaTda::fill_single_e(RhoCandList& org, TH1* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->Energy());
}

void AnaTda::fill_single_e(RhoCandList& l1, RhoCandList& l2, TH1* dest) {
  fill_single_e(l1,dest);
  fill_single_e(l2,dest);
}

void AnaTda::fill_single_mom(RhoCandList& org, TH1* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->P3().Mag());
}

void AnaTda::fill_single_mom(RhoCandList& l1, RhoCandList& l2, TH1* dest) {
  fill_single_mom(l1,dest);
  fill_single_mom(l2,dest);
}

void AnaTda::fill_single_mom_the(RhoCandList& org, TH2F* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->P3().Mag(), org[j]->P3().Theta());
}

void AnaTda::fill_single_mom_the(RhoCandList& l1, RhoCandList& l2, TH2F* dest) {
  fill_single_mom_the(l1,dest);
  fill_single_mom_the(l2,dest);
}

void AnaTda::fill_single_dists() {
  // Single distributions
  h_num_g->Fill(g.GetLength());
  h_num_epm->Fill(ep.GetLength()+em.GetLength());
  h_num_pipm->Fill(pip.GetLength()+pim.GetLength());

  fill_single_e(g, h_e_g[all]);
  fill_single_mom(ep, em, h_mom_epm[all]);
  fill_single_mom_the(ep, em, h_mom_the_epm[all]);

  for (int i = 0; i < pip.GetLength(); ++i)
    h_mom_the_pipm->Fill(pip[i]->P3().Mag(), pip[i]->P3().Theta() );
  for (int i = 0; i < pim.GetLength(); ++i)
    h_mom_the_pipm->Fill(pim[i]->P3().Mag(), pim[i]->P3().Theta() );
}

void AnaTda::truth_match(RhoCandList& org, RhoCandList& dest, const int &pdg) {
  org.SetType(pdg);
  for (int i = 0; i < org.GetLength(); ++i)
    if (fAnalysis->McTruthMatch(org[i])) dest.Append(org[i]);
}

void AnaTda::truth_match_singles() {
  truth_match(ep, ep_tr, -11);
  truth_match(em, em_tr, 11);
  truth_match(pip, pip_tr, pdg_pip);
  truth_match(pim, pim_tr, pdg_pim);
  truth_match(g, g_tr, 22);
}

void AnaTda::fill_single_dists_tr() {
  // Single distributions (After truth match)
  h_num_g_tr->Fill(g_tr.GetLength());
  h_num_epm_tr->Fill(ep_tr.GetLength()+em_tr.GetLength());
  h_num_pipm_tr->Fill(pip_tr.GetLength()+pim_tr.GetLength());

  fill_single_e(g_tr, h_e_g[tr]);
  fill_single_mom(ep_tr, em_tr, h_mom_epm[tr]);
  fill_single_mom_the(ep_tr, em_tr, h_mom_the_epm[tr]);

  for (int i = 0; i < pip_tr.GetLength(); ++i)
    h_mom_the_pipm_tr->Fill(pip_tr[i]->P3().Mag(), pip_tr[i]->P3().Theta() );
  for (int i = 0; i < pim_tr.GetLength(); ++i)
    h_mom_the_pipm_tr->Fill(pim_tr[i]->P3().Mag(), pim_tr[i]->P3().Theta() );
}

void AnaTda::truth_match_residuals() {
  // Find the two "primary" photons that decayed from the pi0s
  static int ndalitz = 0;
  for (int j=0;j<mcList.GetLength();++j) {
    const int pdg = mcList[j]->PdgCode();
    RhoCandidate *mc_p = mcList[j]->TheMother();
    const int p_tnum = mc_p?mc_p->GetTrackNumber():-1;
    const int p_pdg = mc_p?mc_p->PdgCode():-1;
    if (pdg == pdg_pi0 and ((bg_mc and p_tnum == -1) or (!bg_mc and p_tnum == 0)) ) {

      // we got primary pion
      if (mcList[j]->NDaughters()!=2) {
	ndalitz++;
      	cout << "Nevt= " << nevt << " Pion with more than two daughters (prolly dalitz): " << mcList[j]->NDaughters() << " ndalitz= " << ndalitz << endl;
      	continue;
      }
      // Histogram the phi and theta residuals of all reconstructed neutrals
      for (int idaughter=0; idaughter<2; ++idaughter) {
	RhoCandidate *_g = mcList[j]->Daughter(idaughter);
	for (int igrec = 0; igrec<g.GetLength(); ++igrec) {
	  const double dph = _g->P3().Phi() - g[igrec]->P3().Phi();
	  const double dth = _g->P3().Theta() - g[igrec]->P3().Theta();
	  for (int i=0; i<4; ++i) h_resid_phth[i]->Fill(dph, dth);
	}
      }

    } else if (pdg==pdg_pip and p_tnum == -1) {

      // we got primary charged pion+
      RhoCandidate *_pip = mcList[j];
      for (int ipiprec = 0; ipiprec<pip.GetLength(); ++ipiprec) {
	const double dmom = _pip->P3().Mag() - pip[ipiprec]->P3().Mag();
	const double dph = _pip->P3().Phi() - pip[ipiprec]->P3().Phi();
	const double dth = _pip->P3().Theta() - pip[ipiprec]->P3().Theta();
	for (int i=0; i<4; ++i) {
	  h_resid_pip_mom[i]->Fill(dmom);
	  h_resid_pip_phth[i]->Fill(dph, dth);
	}
      }
      for (int ieprec = 0; ieprec<ep.GetLength(); ++ieprec) {
	const double dmom = _pip->P3().Mag() - ep[ieprec]->P3().Mag();
	const double dph = _pip->P3().Phi() - ep[ieprec]->P3().Phi();
	const double dth = _pip->P3().Theta() - ep[ieprec]->P3().Theta();
	for (int i=0; i<4; ++i) {
	  h_resid_pip_elec_hyp_mom[i]->Fill(dmom);
	  h_resid_pip_elec_hyp_phth[i]->Fill(dph, dth);
	}
      }
    }  else if (pdg==pdg_pim and p_tnum == -1) {
      // we got primary charged pion-
      RhoCandidate *_pim = mcList[j];
      for (int ipimrec = 0; ipimrec<pim.GetLength(); ++ipimrec) {
	const double dmom = _pim->P3().Mag() - pim[ipimrec]->P3().Mag();
	const double dph = _pim->P3().Phi() - pim[ipimrec]->P3().Phi();
	const double dth = _pim->P3().Theta() - pim[ipimrec]->P3().Theta();
	for (int i=0; i<4; ++i) {
	  h_resid_pim_mom[i]->Fill(dmom);
	  h_resid_pim_phth[i]->Fill(dph, dth);
	}
      }
      for (int iemrec = 0; iemrec<em.GetLength(); ++iemrec) {
	const double dmom = _pim->P3().Mag() - em[iemrec]->P3().Mag();
	const double dph = _pim->P3().Phi() - em[iemrec]->P3().Phi();
	const double dth = _pim->P3().Theta() - em[iemrec]->P3().Theta();
	for (int i=0; i<4; ++i) {
	  h_resid_pim_mom[i]->Fill(dmom);
	  h_resid_pim_phth[i]->Fill(dph, dth);
	}
      }
    } else if (pdg==pdg_ep and p_pdg == pdg_jpsi) {

      // we got primary e+
      RhoCandidate *_ep = mcList[j];
      for (int ieprec = 0; ieprec<ep.GetLength(); ++ieprec) {
	const double dmom = _ep->P3().Mag() - ep[ieprec]->P3().Mag();
	const double dph = _ep->P3().Phi() - ep[ieprec]->P3().Phi();
	const double dth = _ep->P3().Theta() - ep[ieprec]->P3().Theta();
	for (int i=0; i<4; ++i) {
	  h_resid_ep_mom[i]->Fill(dmom);
	  h_resid_ep_phth[i]->Fill(dph, dth);
	}
      }
    } else if (pdg==pdg_em and p_pdg == pdg_jpsi) {
      // we got primary e-
      RhoCandidate *_em = mcList[j];
      for (int iemrec = 0; iemrec<em.GetLength(); ++iemrec) {
	const double dmom = _em->P3().Mag() - em[iemrec]->P3().Mag();
	const double dph = _em->P3().Phi() - em[iemrec]->P3().Phi();
	const double dth = _em->P3().Theta() - em[iemrec]->P3().Theta();
	for (int i=0; i<4; ++i) {
	  h_resid_em_mom[i]->Fill(dmom);
	  h_resid_em_phth[i]->Fill(dph, dth);
	}
      }
    }

  }
}

void AnaTda::make_pair_lists() {
  // candidate lists for pairwise truth match (full hierarchy)
  epem.Combine(ep, em);
  epem.SetType(pdg_jpsi);
  epem_tr.Combine(ep_tr, em_tr);
  gg.Combine(g, g);
  gg_tr.Combine(g_tr, g_tr);
}

void AnaTda::fill_pair_mass(RhoCandList& org, TH1F* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->M());
}

void AnaTda::fill_pair_oa(RhoCandList& org, TH1F* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(oa(org[j]->Daughter(0),org[j]->Daughter(1)));
}

void AnaTda::fill_pair_oa(RhoCandidate* _g1, RhoCandidate* _g2, const int &itype) {
  const double _oa = oa(_g1,_g2);
  h_oa_gg[itype]->Fill(_oa);
  h_oa_gg_vs_min_e_g[itype]->Fill(_oa, fmin(_g1->Energy(), _g2->Energy()) );
  h_oa_gg_vs_avg_e_g[itype]->Fill(_oa, 0.5*(_g1->Energy() + _g2->Energy()) );
  h_oa_gg_vs_asym_e_g[itype]->Fill(_oa, fabs(_g1->Energy() - _g2->Energy())/(_g1->Energy() + _g2->Energy()) );
}

void AnaTda::fill_pair_oa(RhoCandList& org, const int &itype) {
  for (int j = 0; j < org.GetLength(); ++j) fill_pair_oa(org[j]->Daughter(0), org[j]->Daughter(1), itype);
}

void AnaTda::fill_pair_oa_mc(RhoCandList& org, const int &itype) {
  for (int j = 0; j < org.GetLength(); ++j) fill_pair_oa(org[j]->Daughter(0)->GetMcTruth(), org[j]->Daughter(1)->GetMcTruth(), itype);
}

void AnaTda::fill_pair_dists() {
  fill_pair_mass(epem, h_m_epem[all]);
  fill_pair_mass(epem_tr, h_m_epem[tr]);
  //fill_pair_mass(pippim, h_m_pippim);
  //fill_pair_mass(pippim_tr, h_m_pippim_tr);
  fill_pair_mass(gg, h_m_gg[all]);
  fill_pair_mass(gg_tr, h_m_gg[tr]);
  //fill_pair_oa(gg, h_oa_gg[all]);
  //fill_pair_oa(gg_tr, h_oa_gg[tr]);
  fill_pair_oa(gg, all);
  fill_pair_oa(gg_tr, tr);
}

inline
double AnaTda::oa(RhoCandidate* c1, RhoCandidate* c2) {
  return c1->P3().Angle(c2->P3());
}

inline
double AnaTda::mass(RhoCandidate* c1, RhoCandidate* c2) {
  return (c1->P4() + c2->P4()).M();
}

void AnaTda::primary_match(const int& pdgm, RhoCandidate* r, RhoCandidate* m) {
  for (int j=0;j<mcList.GetLength();++j) {
    RhoCandidate *mcmother = mcList[j]->TheMother();
    if (!mcmother) {
      *r = *mcList[j];
    } else if ( r and mcmother->GetTrackNumber() == r->GetTrackNumber() ) {
      if ( mcList[j]->PdgCode() == pdgm ) {
        *m = *mcList[j];
        break;
      }
    }
  }
}

bool AnaTda::primary_match_single(const int& pdgm, RhoCandidate* c, double &dist) {
  c->SetType(pdgm); // just to make sure in case it was forgotten .. just careful!
  if (!fAnalysis->McTruthMatch(c)) return false;
  RhoCandidate _root, _m;
  primary_match(pdgm,&_root,&_m);
  if (_root.GetTrackNumber() != 0 or _m.GetTrackNumber() == -1 ) return false;
  dist = c->Energy()-_m.Energy();
  return c->GetMcTruth()->GetTrackNumber() == _m.GetTrackNumber();
}

inline
double AnaTda::dist_photon_match(RhoCandidate *rec, RhoCandidate *mc) {
  const double nsig = bg_mc? 10.0 : 3.0; // TODO - this is hacky
  const double dth = (rec->P3().Theta()-mc->P3().Theta())/8.0e-3;
  const double dph = (rec->P3().Phi()-mc->P3().Phi())/8.22e-3;
  //cout << "phot dist= " << TMath::Hypot(dph, dth) << endl;
  return TMath::Hypot(dph, dth) / nsig; // TODO - why this has to be blown up for photon finding efficiency to be good?
}

inline
double AnaTda::dist_chpi_match(RhoCandidate *rec, RhoCandidate *mc) {
  const double dmom = (rec->P3().Mag()-mc->P3().Mag())/1.3e-2;
  const double dth = (rec->P3().Theta()-mc->P3().Theta())/1.54e-3;
  const double dph = (rec->P3().Phi()-mc->P3().Phi())/3.95e-3;
  return TMath::Hypot(dmom, TMath::Hypot(dth, dph)) / 3.0;
}

double AnaTda::primary_match_pair(const int& pdgm, RhoCandidate* c, int &tnp, int &tn0, int &tn1) {
  RhoCandidate *p_m, *p_d0, *p_d1;

  if (bg_mc && pdgm==pdg_jpsi) {
    int nchpi = 0;
    for (int j=0;j<mcList.GetLength();++j) {
      if (mcList[j]->PdgCode()!=pdg_pip and mcList[j]->PdgCode()!=pdg_pim) continue;
      RhoCandidate *mcmother = mcList[j]->TheMother();
      if (!mcmother) {
	if (nchpi==0) p_d0 = mcList[j];
	else p_d1 = mcList[j];
	nchpi++;
      }
      if (nchpi>=2) break;
    }
    if (!p_d0 or !p_d1) return -9999.;
  } else {

    for (int j=0;j<mcList.GetLength();++j) {
      if (mcList[j]->PdgCode()!=pdgm) continue;
      RhoCandidate *mcmother = mcList[j]->TheMother();
      if ((bg_mc and mcmother) or (!bg_mc and mcmother->GetTrackNumber()!=0 ) ) { // For some reason, track0 is initial state in Signal simulation
	cout << "WARNING: First particle with right pdg code " << pdgm << " should be primary but is not" << endl;
	break;
      }
      if (mcList[j]->NDaughters()!=2) {
	//cout << "WARNING: Primary has " << mcList[j]->NDaughters() << " daughters, expecting 2" << endl;
	break;
      }
      p_m = mcList[j];
      p_d0 = mcList[j]->Daughter(0);
      p_d1 = mcList[j]->Daughter(1);
      break;
    }
    if (!p_d0 or !p_d1){
      return -9999.;}
  }

  if (pdgm==pdg_jpsi) {
    if (!bg_mc) {
      // For electrons actually from J/psis use Framwork match to calculate distance
      RhoCandidate *d0 = c->Daughter(0)->GetMcTruth();
      RhoCandidate *d1 = c->Daughter(1)->GetMcTruth();
      if (!d0 or !d1 ) {
	return -9999.;
      } else {
	const int r_tn0 = d0->GetTrackNumber();
	const int r_tn1 = d1->GetTrackNumber();
	const int p_tn0 = p_d0->GetTrackNumber();
	const int p_tn1 = p_d1->GetTrackNumber();
	if (r_tn0==p_tn0 and r_tn1==p_tn1) {
	  tn0 = p_tn0; tn1 = p_tn1; tnp = p_m->GetTrackNumber();
	  return TMath::Hypot(d0->Energy()-p_d0->Energy(), d1->Energy()-p_d1->Energy());
	} else if (r_tn0==p_tn1 and r_tn1==p_tn0) {
	  tn0 = p_tn1; tn1 = p_tn0; tnp = p_m->GetTrackNumber();
	  return TMath::Hypot(d0->Energy()-p_d1->Energy(), d1->Energy()-p_d0->Energy());
	} else {
	  return -9999.;
	}
      }
    } else {

      // For BG MC, some other criterion should be used because nothing comes from J/psi
      const double d00 = dist_chpi_match(c->Daughter(0),p_d0);
      const double d11 = dist_chpi_match(c->Daughter(1),p_d1);
      const double d01 = dist_chpi_match(c->Daughter(0),p_d1);
      const double d10 = dist_chpi_match(c->Daughter(1),p_d0);
      const double m0 = TMath::Hypot(d00,d11);
      const double m1 = TMath::Hypot(d01,d10);
      if (m0 < m1) {
	tn0 = p_d0->GetTrackNumber();
	tn1 = p_d1->GetTrackNumber();
	return m0;
      } else {
	tn0 = p_d1->GetTrackNumber();
	tn1 = p_d0->GetTrackNumber();
	return m1;
      }
    }

  } else if (pdgm==pdg_pi0) {
    // For photons, verify that phi and theta are within tolerance
    double d00 = dist_photon_match(c->Daughter(0),p_d0);
    double d11 = dist_photon_match(c->Daughter(1),p_d1);
    double d01 = dist_photon_match(c->Daughter(0),p_d1);
    double d10 = dist_photon_match(c->Daughter(1),p_d0);
    if (d00<1.0 && d11<1.0) {
      tnp = p_m->GetTrackNumber();
      tn0 = p_d0->GetTrackNumber();
      tn1 = p_d1->GetTrackNumber();
      return TMath::Hypot(d00,d11);
    } else if (d01<1.0 && d10<1.0) {
      tnp = p_m->GetTrackNumber();
      tn0 = p_d1->GetTrackNumber();
      tn1 = p_d0->GetTrackNumber();
      return TMath::Hypot(d01,d10);
    } else {
      return -9999.;
    }
  } else {
    return -9999.;
  }

  cout << "This shouldn't happen. All cases should have been treated by now. if not, go debug" << endl;
}

inline
double AnaTda::dist_pi0_pair_match(RhoCandidate *c) {
  return TMath::Hypot(dist_photon_match(c->Daughter(0),c->Daughter(0)->GetMcTruth()),
		      dist_photon_match(c->Daughter(1),c->Daughter(1)->GetMcTruth()));
}

inline
double AnaTda::dist_jpsi_pair_match(RhoCandidate *c) {
  return TMath::Hypot(c->Daughter(0)->Energy()-c->Daughter(0)->GetMcTruth()->Energy(),
		      c->Daughter(1)->Energy()-c->Daughter(1)->GetMcTruth()->Energy());
}

//void AnaTda::find_primary_gg() {
//
//  if (verb)
//    cout << "=============== FIND PRIMARY NEUTRAL PAIR (ng= "  << g.GetLength() << ") ================" << endl;
//
//  RhoCandidate *_match;
//  double dist_min = 1e9;
//  pi0.SetType(pdg_pi0);
//  pi0.Combine(g, g);
//  int tp = -1;
//  int td[2] = {-1, -1};
//  bool set_mc_match = false;
//  for (int j = 0; j < pi0.GetLength(); ++j) {
//    pi0[j]->SetType(pdg_pi0); // just to make sure in case it was forgotten
//    // Let the framework do the matching first. This sets the matching to the full
//    // arborescence if it scucees to find a primary match. Otherwise do it manually for the
//    // case that fails (pi0->gg in BG simulation)
//    fAnalysis->McTruthMatch(pi0[j]);
//    bool framework_found = false;
//    double dist_framework = 0;
//    if (pi0[j]->GetMcTruth() and pi0[j]->Daughter(0)->GetMcTruth() and pi0[j]->Daughter(1)->GetMcTruth()){
//      framework_found =true;
//      dist_framework = dist_pi0_pair_match(pi0[j]);
//      if (dist_framework < dist_min) {
//	dist_min = dist_framework;
//	_match = pi0[j];
//	tp = pi0[j]->GetMcTruth()->GetTrackNumber();
//	td[0] = pi0[j]->Daughter(0)->GetMcTruth()->GetTrackNumber();
//	td[1] = pi0[j]->Daughter(1)->GetMcTruth()->GetTrackNumber();
//      }
//    }
//
//    int ip = 0, i0 = 0, i1 = 0;
//    double dist_custom = primary_match_pair(pdg_pi0, pi0[j], ip, i0, i1);
//    bool custom_found = false;
//    if ( dist_custom >= 0 ) {
//      custom_found = true;
//      if (dist_custom < dist_min) {
//        _match = pi0[j];
//        dist_min = dist_custom;
//	set_mc_match = true;
//	tp = ip;
//	td[0] = i0;
//	td[1] = i1;
//      }
//    }
//
//    if (framework_found or custom_found ) {
//      if (verb) {
//	cout << "------------------------------------------------------------" << endl;
//	if (framework_found) {
//	  cout << "Framework gg=(" << pi0[j]->Daughter(0)->Uid() << "," << pi0[j]->Daughter(1)->Uid()
//	       << "): " << pdg_pi0 << " M= " << pi0[j]->GetMcTruth()->GetTrackNumber()
//	       << " d0= " << pi0[j]->Daughter(0)->GetMcTruth()->GetTrackNumber()
//	       << " d1= " << pi0[j]->Daughter(1)->GetMcTruth()->GetTrackNumber()
//	       << " dist= " << dist_framework << endl;
//	} else {
//	  cout << "Framework    " << pdg_pi0 << ": no match" << endl;
//	}
//	if (custom_found)
//	  cout << "Custom:  gg=(" << pi0[j]->Daughter(0)->Uid() << "," << pi0[j]->Daughter(1)->Uid() << "): " << pdg_pi0 <<  " M= "
//	       << ip << " d0= " << i0 << " d1= " << i1 << " dist= " << dist_custom << endl;
//	else
//	  cout << "Custom    " << pdg_pi0 << ": no match" << endl;
//	cout << "------------------------------------------------------------" << endl;
//      }
//    }
//  }
//
//  if (dist_min<1e8) {
//    if (verb) {
//      cout << "===========================================================" << endl;
//      cout << "Custom Best Pair: TrackNum(pi0 " << _match->Daughter(0)->Uid() << "," << _match->Daughter(1)->Uid() << ")  M= "
//	   << _match->GetTrackNumber() << " d0= " << _match->Daughter(0)->GetTrackNumber()
//	   << " d1= " << _match->Daughter(1)->GetTrackNumber() << endl;
//      cout << "Custom Best Pair: TruthMatch(pi0 " << _match->Daughter(0)->Uid() << "," << _match->Daughter(1)->Uid() << ")  M= "
//	 << tp << " d0= " << td[0] << " d1= " << td[1] << " dist= " << dist_min << endl;
//    }
//    if (set_mc_match) {
//      if (verb) cout << "Framwork didn't find match. Manually setting match" << endl;
//      _match->SetMcTruth(mcList[tp]);
//      for (int i=0; i<2; ++i) _match->Daughter(i)->SetMcTruth(mcList[td[i]]);
//    }
//    if (verb) {
//      cout << "Custom  Check: TruthMatch(pi0): M= " << _match->GetMcTruth()->GetTrackNumber()
//	   << " d0= " << _match->Daughter(0)->GetMcTruth()->GetTrackNumber()
//	   << " d1= " << _match->Daughter(1)->GetMcTruth()->GetTrackNumber()  << endl;
//      cout << "===========================================================" << endl;
//    }
//    pi0_true.Append(_match);
//    fill_gamma_from_pi0s();
//  }  else {
//    if (g.GetLength()>1) {
//      cout << "Nevt= " << nevt << " valid reonstructed photon pair, but no priamry match found" << endl;
//    }
//  }
//}

// Much simpler version, but there may be more than one "true" pi0 found per event
void AnaTda::find_primary_gg() {
  if (verb)
    cout << "=============== FIND PRIMARY NEUTRAL PAIR (ng= "  << g.GetLength() << ") ================" << endl;
  pi0.SetType(pdg_pi0);
  pi0.Combine(g, g);
  for (int j = 0; j < pi0.GetLength(); ++j) {
    pi0[j]->SetType(pdg_pi0); // just to make sure in case it was forgotten
    fAnalysis->McTruthMatch(pi0[j]);
    if (pi0[j]->GetMcTruth() and pi0[j]->Daughter(0)->GetMcTruth() and pi0[j]->Daughter(1)->GetMcTruth()){
      pi0_true.Append(pi0[j]);
    }
  }
  fill_gamma_from_pi0s();
}

void AnaTda::fill_gamma_from_pi0s() {
  if (pi0_true.GetLength()!=1) return;
  RhoCandidate* _g1 = pi0_true[0]->Daughter(0);
  RhoCandidate* _g2 = pi0_true[0]->Daughter(1);
  RhoCandidate* _g1_mc = pi0_true[0]->Daughter(0)->GetMcTruth();
  RhoCandidate* _g2_mc = pi0_true[0]->Daughter(1)->GetMcTruth();

  h_e_g[pm]->Fill(_g1->Energy());
  h_e_g[pm]->Fill(_g2->Energy());
  h_e_g[pm_mc]->Fill(_g1_mc->Energy());
  h_e_g[pm_mc]->Fill(_g2_mc->Energy());

  h_m_gg[pm_mc]->Fill(mass(_g1_mc,_g2_mc));
  fill_pair_oa_mc(pi0_true,pm_mc);

  h_m_gg[pm]->Fill(mass(_g1,_g2));
  fill_pair_oa(pi0_true,pm);
}

void AnaTda::pi0_analysis_cut() {
  pi0_analysis_cut(pi0, pi0_ana, true);
  fill_pi0_analysis_hists(pi0_ana, ana);
  if (pi0_true.GetLength()==1){
    pi0_analysis_cut(pi0_true, pi0_pm_ana, true);
    fill_pi0_analysis_hists(pi0_pm_ana, pm_ana);
  }
  for (int i=0; i < npi0ana; ++i) {
    pi0_analysis_cut(pi0, pi0_ana_[i], i==1);
    fill_pi0_analysis_hists(pi0_ana_[i], i, ana);
    if (pi0_true.GetLength()==1) {
      pi0_analysis_cut(pi0_true, pi0_pm_ana_[i], i==1);
      fill_pi0_analysis_hists(pi0_pm_ana_[i], i, pm_ana);
    }
  }
}

bool AnaTda::oa_vs_avg_cut(const double& _oa, const double &_avg){
  bool acc = (_avg > (lw[iplab][2] + (lw[iplab][0]/(_oa-lw[iplab][1]))));
  if (_oa>up[iplab][1])
    acc = acc && (_avg < (up[iplab][2] + (up[iplab][0] / (_oa-up[iplab][1]))));
  return acc;
}

void AnaTda::pi0_analysis_cut(RhoCandList& org, RhoCandList& dest, bool mcut = false) {
  for (int i = 0; i < org.GetLength(); ++i) {
    RhoCandidate *_g1 = org[i]->Daughter(0);
    RhoCandidate *_g2 = org[i]->Daughter(1);
    const double _oa = oa(_g1,_g2);
    const double _e_avg = 0.5 * (_g1->Energy() + _g2->Energy());
    if ( oa_vs_avg_cut(_oa, _e_avg) ) {
      if (!mcut || (mcut && org[i]->M() > pi0mcut_min and org[i]->M() < pi0mcut_max) ) {
        dest.Append(org[i]);
      }
    }
  }
}

void AnaTda::fill_pi0_analysis_hists(RhoCandList& list, const int& itype) {
  fill_pair_mass(list, h_m_gg[itype]);
  fill_pair_oa(list, itype);
}

void AnaTda::fill_pi0_analysis_hists(RhoCandList& list, const int& id, const int& itype) {
  fill_pair_mass(list, h_m_gg_pi0ana[id][itype]);
}

//void AnaTda::find_primary_epem() {
//
//  RhoCandidate *_match;
//  double dist_min = 1e9;
//  jpsi.SetType(pdg_jpsi);
//  jpsi.Combine(ep, em);
//
//  if (verb)
//    cout << "=============== FIND PRIMARY CHARGED PAIR (nep= " << ep.GetLength() << " nem=" << em.GetLength()
//	 << " njpsi= " << jpsi.GetLength() << " ) ================" << endl;
//
//  bool set_mc_match = false;
//  int tp = -1;
//  int td[2] = {-1, -1};
//  for (int j = 0; j < jpsi.GetLength(); ++j) {
//    jpsi[j]->SetType(pdg_jpsi); // just to make sure in case it was forgotten
//    // Let the framework do the matching first. This sets the matching to the full
//    // arborescence if it scucees to find a primary match. Otherwise do it manually for the
//    // case that fails. For J/psi, framework usually succeeds
//    fAnalysis->McTruthMatch(jpsi[j]);
//    bool framework_found = false;
//    double dist_framework = 0;
//    if (jpsi[j]->GetMcTruth() and jpsi[j]->Daughter(0)->GetMcTruth() and jpsi[j]->Daughter(1)->GetMcTruth()){
//      framework_found =true;
//      dist_framework = dist_jpsi_pair_match(jpsi[j]);
//      if (dist_framework < dist_min) {
//	dist_min = dist_framework;
//	_match = jpsi[j];
//	tp = jpsi[j]->GetMcTruth()->GetTrackNumber();
//	td[0] = jpsi[j]->Daughter(0)->GetMcTruth()->GetTrackNumber();
//	td[1] = jpsi[j]->Daughter(1)->GetMcTruth()->GetTrackNumber();
//      }
//    }
//
//    int ip = 0, i0 = 0, i1 = 0;
//    double dist_custom = primary_match_pair(pdg_jpsi, jpsi[j], ip, i0, i1);
//    bool custom_found = false;
//    if ( dist_custom >= 0 ) {
//      custom_found =true;
//      if (dist_custom < dist_min) {
//        _match = jpsi[j];
//        dist_min = dist_custom;
//	set_mc_match = true;
//	tp = ip;
//	td[0] = i0;
//	td[1] = i1;
//      }
//    }
//
//    if (framework_found || custom_found ) {
//      if (verb) {
//
//	cout << "------------------------------------------------------------" << endl;
//	if (framework_found) {
//	  cout << "Framework ee=(" << jpsi[j]->Daughter(0)->Uid() << "," << jpsi[j]->Daughter(1)->Uid()
//	       << "): " << pdg_jpsi << " M= " << jpsi[j]->GetMcTruth()->GetTrackNumber()
//	       << " d0= " << jpsi[j]->Daughter(0)->GetMcTruth()->GetTrackNumber()
//	       << " d1= " << jpsi[j]->Daughter(1)->GetMcTruth()->GetTrackNumber()
//	       << " dist= " << dist_framework << endl;
//	} else {
//	  cout << "Framework    " << pdg_jpsi << ": no match" << endl;
//	}
//	if (custom_found)
//	  cout << "Custom:  ep=(" << jpsi[j]->Daughter(0)->Uid() << "," << jpsi[j]->Daughter(1)->Uid() << "): " << pdg_jpsi << " M= "
//	       << ip << " d0= " << i0 << " d1= " << i1 << " dist= " << dist_custom << endl;
//	else
//	  cout << "Custom    " << pdg_jpsi << ": no match" << endl;
//	cout << "------------------------------------------------------------" << endl;
//      }
//    }
//  }
//
//  if (dist_min<1e8) {
//    if (verb) {
//      cout << "===========================================================" << endl;
//      cout << "Custom Best Pair: TrackNum(psi " << _match->Daughter(0)->Uid() << "," << _match->Daughter(1)->Uid()
//	   << ")  M= " << _match->GetTrackNumber() << " d0= " << _match->Daughter(0)->GetTrackNumber()
//	   << " d1= " << _match->Daughter(1)->GetTrackNumber() << endl;
//      cout << "Custom Best Pair: TruthMatch(psi " << _match->Daughter(0)->Uid() << "," << _match->Daughter(1)->Uid() << ")  M= "
//	   << tp << " d0= " << td[0] << " d1= " << td[1] << " dist= " << dist_min << endl;
//    }
//    if (set_mc_match) {
//      if (verb) cout << "Framwork didn't find match. Manually setting match" << endl;
//      _match->SetMcTruth(mcList[tp]);
//      for (int i=0; i<2; ++i) _match->Daughter(i)->SetMcTruth(mcList[td[i]]);
//    }
//    if (verb){
//      cout << "Custom check: TruthMatch(psi): M= " << _match->GetMcTruth()->GetTrackNumber()
//	   << " d0= " << _match->Daughter(0)->GetMcTruth()->GetTrackNumber()
//	   << " d1= " << _match->Daughter(1)->GetMcTruth()->GetTrackNumber()  << endl;
//      cout << "===========================================================" << endl;
//    }
//    jpsi_true.Append(_match);
//    fill_elecs_from_jpsi();
//  } else {
//    if (jpsi.GetLength()>0) {
//      cout << "Nevt= " << nevt <<  " Valid charged pair reconstructed but coudln't be matched to primary charged pair" << endl;
//    }
//  }
//
//}

void AnaTda::find_primary_epem() {

  jpsi.Combine(ep, em);

  // bg and signal treated differently because by definition bg events don't have true electrons
  if (bg_mc != 0) {
    jpsi.SetType(pdg_jpsi);
    for (int j=0; j < jpsi.GetLength(); ++j) {
      jpsi[j]->SetType(pdg_jpsi);
      fAnalysis->McTruthMatch(jpsi[j]);
      if (jpsi[j]->GetMcTruth() and jpsi[j]->Daughter(0)->GetMcTruth() and jpsi[j]->Daughter(1)->GetMcTruth()){
	jpsi_true.Append(jpsi[j]);
      }
    }

  } else {
    RhoCandidate *_match;
    double dist_min = 1e9;
    bool custom_found = false;
    int tp = -1;
    int td[2] = {-1, -1};
    for (int j = 0; j < jpsi.GetLength(); ++j) {
      int ip = 0, i0 = 0, i1 = 0;
      double dist_custom = primary_match_pair(pdg_jpsi, jpsi[j], ip, i0, i1);
      if ( dist_custom >= 0 ) {
	custom_found =true;
	if (dist_custom < dist_min) {
	  _match = jpsi[j];
	  dist_min = dist_custom;
	  tp = ip;
	  td[0] = i0;
	  td[1] = i1;
	}
      }
    }
    if (custom_found) {
      _match->SetMcTruth(mcList[tp]);
      for (int i=0; i<2; ++i) _match->Daughter(i)->SetMcTruth(mcList[td[i]]);
      jpsi_true.Append(_match);
    }
  }
  fill_elecs_from_jpsi();
}

void AnaTda::find_primary_pippim() {

}

void AnaTda::fill_elecs_from_jpsi() {
  if (jpsi_true.GetLength()!=1) return;
  RhoCandidate* _e = jpsi_true[0]->Daughter(0);
  RhoCandidate* _p = jpsi_true[0]->Daughter(1);
  RhoCandidate* _e_mc = jpsi_true[0]->Daughter(0)->GetMcTruth();
  RhoCandidate* _p_mc = jpsi_true[0]->Daughter(1)->GetMcTruth();
  h_m_epem[pm]->Fill(mass(_e,_p));
  h_m_epem[pm_mc]->Fill(mass(_e_mc,_p_mc));
  h_mom_epm[pm]->Fill(_e->P3().Mag());
  h_mom_epm[pm]->Fill(_p->P3().Mag());
  h_mom_epm[pm_mc]->Fill(_e_mc->P3().Mag());
  h_mom_epm[pm_mc]->Fill(_p_mc->P3().Mag());
  h_mom_the_epm[pm]->Fill(_e->P3().Mag(),_e->P3().Theta());
  h_mom_the_epm[pm]->Fill(_p->P3().Mag(),_p->P3().Theta());
  h_mom_the_epm[pm_mc]->Fill(_e_mc->P3().Mag(),_e_mc->P3().Theta());
  h_mom_the_epm[pm_mc]->Fill(_p_mc->P3().Mag(),_p_mc->P3().Theta());
}

void AnaTda::jpsi_analysis_cut() {
  jpsi_analysis_cut(jpsi, jpsi_ana);
  fill_jpsi_analysis_hists(jpsi_ana,ana);
  if (jpsi_true.GetLength()==1) {
    jpsi_analysis_cut(jpsi_true, jpsi_pm_ana);
    fill_jpsi_analysis_hists(jpsi_pm_ana,pm_ana);
  }
}

void AnaTda::jpsi_analysis_cut(RhoCandList& org, RhoCandList& dest) {
  for (int j = 0; j < org.GetLength(); ++j) {
    const double m = org[j]->M();
    if ( jpsi_mcut_min < m && m < jpsi_mcut_max ) {
      dest.Append(org[j]);
    }
  }
}

void AnaTda::fill_jpsi_analysis_hists(RhoCandList& org, const int& itype) {
  fill_pair_mass(org, h_m_epem[itype]);
}

void AnaTda::Finish() {

  write_kin_fit_hists();
  write_eff_hists();
  const char* root = gDirectory->GetPath();
  gDirectory->mkdir("resid");
  gDirectory->cd("resid");
  for (int i =0; i < 4; ++i) {
    h_resid_phth[i]->Write();
    h_resid_pip_phth[i]->Write();
    h_resid_pip_mom[i]->Write();
    h_resid_pim_phth[i]->Write();
    h_resid_pim_mom[i]->Write();
    h_resid_ep_phth[i]->Write();
    h_resid_ep_mom[i]->Write();
    h_resid_em_phth[i]->Write();
    h_resid_em_mom[i]->Write();
    h_resid_pip_elec_hyp_phth[i]->Write();
    h_resid_pip_elec_hyp_mom[i]->Write();
    h_resid_pim_elec_hyp_phth[i]->Write();
    h_resid_pim_elec_hyp_mom[i]->Write();
  }
  h_m_pippim->Write();
  h_m_pippim_tr->Write();
  gDirectory->cd(root);

  gDirectory->mkdir("single");
  gDirectory->cd("single");
  h_num_g->Write();
  h_num_epm->Write();
  h_num_pipm->Write();
  h_mom_the_pipm->Write();
  h_num_g_tr->Write();
  h_num_epm_tr->Write();
  h_num_pipm_tr->Write();
  h_mom_the_pipm_tr->Write();
  for (int it=0; it<nhist; ++it) h_e_g[it]->Write();
  for (int it=0; it<nhist; ++it) h_mom_epm[it]->Write();
  gDirectory->cd(root);

  gDirectory->mkdir("man_kin_fit");
  gDirectory->cd("man_kin_fit");
  for (int it=4; it<nhist; ++it) write_manual_kin_fit_hists(it);
  gDirectory->cd(root);

  //for (int it=0; it<nhist; ++it) {
  gDirectory->mkdir("gg");
  gDirectory->cd("gg");
  for (int it=0; it<nhist; ++it) h_oa_gg[it]->Write();
  for (int it=0; it<nhist; ++it) h_oa_gg_vs_min_e_g[it]->Write();
  for (int it=0; it<nhist; ++it) h_oa_gg_vs_avg_e_g[it]->Write();
  for (int it=0; it<nhist; ++it) h_oa_gg_vs_asym_e_g[it]->Write();
  for (int it=0; it<nhist; ++it) h_m_gg[it]->Write();
  gDirectory->cd(root);

  gDirectory->mkdir("epem");
  gDirectory->cd("epem");
  for (int it=0; it<nhist; ++it) h_m_epem[it]->Write();
  for (int it=0; it<nhist; ++it) h_mom_the_epm[it]->Write();
  gDirectory->cd(root);

  gDirectory->mkdir("pi0_ana");
  gDirectory->cd("pi0_ana");
  for (int it=0; it<nhist; ++it) {
    for (int ia = 0; ia< npi0ana; ++ia) {
      if (it < 4) continue;
      h_m_gg_pi0ana[ia][it]->Write();
    }
  }
  gDirectory->cd(root);

  write_full_sys_hists();

}

void AnaTda::def_full_sys_hists() {

  const char *dth = (bg_mc?"#theta_{#gamma#gamma}+#theta_{#pi^{+}#pi^{-}}":"#theta_{#gamma#gamma}+#theta_{e^{+}e^{-}}");
  const char *dph = (bg_mc?"#phi_{#gamma#gamma}-#phi_{#pi^{+}#pi^{-}}":"#phi_{#gamma#gamma}-#phi_{e^{+}e^{-}}");
  const char *mtot = (bg_mc?"M^{inv}_{#gamma#gamma - #pi^{+}#pi^{-}}":"M^{inv}_{#gamma#gamma - e^{+}e^{-}}");
  const char *mgg = "M^{inv}_{#gamma#gamma}";
  const char *mepem = (bg_mc?"M^{inv}_{#pi^{+}#pi^{-}}":"M^{inv}_{e^{+}e^{-}}");
  //enum {all=0, tr=1, pm=2, pm_mc=3, ana=4, pm_ana=5};
  const char *n[nhist] = {"all_rec", "true_rec", "truepi0jpsi_rec", "truepi0jpsi_mc", "ana", "pm_ana", "mconst"};
  //const char *t[nhist] = {"all (Reco)", "truth match (Reco)", "from true J/#psi-#pi^{0} (Reco)",
  //     "from true J/#psi-#pi^{0} (MC)", "after analysis cuts", "from true J/#psi-#pi^{0} after ana. cut"};
  const char *t[nhist] = {"all (Reco)", "truth match (Reco)", (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} (Reco)":"from true J/#psi-#pi^{0} (Reco)"),
			  (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} (MC)":"from true J/#psi-#pi^{0} (MC)"), "after analysis cuts",
			  (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} after ana. cut":"from true J/#psi-#pi^{0} after ana. cut"),
			  (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} after mass constraint":"from true J/#psi-#pi^{0} after mass constraint"),
  };

  double _m_tot_min=iplab==0?1.5:(iplab==1?2.0:2.5);
  double _m_tot_max=iplab==0?4.5:(iplab==1?5.5:6.5);
  for (int it=2; it<nhist; ++it) {
    h_dth_vs_mass_gg_epair[it] = new
      TH2F(Form("h_dth_vs_mass_gg_epair_%s",n[it]),
	   Form("%s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta[rad]",dth,mtot,t[it], mtot),
	   200, _m_tot_min, _m_tot_max, 200, TMath::Pi()/2., 3*TMath::Pi()/2.);
    h_dth_vs_dph_gg_epair[it] = new
      TH2F(Form("h_dth_vs_dph_gg_epair_%s",n[it]),
	   Form("%s vs. %s (%s);#Delta#phi[rad];#Delta#theta[rad]",dth,dph,t[it]),
	   //200, 2.94, 3.34, 200, 2.94, 3.34);
	   200, TMath::Pi()/2., 3*TMath::Pi()/2., 200, TMath::Pi()/2., 3*TMath::Pi()/2.);
    h_dth_gg_epair_vs_mass_gg[it] = new
      TH2F(Form("h_dth_gg_epair_vs_mass_gg_%s",n[it]),
	   Form("%svs. %s (%s);%s[GeV/c^{2}];#Delta#theta[rad]",dth, mgg,t[it], mgg),
	   200, 0, 0.2, 200, TMath::Pi()/2., 3*TMath::Pi()/2.);
    h_mass_gg_epair_vs_mass_gg[it] = new
      TH2F(Form("h_mass_gg_epair_vs_mass_gg_%s",n[it]),
	   Form("%s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta[rad]",dth,mgg,t[it],mgg),
	   200, 0, 0.2, 200, _m_tot_min, _m_tot_max);
  }
}

void AnaTda::def_manual_kin_fit_hists(const int &itype ) {
  //const char *dth = "#theta_{#gamma#gamma}+#theta_{e^{+}e^{-}}";
  //const char *mtot = "M^{inv}_{#gamma#gamma - e^{+}e^{-}}";
  //const char *mgg = "M^{inv}_{#gamma#gamma}";
  //const char *mepem = "M^{inv}_{e^{+}e^{-}}";

  const char *dth = (bg_mc?"#theta_{#gamma#gamma}+#theta_{#pi^{+}#pi^{-}}":"#theta_{#gamma#gamma}+#theta_{e^{+}e^{-}}");
  const char *dph = (bg_mc?"#phi_{#gamma#gamma}-#phi_{#pi^{+}#pi^{-}}":"#phi_{#gamma#gamma}-#phi_{e^{+}e^{-}}");
  const char *mtot = (bg_mc?"M^{inv}_{#gamma#gamma - #pi^{+}#pi^{-}}":"M^{inv}_{#gamma#gamma - e^{+}e^{-}}");
  const char *mgg = "M^{inv}_{#gamma#gamma}";
  const char *mepem = (bg_mc?"M^{inv}_{#pi^{+}#pi^{-}}":"M^{inv}_{e^{+}e^{-}}");

  //enum {all=0, tr=1, pm=2, pm_mc=3, ana=4, pm_ana=5};
  const char *n[nhist] = {"all_rec", "true_rec", "truepi0jpsi_rec", "truepi0jpsi_mc", "ana", "pm_ana", "mconst"};
  //const char *t[nhist] = {"all (Reco)", "truth match (Reco)", "from true J/#psi-#pi^{0} (Reco)",
  //     "from true J/#psi-#pi^{0} (MC)", "after analysis cuts", "from true J/#psi-#pi^{0} after ana. cut"};
  const char *t[nhist] =
    {"all (Reco)", "truth match (Reco)", (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} (Reco)":"from true J/#psi-#pi^{0} (Reco)"),
     (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} (MC)":"from true J/#psi-#pi^{0} (MC)"), "after analysis cuts",
     (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} after ana. cut":"from true J/#psi-#pi^{0} after ana. cut"),
     (bg_mc?"from true #pi^{+}#pi^{-}#pi^{0} after mass constraint":"from true J/#psi-#pi^{0} after mass constraint") };

  // cts = closest-to-s, btb = most-back-to-back
  h_dth_vs_mass_gg_epair_btb[itype] = new
    TH2F(Form("h_dth_vs_mass_gg_epair_btb_%s",n[itype]),
	 Form("most back-to-back %s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta", dth, mtot, t[itype], mtot),
	 200, 2.9, 3.7, 200, TMath::Pi()/2., 3*TMath::Pi()/2.);

  h_dth_vs_mass_gg_epair_cts[itype] = new
    TH2F(Form("h_dth_vs_mass_gg_epair_cts_%s",n[itype]),
	 Form("closest to #sqrt{s} %s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta", dth, mtot, t[itype], mtot),
	 200, 2.9, 3.7, 200, TMath::Pi()/2., 3*TMath::Pi()/2.);

  h_m_gg_btb[itype] = new
    TH1F(Form("h_m_gg_btb_%s",n[itype]),
	 Form("most back-to-back %s (%s);%s[GeV/c^{2}]", mgg, t[itype], mgg),
	 100, 0, 0.2);

  h_m_gg_cts[itype] = new
    TH1F(Form("h_m_gg_cts_%s",n[itype]),
	 Form("closest to #sqrt{s} %s (%s); %s[GeV/c^{2}]", mgg, t[itype], mgg),
	 100, 0, 0.2);

  h_m_epem_btb[itype] = new
    TH1F(Form("h_m_epem_btb_%s",n[itype]),
	 Form("most back-to-back %s (%s);%s[GeV/c^{2}]", mepem, t[itype], mepem),
	 100, 0, 0.2);

  h_m_epem_cts[itype] = new
    TH1F(Form("h_m_epem_cts_%s",n[itype]),
	 Form("closest to #sqrt{s} %s (%s); %s[GeV/c^{2}]", mepem, t[itype], mepem),
	 100, 0, 0.2);

}

void AnaTda::calc_kin(
  RhoCandidate* _gg, RhoCandidate *_epem,
  double &m, double &mgg, double &dth, double &dph,
  double &thpi0, double &thep) {

  TLorentzVector p4gg = _gg->P4();

  thpi0 = p4gg.Vect().Theta();

  TLorentzVector p4epair = _epem->P4();

  TLorentzVector p4ep = _epem->Daughter(0)->Charge()>0? _epem->Daughter(0)->P4(): _epem->Daughter(1)->P4();

  m = (p4gg+p4epair).M();
  mgg = (p4gg).M();
  p4gg.Boost(boost_to_cm);
  p4epair.Boost(boost_to_cm);
  p4ep.Boost(boost_to_cm);
  dth = fabs(p4gg.Vect().Theta() + p4epair.Vect().Theta());
  dph = p4gg.Vect().Phi() - p4epair.Vect().Phi();
  if (dph<0) dph *= -1;

  thep = p4ep.Vect().Theta();

  //cout << "dph= " << dph << endl;
}

void AnaTda::calc_kin_from_daughters(
  RhoCandidate* _gg0,   RhoCandidate* _gg1,
  RhoCandidate *_epem0, RhoCandidate *_epem1,
  double &m, double &mgg, double &dth, double &dph) {
  TLorentzVector p4gg = _gg0->P4()+_gg1->P4();
  TLorentzVector p4epair = _epem0->P4()+_epem1->P4();
  m = (p4gg+p4epair).M();
  mgg = (p4gg).M();
  p4gg.Boost(boost_to_cm);
  p4epair.Boost(boost_to_cm);
  dth = fabs(p4gg.Vect().Theta() + p4epair.Vect().Theta());
  dph = p4gg.Vect().Phi() - p4epair.Vect().Phi();
  if (dph<0) dph *= -1;
  //cout << "dph= " << dph << endl;
}


void AnaTda::pi0jpsi_true_kinematics(RhoCandList& org_gg, RhoCandList& org_epem) {

  double m, mgg, dth, dph;

  if (verb) cout << "========= FILL PI0JPSI TRUE KINEMATICS HISTS ========" << endl;
  if (verb) cout << "org_epem.GetLength= " << org_epem.GetLength() << "org_gg.GetLength= " << org_gg.GetLength() << endl;
  if (org_epem.GetLength()!=1 or org_gg.GetLength()!=1) return;

  if (verb){
    cout << "TrackNum(pi0): M= " << org_gg[0]->GetTrackNumber()
	 << " d0= " << org_gg[0]->Daughter(0)->GetTrackNumber() << " d1= " << org_gg[0]->Daughter(1)->GetTrackNumber() << endl;
    cout << "TrackNum(psi): M= " << org_epem[0]->GetTrackNumber()
	 << " d0= " << org_epem[0]->Daughter(0)->GetTrackNumber() << " d1= " << org_epem[0]->Daughter(1)->GetTrackNumber() << endl;
  }

  //calc_kin(org_gg[0], org_epem[0], m, mgg, dth, dph);
  calc_kin_from_daughters(org_gg[0]->Daughter(0), org_gg[0]->Daughter(1),
			  org_epem[0]->Daughter(0), org_epem[0]->Daughter(1), m, mgg, dth, dph);

  h_dth_gg_epair_vs_mass_gg[pm]->Fill(mgg, dth);
  h_mass_gg_epair_vs_mass_gg[pm]->Fill(mgg, m);

  h_dth_vs_mass_gg_epair[pm]->Fill(m, dth);
  h_dth_vs_dph_gg_epair[pm]->Fill(dph, dth);

  if (verb) {
    cout << "pi0: ";
    org_gg[0]->PrintOn(cout);
    cout<<endl;
    cout << "psi: ";
    org_epem[0]->PrintOn(cout);
    cout<<endl;
    cout << "pi0_mc: ";
    org_gg[0]->GetMcTruth()->PrintOn(cout);
    cout<<endl;
    cout << "psi_mc: ";
    org_epem[0]->GetMcTruth()->PrintOn(cout);
    cout<<endl;

    cout << "TruthMatch(pi0): M= " << (org_gg[0]->GetMcTruth()?org_gg[0]->GetMcTruth()->GetTrackNumber(): -9999)
	 << " d0= " << (org_gg[0]->Daughter(0)->GetMcTruth()?org_gg[0]->Daughter(0)->GetMcTruth()->GetTrackNumber(): -9999)
	 << " d1= " << (org_gg[0]->Daughter(1)->GetMcTruth()?org_gg[0]->Daughter(1)->GetMcTruth()->GetTrackNumber(): -9999) << endl;
    cout << "TruthMatch(psi): M= " << (org_epem[0]->GetMcTruth()?org_epem[0]->GetMcTruth()->GetTrackNumber(): -9999)
	 << " d0= " << (org_epem[0]->Daughter(0)->GetMcTruth()?org_epem[0]->Daughter(0)->GetMcTruth()->GetTrackNumber(): -9999)
	 << " d1= " << (org_epem[0]->Daughter(1)->GetMcTruth()?org_epem[0]->Daughter(1)->GetMcTruth()->GetTrackNumber(): -9999) << endl;
  }

  calc_kin_from_daughters(org_gg[0]->Daughter(0)->GetMcTruth(), org_gg[0]->Daughter(1)->GetMcTruth(),
			  org_epem[0]->Daughter(0)->GetMcTruth(), org_epem[0]->Daughter(1)->GetMcTruth(), m, mgg, dth, dph);
  h_dth_gg_epair_vs_mass_gg[pm_mc]->Fill(mgg, dth);
  h_mass_gg_epair_vs_mass_gg[pm_mc]->Fill(mgg, m);
  h_dth_vs_mass_gg_epair[pm_mc]->Fill(m, dth);
  h_dth_vs_dph_gg_epair[pm_mc]->Fill(dph, dth);
}

inline
bool AnaTda::passes_kin_cut(const double &m, const double &dth) {
  return (  m > etot_min  and  m < etot_max  and dth > dth_min and dth < dth_max );
}

void AnaTda::pi0jpsi_efficiency(RhoCandList& org_gg, RhoCandList& org_epem, const int& tt) {

  // Start with pi0_ana and jpsi_ana

  // pi0->All gg pairs,  jpsi->All ep/em or pip/pim pairs
  //REF: pi0_true->All primary matched gg pairs,  jpsi_true->All primary matched ep/em or pip/pim pairs

  // pi0_ana->gg pairs with ana cuts,  jpsi_ana->All ep/em or pip/pim pairs with ana cuts
  //REF: pi0_pm_ana->primary matched gg pairs with ana cuts,  jpsi_pm_ana->Primary matched ep/em or pip/pim pairs with ana cuts

  // If finding was ambiguous, skip event
  if (org_gg.GetLength()!=1) return;
  if (org_epem.GetLength()!=1) return;

  double m, mgg, dth, dph, thpi0, thep;
  int nexcl = 0;
  h_num_gg[tt]->Fill(org_gg.GetLength());
  h_num_epem[tt]->Fill(org_epem.GetLength());
  for (int i = 0; i < org_epem.GetLength(); ++i) {
    for (int j = 0; j < org_gg.GetLength(); ++j) {
      if (tt == mconst) {
	calc_kin(org_gg[j], org_epem[i]->GetFit(), m, mgg, dth, dph, thpi0, thep);
	if (passes_kin_cut(m,dth)){
	  h_eff_thpi0[tt]->Fill(thpi0);
	  h_eff_thep[tt]->Fill(thep);
	}
      } else {
	calc_kin(org_gg[j], org_epem[i], m, mgg, dth, dph, thpi0, thep);
	if (tt==eff_kin || tt==eff_excl) {
	  if (passes_kin_cut(m,dth)){
	    if (tt==eff_kin){
	      h_eff_thpi0[tt]->Fill(thpi0);
	      h_eff_thep[tt]->Fill(thep);
	    }
	    nexcl++;
	  }
	} else {
	  h_eff_thpi0[tt]->Fill(thpi0);
	  h_eff_thep[tt]->Fill(thep);
	}
      }
    }
  }

  if (tt==eff_excl && nexcl==1) {
    for (int i = 0; i < org_epem.GetLength(); ++i) {
      for (int j = 0; j < org_gg.GetLength(); ++j) {
	calc_kin(org_gg[j], org_epem[i], m, mgg, dth, dph, thpi0, thep);
	if (passes_kin_cut (m, dth)) {
	  h_eff_thpi0[tt]->Fill(thpi0);
	  h_eff_thep[tt]->Fill(thep);
	}
      }
    }
  }
}

void AnaTda::def_eff_hists() {
  //enum {eff_ref, eff_pi0sel, eff_jpsisel, eff_kin, eff_excl};
  const char *t[] = {"Reference","#pi^{0} selection","J/#psi selection",
		     Form("Kin. Cut (%3.1f<m_{tot}<%3.1f && %3.1f<#Delta(#theta)<%3.1f)", etot_min,etot_max,dth_min,dth_max),
		     "Exclusive", "Kin. Cut After Mass Const."};
  for (int i=0; i<neffhist; ++i) {
    h_eff_thpi0[i] = new TH1F(Form("h_eff_thpi0_%d", i),Form("#theta^{lab}_{#pi^{0}} %s",t[i]),200,0,TMath::Pi());
    h_eff_thep[i] = new TH1F(Form("h_eff_thep_%d", i),Form("#theta^{CM}_{e^{+}} %s",t[i]),200,0,TMath::Pi());
    h_num_gg[i] = new TH1F(Form("h_num_gg_%d",i), Form("Num of #gamma#gamma pairs %s", t[i]), 10,0,10);
    h_num_epem[i] = new TH1F(Form("h_num_epem_%d",i), Form("Num of charged track pairs %s", t[i]), 10,0,10);
  }
}

void AnaTda::write_eff_hists(){
  const char* root = gDirectory->GetPath();
  gDirectory->mkdir("eff");
  gDirectory->cd("eff");
  for (int i=0; i<neffhist; ++i) {
    h_eff_thpi0[i]->Write();
    h_eff_thep[i]->Write();
    h_num_gg[i]->Write();
    h_num_epem[i]->Write();
  }
  gDirectory->cd(root);
}

void AnaTda::pi0jpsi_kin_fit(RhoCandList& org_gg, RhoCandList& org_epem) {
  RhoCandList full_sys;
  full_sys.Combine(org_gg,org_epem);
  for (int i=0; i<org_epem.GetLength(); ++i) {
    PndKinFitter fitter(org_epem[i]);
    //fitter.Add4MomConstraint(ini);
    fitter.AddMassConstraint(3.096);
    fitter.Fit();
    h_m_epem[mconst]->Fill(org_epem[i]->GetFit()->M());
    jpsi_mconst.Append(org_epem[i]->GetFit());
  }
  pi0jpsi_kinematics(org_gg, jpsi_mconst, mconst);

  //double m, mgg, dth, dph, thpi0, thep;
  //for (int i = 0; i < org_epem.GetLength(); ++i) {
  //  for (int j = 0; j < org_gg.GetLength(); ++j) {
  //    calc_kin(org_gg[j], org_epem[i]->GetFit(), m, mgg, dth, dph, thpi0, thep);
  //    h_dth_vs_mass_gg_epair[mconst]->Fill(m, dth);
  //    h_dth_vs_dph_gg_epair[mconst]->Fill(dph, dth);
  //  }
  //}

}

void AnaTda::pi0jpsi_kinematics(RhoCandList& org_gg, RhoCandList& org_epem, const int& tt) {
  double m, mgg, dth, dph, thpi0, thep;
  for (int i = 0; i < org_epem.GetLength(); ++i) {
    for (int j = 0; j < org_gg.GetLength(); ++j) {
      calc_kin(org_gg[j], org_epem[i], m, mgg, dth, dph, thpi0, thep);
      h_dth_gg_epair_vs_mass_gg[tt]->Fill(mgg, dth);
      h_mass_gg_epair_vs_mass_gg[tt]->Fill(mgg, m);
      h_dth_vs_mass_gg_epair[tt]->Fill(m, dth);
      h_dth_vs_dph_gg_epair[tt]->Fill(dph, dth);
    }
  }
}

void AnaTda::pi0_kinematic_selection(RhoCandList& org_gg, RhoCandList& org_epem, const int& tt) {
  // pi0 selection based on kinematics (opposite to lepton/pion
  // pair in CM, blieveable invariant mass - close to sqrt(s))
  // start with the list of all gg pairs and plot cm-OA wrt
  // lepton pair, and inv-mass of full system

  // cts = closest-to-s, btb = most-back-to-back
  int i_btb_e = -1, i_cts_e = -1;
  int j_btb_g = -1, j_cts_g = -1;
  double diff_2pi_min = 1e9, diff_s_min = 1e9;
  double m, mgg, dth, dph, thpi0, thep;

  for (int i = 0; i < org_epem.GetLength(); ++i) {
    for (int j = 0; j < org_gg.GetLength(); ++j) {
      calc_kin(org_gg[j], org_epem[i], m, mgg, dth, dph, thpi0, thep);
      const double diff_2pi = TMath::Pi() - dth;
      const double diff_s = fabs(sqrt_s - m);
      if (diff_2pi < diff_2pi_min) {
        i_btb_e = i;
        j_btb_g = j;
        diff_2pi_min = diff_2pi;
      }
      if (diff_s < diff_s_min) {
        i_cts_e = i;
        j_cts_g = j;
        diff_s_min = diff_s;
      }
    }
  }

  // now that we have closest pair in OA and M to what we want, fill
  if (i_btb_e >= 0 and j_btb_g >= 0 ) {
    calc_kin(org_gg[j_btb_g], org_epem[i_btb_e], m, mgg, dth, dph, thpi0, thep);
    h_dth_vs_mass_gg_epair_btb[tt]->Fill(m, dth);
    h_m_gg_btb[tt]->Fill(org_gg[j_btb_g]->M());
    h_m_epem_btb[tt]->Fill(org_gg[j_btb_g]->M());
    //pi0_btb.Append(gg[i_btb]);
  }

  if (i_cts_e >= 0 and j_cts_g >= 0) {
    calc_kin(gg[j_cts_g], org_epem[i_cts_e], m, mgg, dth, dph, thpi0, thep);
    h_dth_vs_mass_gg_epair_cts[tt]->Fill(m, dth);
    h_m_gg_cts[tt]->Fill(gg[j_cts_g]->M());
    h_m_epem_cts[tt]->Fill(gg[j_cts_g]->M());
      //pi0_cts.Append(gg[i_cts]);
  }

}

void AnaTda::write_full_sys_hists() {
  const char* root = gDirectory->GetPath();
  gDirectory->mkdir("full_sys");
  gDirectory->cd("full_sys");
  for (int it=2; it<nhist; ++it) {
    h_dth_vs_mass_gg_epair[it]->Write();
    h_dth_vs_dph_gg_epair[it]->Write();
    h_dth_gg_epair_vs_mass_gg[it]->Write();
    h_mass_gg_epair_vs_mass_gg[it]->Write();
  }
  gDirectory->cd(root);
}

void AnaTda::write_manual_kin_fit_hists(const int& itype){
  h_dth_vs_mass_gg_epair_btb[itype]->Write();      // btb = most-back-to-back
  h_dth_vs_mass_gg_epair_cts[itype]->Write();      // cts = closest-to-s
  h_m_gg_btb[itype]->Write();     // btb = most-back-to-back
  h_m_gg_cts[itype]->Write();     // cts = closest-to-s
  h_m_epem_btb[itype]->Write();     // btb = most-back-to-back
  h_m_epem_cts[itype]->Write();     // cts = closest-to-s
}

// --- below are some old bad ideas
void AnaTda::def_kin_fit_hists(const int& itype, const int& hyp) {

  h_m_epem_bef_mass_fit = new TH1F("h_m_epem_bef_mass_fit","h_m_epem_bef_mass_fit",200,0,5);
  h_m_epem_aft_mass_fit = new TH1F("h_m_epem_aft_mass_fit","h_m_epem_aft_mass_fit",200,0,5);

  //const char *n[] = {"","_btb","_cts"};
  //const char *nhyp[] = {"epem", "pippim"};
  //const char *thyp[] = {"e^{+}e^{-}", "#pi^{+}#pi^{-}"};
  //
  //h_4c_chi2[itype][hyp] = new TH1F(
  //  Form("h_4c_chi2%s_%spi0", n[itype], nhyp[hyp]),
  //  Form("#chi^{2} of 4C fit on %s#pi^{0} system; #chi^{2}", thyp[hyp]),
  //  100, 0, 200);
  //
  //h_4c_prob[itype][hyp] = new TH1F(
  //  Form("h_4c_prob%s_%spi0", n[itype], nhyp[hyp]),
  //  Form("Probability of #chi^{2} of 4C fit on %s#pi^{0} system; prob", thyp[hyp]),
  //  100, -0.1, 1.1);
  //
  //h_4c_m[itype][hyp] = new TH1F(
  //  Form("h_4c_m%s_%s_4c", n[itype], nhyp[hyp]),
  //  Form("Invariant mass of %s After 4C Fit;M_{%s}[GeV/c^{2}]", thyp[hyp], thyp[hyp]),
  //  100, 0, 4.5);
  //
  //h_4c_prob_vs_m[itype][hyp] = new TH2F(
  //  Form("h_4c_prob_vs_m%s_%s_4c", n[itype], nhyp[hyp]),
  //  Form("Prob. vs. Invariant mass of %s After 4C Fit;M_{%s}[GeV/c^{2}]; prob", thyp[hyp], thyp[hyp]),
  //  100, 2.9, 3.7, 100, 0, 1);
}

void AnaTda::def_kin_fit_hists() {
  for (int it = 0; it < 3; ++it)
    for (int ih = 0; ih < 2; ++ih)
      def_kin_fit_hists(it, ih);
}

void AnaTda::write_kin_fit_hists() {

  h_m_epem_bef_mass_fit->Write();
  h_m_epem_aft_mass_fit->Write();

  //for (int it = 0; it < 3; ++it) {
  //  for (int ih = 0; ih < 2; ++ih) {
  //    h_4c_chi2[it][ih]->Write();
  //    h_4c_prob[it][ih]->Write();
  //    h_4c_m[it][ih]->Write();
  //    h_4c_prob_vs_m[it][ih]->Write();
  //  }
  //}
}

void AnaTda::kin_fit_full_sys(RhoCandList& org, const int& it, const int& ih) {
  if (org.GetLength() == 1) {
    PndKinFitter fitter(org[0]);
    fitter.Add4MomConstraint(ini);
    fitter.Fit();
    const double chi2 = fitter.GetChi2();
    const double prob = fitter.GetProb();
    h_4c_chi2[it][ih]->Fill(chi2);
    h_4c_prob[it][ih]->Fill(prob);
    h_4c_prob_vs_m[it][ih]->Fill(org[0]->M(), prob);
  } else {
    if (org.GetLength() != 0) {
      cout << "Not possible btb" << endl;
    }
  }
}

void AnaTda::kin_fit_all() {
  kin_fit_epem_pi0_btb();
  kin_fit_epem_pi0_cts();
  kin_fit_pippim_pi0_btb();
  kin_fit_pippim_pi0_cts();
}


void AnaTda::kin_fit_epem_pi0_btb() {
  epem_mcut_pi0_btb.Combine(epem_mcut, pi0_btb);
  kin_fit_full_sys(epem_mcut_pi0_btb, 1, 0);
}

void AnaTda::kin_fit_epem_pi0_cts() {
  epem_mcut_pi0_cts.Combine(epem_mcut, pi0_cts);
  kin_fit_full_sys(epem_mcut_pi0_cts, 2, 0);
}

void AnaTda::kin_fit_pippim_pi0_btb() {
  pippim_mcut_pi0_btb.Combine(pippim_mcut, pi0_btb);
  kin_fit_full_sys(pippim_mcut_pi0_btb, 1, 1);
}

void AnaTda::kin_fit_pippim_pi0_cts() {
  pippim_mcut_pi0_cts.Combine(pippim_mcut, pi0_cts);
  kin_fit_full_sys(pippim_mcut_pi0_cts, 2, 1);
}

void AnaTda::def_tutorial_hists() {
  // *** create some histograms
  hjpsim_all = new TH1F("hjpsim_all", "J/#psi mass (all)", 100, 0, 4.5);
  hpi0m_all  = new TH1F("hpi0m_all", "#pi^{0} mass (all)", 100, 0, 0.2);
  hjpsipi0m_all  = new TH1F("hjpsipi0m_all", "J/#psi-#pi^{0} mass (all)", 100, 0, 4.5);

  hjpsim_ftm = new TH1F("hjpsim_ftm", "J/#psi mass (full truth match)", 100, 0, 4.5);
  hpi0m_ftm  = new TH1F("hpi0m_ftm", "#pi^{0} mass (full truth match)", 100, 0, 0.2);
  hjpsipi0m_ftm  = new TH1F("hjpsipi0m_ftm", "J/#psi-#pi^{0} mass (full truth match)", 100, 0, 4.5);

  hjpsim_nm = new TH1F("hjpsim_nm", "J/#psi mass (no truth match)", 100, 0, 4.5);
  hpi0m_nm  = new TH1F("hpi0m_nm", "#pi^{0} mass (no truth match)", 100, 0, 0.2);
  hjpsipi0m_nm  = new TH1F("hjpsipi0m_nm", "J/#psi-#pi^{0} mass (no truth match)", 100, 0, 4.5);

  hjpsim_diff = new TH1F("hjpsim_diff", "J/#psi mass diff to truth", 100, -2, 2);
  hpi0m_diff  = new TH1F("hpi0m_diff", "#pi^{0} mass diff to truth", 100, -2, 2);
  hjpsipi0m_diff  = new TH1F("hjpsipi0m_diff", "J/#psi-#pi^{0} mass diff to truth", 100, -2, 2);
}

void AnaTda::write_tut_hists() {
  hjpsim_all->Write();
  hpi0m_all->Write();
  hjpsipi0m_all->Write();

  hjpsim_ftm->Write();
  hpi0m_ftm->Write();
  hjpsipi0m_ftm->Write();

  hjpsim_nm->Write();
  hpi0m_nm->Write();
  hjpsipi0m_nm->Write();

  hjpsim_diff->Write();
  hpi0m_diff->Write();
  hjpsipi0m_diff->Write();
}

void AnaTda::kin_fit_pi0_nearest(RhoCandList& org, const int& itype) {
  for (int j = 0; j < org.GetLength(); ++j) {
    //if (itype==0)
    //  h_m_epem_pi0n->Fill(org[j]->M());
    //else
    //  h_m_pippim_pi0n->Fill(org[j]->M());
    PndKinFitter fitter(org[j]);  // instantiate the kin fitter
    fitter.Add4MomConstraint(ini);  //  set 4 constraint
    fitter.Fit();  //  do fit
    double chi2_4c = fitter.GetChi2();  //  get chi2 of fit
    double prob_4c = fitter.GetProb();  //  access probability of fit
    h_4c_chi2[0][itype]->Fill(chi2_4c);
    h_4c_prob[0][itype]->Fill(prob_4c);
    if ( prob_4c > 0.01 ) {
      RhoCandidate *_fit = org[j]->Daughter(0)->GetFit();
      h_4c_m[0][itype]->Fill(_fit->M());
    }
  }
}

void AnaTda::kin_fit_pi0_nearest_all() {
  kin_fit_epem_pi0_nearest();
  kin_fit_pippim_pi0_nearest();
}

void AnaTda::kin_fit_epem_pi0_nearest() {
  epem_pi0nearest.Combine(epem, pi0nearest);
  kin_fit_pi0_nearest(epem_pi0nearest, 0);
}

void AnaTda::kin_fit_pippim_pi0_nearest() {
  pippim_pi0nearest.Combine(pippim, pi0nearest);
  kin_fit_pi0_nearest(pippim_pi0nearest, 1);
}

void AnaTda::pdgm_nearest_pi0s() {
  double dm_min = 1e10;
  int min_j = -1;
  for (int j = 0; j < pi0.GetLength(); ++j) {
    const double dm = fabs(pi0[j]->M() - m0_pi0);
    if ( dm < dm_min ) {
      dm_min = dm;
      min_j = j;
    }
  }
  if (min_j >= 0) {
    pi0nearest.Append(pi0[min_j]);
  }
}

//void AnaTda::jpsi_truth_match() {
//  // following is relevant only for signal
//  jpsi.Combine(ep, em);
//  jpsi.SetType(pdg_jpsi);
//  for (int j = 0; j < jpsi.GetLength(); ++j) {
//    hjpsim_all->Fill(jpsi[j]->M());
//    if (fAnalysis->McTruthMatch(jpsi[j])) {
//      hjpsim_ftm->Fill(jpsi[j]->M());
//      hjpsim_diff->Fill(jpsi[j]->GetMcTruth()->M() - jpsi[j]->M());
//      jpsi_true.Append(jpsi[j]);
//    } else {
//     hjpsim_nm->Fill(jpsi[j]->M());
//    }
//  }
//}

//void AnaTda::pi0_truth_match() {
//  pi0.Combine(g, g);
//  pi0.SetType(pdg_pi0);
//  for (int j = 0; j < pi0.GetLength(); ++j) {
//    hpi0m_all->Fill(pi0[j]->M() );
//    double edist = 1e9;
//    if ( primary_match_pi0(pi0[j], edist) ) {
//        hpi0m_ftm->Fill( pi0[j]->M() );
//        hpi0m_diff->Fill(pi0[j]->GetMcTruth()->M() - pi0[j]->M() );
//        pi0_true.Append(pi0[j]);
//        fill_gamma_from_pi0s(pi0[j]);
//    } else {
//      //cout << "no match :( mgg = " << mgg << " mpi0 = " << mpi0 << endl;
//       hpi0m_nm->Fill(pi0[j]->M() );
//    }
//  }
//}



ClassImp(AnaTda)

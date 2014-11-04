// The header file
#include "AnaTda.h"

// C++ headers
#include <vector>
#include <string>
#include <iostream>

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

AnaTda::AnaTda(const int &brem)
  :FairTask("Radiation Length Profiler"), nevt(0),
  //pi0_ana_(npi0ana), pi0_pm_ana_(npi0ana),
  pi0oacut(npi0ana), pi0ecut_min(npi0ana), pi0ecut_max(npi0ana) {
  fBremCorr = (brem != 0);

  pi0oacut[0] = pi0oacut[1] = pi0oacut[2] = pi0oacut[3] = 0.2; //
  pi0oacut[4] = pi0oacut[5] = pi0oacut[6] = pi0oacut[7] = 0.2;

  pi0mcut_min[0] = 0.005;
  pi0mcut_min[1] = 0.01;
  pi0mcut_min[2] = 0.015;
  pi0mcut_min[3] = 0.02;
  pi0mcut_min[4] = 0.025;
  pi0mcut_min[5] = 0.05;
  pi0mcut_min[6] = 0.075;
  pi0mcut_min[7] = 0.1;

  pi0ecut_max[0] = 1e6;
  pi0ecut_max[1] = 1e6;
  pi0ecut_max[2] = 1e6;
  pi0ecut_max[3] = 1e6;
  pi0ecut_max[4] = 1e6;
  pi0ecut_max[5] = 1e6;
  pi0ecut_max[6] = 1e6;
  pi0ecut_max[7] = 1e6;

  pi0mcut_min = 0.1;
  pi0mcut_max = 0.16;

  //jpsi_mcut_min = 2.96;
  //jpsi_mcut_max = 3.22;
  jpsi_mcut_min = 2.5;
  jpsi_mcut_max = 3.7;
}

InitStatus AnaTda::Init() {
  fAnalysis = new PndAnalysis();
  def_hists();
  set_selectors();
  initial_state();
  cout << "AnaTda::Init done" << endl;
  return kSUCCESS;
}

void AnaTda::def_hists() {
  def_gamma_from_pi0_hists();
  def_elecs_from_jpsi_hists();
  def_tutorial_hists();
  def_pair_hists();
  def_single_hists();
  //def_kin_fit_hists();
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
  ini = TLorentzVector(0, 0, 5.513, 6.53023); //(0, 0, 6.231552, 7.240065);
  const double mass_prot= 0.938;
  const double p_antip = 5.513;
  const double E_antip = TMath::Hypot(mass_prot, p_antip);
  const double beta_cm = p_antip/(E_antip + mass_prot);
  cout << "betac_cm = " << beta_cm << endl;
  boost_to_cm.SetZ(-beta_cm);
  boost_to_cm.Print();
  sqrt_s = TMath::Sqrt(2*mass_prot*(mass_prot+E_antip));
  cout << "E_cm = " << sqrt_s << endl;
}

void AnaTda::Exec(Option_t* opt) {
  if (++nevt%100 == 0)
    cout << "===== AnaTda::Exec -- Event " << nevt << " ====="<< endl;

  cleanup_rho_cand_lists();

  fAnalysis->GetEvent();

  get_singles_lists();

  print_mc_list();

  fill_single_dists();
  truth_match_singles();
  fill_single_dists_tr();

  make_pair_lists();
  fill_pair_dists();

  pi0_truth_match();
  pi0_analysis_cut(); // this preps pi0_ana

  //for (int i=0; i<g1.GetLength(); ++i) {
  //  cout << "g1.GetTrackNumber= " << g1[i]->GetTrackNumber();
  //  if (g1[i]->GetMcTruth() != 0) {
  //    cout << "  -> GetMcTruth.GetTrackNumber= " << g1[i]->GetMcTruth()->GetTrackNumber();
  //  }
  //  cout << endl;
  //}
  //for (int i=0; i<g2.GetLength(); ++i) {
  //  cout << "g2.GetTrackNumber= " << g2[i]->GetTrackNumber();
  //  if (g2[i]->GetMcTruth() != 0) {
  //    cout << "  -> GetMcTruth.GetTrackNumber= " << g2[i]->GetMcTruth()->GetTrackNumber();
  //  }
  //  cout << endl;
  //}
  //return;

  jpsi_truth_match();
  jpsi_analysis_cut();

  all_ana.Combine(pi0_ana,jpsi_ana);
  all_pm_ana.Combine(pi0_pm_ana,jpsi_pm_ana);

  //cout << "pi0_pm_ana.Length= " << pi0_pm_ana.GetLength()
  //     << "jpsi_pm_ana.Length= " << jpsi_pm_ana.GetLength()
  //     << "all_pm_ana.Length= " << all_pm_ana.GetLength() << endl;

  //cout << "pi0_ana.Length= " << pi0_ana.GetLength()
  //     << "jpsi_ana.Length= " << jpsi_ana.GetLength()
  //     << "all_ana.Length= " << all_ana.GetLength() << endl;

  pi0_kinematic_selection(pi0_ana,jpsi_ana,ana);
  pi0_kinematic_selection(pi0_pm_ana,jpsi_pm_ana,pm_ana);

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

void AnaTda::def_gamma_from_pi0_hists() {
  const char *n[nhist] = {"all_rec", "true_rec", "truepi0_rec", "truepi0_mc", "ana", "pm_ana"};
  const char *t[nhist] = {"all (Reco)", "truth match (Reco)", "from true #pi^{0} (Reco)",
     "from true #pi^{0} (MC)", "after analysis cuts", "from true #pi^{0} after ana. cut"};
  for (int it = 0; it < nhist; ++it) {
    h_oa_gg[it] = new TH1F(Form("h_oa_gg_%s", n[it]),
      Form("Opening angle of #gamma-#gamma pairs %s;OA[rad]", t[it]),
      100, 0, TMath::Pi()/2.);
    h_m_gg[it] = new TH1F(Form("h_m_gg_%s", n[it]),
      Form("Mass of #gamma-#gamma pairs %s;M_{#gamma#gamma}[GeV/c^{2}]", t[it]),
      100, 0, 0.25);
    for (int ia = 0; ia < npi0ana; ++ia) {
      if (it < 4) continue;
      h_m_gg_pi0ana[ia][it] = new TH1F(Form("h_m_gg_%s_%d", n[it], ia),
        Form("Mass of #gamma-#gamma pairs %s (Ana. Cut %d);M_{#gamma#gamma}[GeV/c^{2}]", t[it], ia),
        100, 0, 0.4);
    }
    h_e_g[it] = new TH1F(Form("h_e_g_%s", n[it]),
      Form("Energy of #gamma %s;E[GeV]", t[it]),
      100, 0, 1.6);
  }
}

void AnaTda::def_elecs_from_jpsi_hists() {
  const char *n[nhist] = {"all_rec", "true_rec", "truejpsi_rec", "truejpsi_mc", "ana", "pm_ana"};
  const char *t[nhist] = {"all (Reco)", "truth match (Reco)", "from true J/#psi (Reco)",
     "from true J/#psi (MC)", "after analysis cuts", "from true J/#psi after ana. cut"};
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
  g1.Cleanup();
  g2.Cleanup();

  pip_tr.Cleanup();
  pim_tr.Cleanup();
  ep_tr.Cleanup();
  em_tr.Cleanup();
  g1_tr.Cleanup();
  g2_tr.Cleanup();

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

  all_ana.Cleanup();
  all_pm_ana.Cleanup();

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
  for (int j=0;j<mcList.GetLength();++j) {
    RhoCandidate *mcmother = mcList[j]->TheMother();
    int muid = -1;
    if (mcmother) muid = mcmother->GetTrackNumber();
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
  fAnalysis->FillList(g1, "Neutral");
  fAnalysis->FillList(g2, "Neutral");
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
  h_num_g->Fill(g1.GetLength());
  h_num_epm->Fill(ep.GetLength()+em.GetLength());
  h_num_pipm->Fill(pip.GetLength()+pim.GetLength());

  fill_single_e(g1, g2, h_e_g[all]);
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
  truth_match(pip, pip_tr, 211);
  truth_match(pim, pim_tr, -211);
  truth_match(g1, g1_tr, 22);
  truth_match(g2, g2_tr, 22);
}

void AnaTda::fill_single_dists_tr() {
  // Single distributions (After truth match)
  h_num_g_tr->Fill(g1_tr.GetLength());
  h_num_epm_tr->Fill(ep_tr.GetLength()+em_tr.GetLength());
  h_num_pipm_tr->Fill(pip_tr.GetLength()+pim_tr.GetLength());

  fill_single_e(g1, g2, h_e_g[tr]);
  fill_single_mom(ep, em, h_mom_epm[tr]);
  fill_single_mom_the(ep, em, h_mom_the_epm[tr]);

  for (int i = 0; i < pip_tr.GetLength(); ++i)
    h_mom_the_pipm_tr->Fill(pip_tr[i]->P3().Mag(), pip_tr[i]->P3().Theta() );
  for (int i = 0; i < pim_tr.GetLength(); ++i)
    h_mom_the_pipm_tr->Fill(pim_tr[i]->P3().Mag(), pim_tr[i]->P3().Theta() );
}

void AnaTda::make_pair_lists() {
  // candidate lists for pairwise truth match (full hierarchy)
  epem.Combine(ep, em);
  epem.SetType(443);
  epem_tr.Combine(ep_tr, em_tr);
  //pippim.Combine(pip, pim);
  //pippim_tr.Combine(pip_tr, pim_tr);
  gg.Combine(g1, g2);
  gg_tr.Combine(g1_tr, g2_tr);
}

void AnaTda::fill_pair_mass(RhoCandList& org, TH1F* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->M());
}

void AnaTda::fill_pair_oa(RhoCandList& org, TH1F* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(oa(org[j]->Daughter(0),org[j]->Daughter(1)));
}

void AnaTda::fill_pair_dists() {
  fill_pair_mass(epem, h_m_epem[all]);
  fill_pair_mass(epem_tr, h_m_epem[tr]);
  //fill_pair_mass(pippim, h_m_pippim);
  //fill_pair_mass(pippim_tr, h_m_pippim_tr);
  fill_pair_mass(gg, h_m_gg[all]);
  fill_pair_mass(gg_tr, h_m_gg[tr]);
  fill_pair_oa(gg, h_oa_gg[all]);
  fill_pair_oa(gg_tr, h_oa_gg[tr]);
}

inline
double AnaTda::oa(RhoCandidate* c1, RhoCandidate* c2) {
  return c1->P3().Angle(c2->P3());
}

inline
double AnaTda::mass(RhoCandidate* c1, RhoCandidate* c2) {
  return (c1->P4() + c2->P4()).M();
}

void AnaTda::primary_match(const int& pdgm, RhoCandidate* r, RhoCandidate* m, RhoCandidate *d1, RhoCandidate *d2) {
  for (int j=0;j<mcList.GetLength();++j) {
    RhoCandidate *mcmother = mcList[j]->TheMother();
    if (!mcmother) {
      *r = *mcList[j];  // deep copy required here!
    } else if ( r and mcmother->GetTrackNumber() == r->GetTrackNumber() ) {
      if ( mcList[j]->PdgCode() == pdgm ) {
        *m = *mcList[j];
        *d1 = *mcList[j]->Daughter(0);
        *d2 = *mcList[j]->Daughter(1);
        break;
      }
    }
  }
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

bool AnaTda::primary_match_pair(const int& pdgm, RhoCandidate* c, double &dist) {
  c->SetType(pdgm); // just to make sure in case it was forgotten
  if (!fAnalysis->McTruthMatch(c)) return false;
  RhoCandidate _root, _m, _d1, _d2;
  primary_match(pdgm,&_root,&_m,&_d1,&_d2);
  if (_root.GetTrackNumber() != 0 or _m.GetTrackNumber()==-1
    or _d1.GetTrackNumber()==-1 or _d2.GetTrackNumber()==-1) return false;
  dist = hypot(c->Daughter(0)->Energy()-_d1.Energy(), c->Daughter(1)->Energy()-_d2.Energy());
  return ( (c->GetMcTruth()->GetTrackNumber() == _m.GetTrackNumber())
    and (c->Daughter(0)->GetMcTruth()->GetTrackNumber() == _d1.GetTrackNumber())
    and (c->Daughter(1)->GetMcTruth()->GetTrackNumber() == _d2.GetTrackNumber()));
}

bool AnaTda::primary_match_pi0(RhoCandidate *c, double &dist) {
  return primary_match_pair(111, c, dist);
}

void AnaTda::pi0_truth_match() {
  RhoCandidate *_match;
  double dist_min = 1e9;
  pi0.Combine(g1, g2);
  pi0.SetType(111);
  for (int j = 0; j < pi0.GetLength(); ++j) {
    double dist = 0;
    if ( primary_match_pi0(pi0[j], dist) ) {
      if (dist < dist_min) {
        _match = pi0[j];
        dist_min = dist;
      }
    }
  }
  if (dist_min<1e8) {
    pi0_true.Append(_match);
    fill_gamma_from_pi0s();
  }
}

void AnaTda::fill_gamma_from_pi0s() {
  if (pi0_true.GetLength()!=1) return;
  RhoCandidate* _g1 = pi0_true[0]->Daughter(0);
  RhoCandidate* _g2 = pi0_true[0]->Daughter(1);
  RhoCandidate* _g1_mc = pi0_true[0]->Daughter(0)->GetMcTruth();
  RhoCandidate* _g2_mc = pi0_true[0]->Daughter(1)->GetMcTruth();
  h_oa_gg[pm]->Fill(oa(_g1,_g2));
  h_oa_gg[pm_mc]->Fill(oa(_g1_mc,_g2_mc));
  h_m_gg[pm]->Fill(mass(_g1,_g2));
  h_m_gg[pm_mc]->Fill(mass(_g1_mc,_g2_mc));
  h_e_g[pm]->Fill(_g1->Energy());
  h_e_g[pm]->Fill(_g2->Energy());
  h_e_g[pm_mc]->Fill(_g1_mc->Energy());
  h_e_g[pm_mc]->Fill(_g2_mc->Energy());
}

void AnaTda::pi0_analysis_cut() {
  pi0_analysis_cut(pi0, pi0_ana, 0, true);
  fill_pi0_analysis_hists(pi0_ana, ana);
  pi0_analysis_cut(pi0_true, pi0_pm_ana, 0, true);
  fill_pi0_analysis_hists(pi0_pm_ana, pm_ana);
  for (int i=0; i < npi0ana; ++i) {
    pi0_analysis_cut(pi0, pi0_ana_[i], i%(npi0ana/2), i>=npi0ana/2);
    fill_pi0_analysis_hists(pi0_ana_[i], i, ana);
    pi0_analysis_cut(pi0_true, pi0_pm_ana_[i], i%(npi0ana/2), i>=npi0ana/2);
    fill_pi0_analysis_hists(pi0_pm_ana_[i], i, pm_ana);
  }
}

void AnaTda::pi0_analysis_cut(RhoCandList& org, RhoCandList& dest, const int& id, bool mcut = false) {
  for (int i = 0; i < org.GetLength(); ++i) {
    RhoCandidate *_g1 = org[i]->Daughter(0);
    RhoCandidate *_g2 = org[i]->Daughter(1);
    if ( oa(_g1,_g2) > pi0oacut[id]
      and _g1->Energy() > pi0ecut_min[id] and _g1->Energy() < pi0ecut_max[id]
      and _g2->Energy() > pi0ecut_min[id] and _g2->Energy() < pi0ecut_max[id]) {
      if (!mcut || (mcut && org[i]->M() > pi0mcut_min and org[i]->M() < pi0mcut_max) ) {
        dest.Append(org[i]);
      }
    }
  }
}

void AnaTda::fill_pi0_analysis_hists(RhoCandList& list, const int& type) {
  fill_pair_mass(list, h_m_gg[type]);
  fill_pair_oa(list, h_oa_gg[type]);
}

void AnaTda::fill_pi0_analysis_hists(RhoCandList& list, const int& id, const int& type) {
  fill_pair_mass(list, h_m_gg_pi0ana[id][type]);
}

bool AnaTda::primary_match_jpsi(RhoCandidate *c, double &dist) {
  return primary_match_pair(443, c, dist);
}

void AnaTda::jpsi_truth_match() {
  RhoCandidate *_match;
  double dist_min = 1e9;
  jpsi.Combine(ep, em);
  jpsi.SetType(443);
  for (int j = 0; j < jpsi.GetLength(); ++j) {
    double dist = 0;
    if ( primary_match_jpsi(jpsi[j], dist) ) {
      if (dist < dist_min) {
        _match = jpsi[j];
        dist_min = dist;
      }
    }
  }
  if (dist_min<1e8) {
    jpsi_true.Append(_match);
    fill_elecs_from_jpsi();
  }
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
  jpsi_analysis_cut(jpsi_true, jpsi_pm_ana);
  fill_jpsi_analysis_hists(jpsi_pm_ana,pm_ana);
}

void AnaTda::jpsi_analysis_cut(RhoCandList& org, RhoCandList& dest) {
  for (int j = 0; j < org.GetLength(); ++j) {
    const double m = org[j]->M();
    //if ( jpsi_mcut_min < m && m < jpsi_mcut_max ) {
      dest.Append(org[j]);
      //}
  }
}

void AnaTda::fill_jpsi_analysis_hists(RhoCandList& org, const int& type) {
  fill_pair_mass(org, h_m_epem[type]);
}

void AnaTda::Finish() {

  h_m_pippim->Write();
  h_m_pippim_tr->Write();

  h_num_g->Write();
  h_num_epm->Write();
  h_num_pipm->Write();
  h_mom_the_pipm->Write();
  h_num_g_tr->Write();
  h_num_epm_tr->Write();
  h_num_pipm_tr->Write();
  h_mom_the_pipm_tr->Write();

  for (int it=4; it<nhist; ++it) write_manual_kin_fit_hists(it);

  //for (int it=0; it<nhist; ++it) {
  for (int it=0; it<nhist; ++it) h_oa_gg[it]->Write();
  for (int it=0; it<nhist; ++it) h_m_gg[it]->Write();
  for (int it=0; it<nhist; ++it) {
    for (int ia = 0; ia< npi0ana; ++ia) {
      if (it < 4) continue;
      h_m_gg_pi0ana[ia][it]->Write();
    }
  }
  for (int it=0; it<nhist; ++it) h_e_g[it]->Write();

  for (int it=0; it<nhist; ++it) h_m_epem[it]->Write();
  for (int it=0; it<nhist; ++it) h_mom_epm[it]->Write();
  for (int it=0; it<nhist; ++it) h_mom_the_epm[it]->Write();

}




void AnaTda::def_manual_kin_fit_hists(const int &type ) {
  const char *dth = "#theta_{#gamma#gamma}+#theta_{e^{+}e^{-}}";
  const char *mtot = "M^{inv}_{#gamma#gamma - e^{+}e^{-}}";
  const char *mgg = "M^{inv}_{#gamma#gamma}";
  const char *mepem = "M^{inv}_{e^{+}e^{-}}";
  const char *n[nhist] = {"all_rec", "true_rec", "truepi0jpsi_rec", "true_mc", "ana", "pm_ana"};
  const char *t[nhist] = {"all (Reco)", "truth match (Reco)", "from true J/#psi-#pi^{0} (Reco)",
       "from true J/#psi-#pi^{0} (MC)", "after analysis cuts", "from true J/#psi-#pi^{0} after ana. cut"};

  //h_dth_gg_epair[type] = new TH1F(Form("h_dth_gg_epair_%s",n[type]),
  //  Form("%s (%s);#Delta#theta[rad]",dth,t[type]),
  //  100, 0, TMath::Pi());

  //h_mass_gg_epair[type] = new TH1F(Form("h_mass_gg_epair_%s",n[type]),
  //  Form("%s (%s);#Delta#theta[rad]",dth,t[type]),
  //  100, 0, 5.);

  h_dth_vs_mass_gg_epair[type] = new
    TH2F(Form("h_dth_vs_mass_gg_epair_%s",n[type]),
	 Form("%s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta[rad]",dth,mtot,t[type], mtot),
	 200, 2.9, 3.7, 200, -TMath::Pi()/2., TMath::Pi()/2.);

  h_dth_gg_epair_vs_mass_gg[type] = new
    TH2F(Form("h_dth_gg_epair_vs_mass_gg_%s",n[type]),
	 Form("%svs. %s (%s);%s[GeV/c^{2}];#Delta#theta[rad]",dth, mgg,t[type], mgg),
	 200, 0, 0.2, 200, -TMath::Pi()/2., TMath::Pi()/2.);

  h_mass_gg_epair_vs_mass_gg[type] = new
    TH2F(Form("h_mass_gg_epair_vs_mass_gg_%s",n[type]),
	 Form("%s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta[rad]",dth,mgg,t[type],mgg),
	 200, 0, 0.2, 200, 2.9, 3.7);

  // cts = closest-to-s, btb = most-back-to-back
  h_dth_vs_mass_gg_epair_btb[type] = new
    TH2F(Form("h_dth_vs_mass_gg_epair_btb_%s",n[type]),
	 Form("most back-to-back %s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta", dth, mtot, t[type], mtot),
	 200, 2.9, 3.7, 200, -TMath::Pi()/2., TMath::Pi()/2.);

  h_dth_vs_mass_gg_epair_cts[type] = new
    TH2F(Form("h_dth_vs_mass_gg_epair_cts_%s",n[type]),
	 Form("closest to #sqrt{s} %s vs. %s (%s);%s[GeV/c^{2}];#Delta#theta", dth, mtot, t[type], mtot),
	 200, 2.9, 3.7, 200, -TMath::Pi()/2., TMath::Pi()/2.);

  h_m_gg_btb[type] = new
    TH1F(Form("h_m_gg_btb_%s",n[type]),
	 Form("most back-to-back %s (%s);%s[GeV/c^{2}]", mgg, t[type], mgg),
	 100, 0, 0.2);

  h_m_gg_cts[type] = new
    TH1F(Form("h_m_gg_cts_%s",n[type]),
	 Form("closest to #sqrt{s} %s (%s); %s[GeV/c^{2}]", mgg, t[type], mgg),
	 100, 0, 0.2);

  h_m_epem_btb[type] = new
    TH1F(Form("h_m_epem_btb_%s",n[type]),
	 Form("most back-to-back %s (%s);%s[GeV/c^{2}]", mepem, t[type], mepem),
	 100, 0, 0.2);

  h_m_epem_cts[type] = new
    TH1F(Form("h_m_epem_cts_%s",n[type]),
	 Form("closest to #sqrt{s} %s (%s); %s[GeV/c^{2}]", mepem, t[type], mepem),
	 100, 0, 0.2);

}

void AnaTda::calc_kin(
  RhoCandidate* _gg, RhoCandidate *_epem,
  double &m, double &mgg, double &dth) {
  TLorentzVector p4gg = _gg->P4();
  TLorentzVector p4epair = _epem->P4();
  m = (p4gg+p4epair).M();
  mgg = (p4gg).M();
  p4gg.Boost(boost_to_cm);
  p4epair.Boost(boost_to_cm);
  dth = fabs(p4gg.Vect().Theta() + p4epair.Vect().Theta());
}

void AnaTda::pi0_kinematic_selection(RhoCandList& org_gg, RhoCandList& org_epem, const int& type) {
  // pi0 selection based on kinematics (opposite to lepton/pion
  // pair in CM, blieveable invariant mass - close to sqrt(s))
  // start with the list of all gg pairs and plot cm-OA wrt
  // lepton pair, and inv-mass of full system

  // cts = closest-to-s, btb = most-back-to-back
  int i_btb_e = -1, i_cts_e = -1;
  int i_btb_g = -1, i_cts_g = -1;
  double diff_2pi_min = 1e9, diff_s_min = 1e9;
  double m, mgg, dth;

  for (int i = 0; i < org_epem.GetLength(); ++i) {
    for (int j = 0; j < org_gg.GetLength(); ++j) {
      calc_kin(org_gg[j], org_epem[i], m, mgg, dth);
      const double diff_2pi = TMath::Pi() - dth;
      const double diff_s = fabs(sqrt_s - m);
      if (diff_2pi < diff_2pi_min) {
        i_btb_e = i;
        i_btb_g = j;
        diff_2pi_min = diff_2pi;
      }
      if (diff_s < diff_s_min) {
        i_cts_e = i;
        i_cts_g = j;
        diff_s_min = diff_s;
      }
      h_dth_gg_epair_vs_mass_gg[type]->Fill(mgg, dth);
      h_mass_gg_epair_vs_mass_gg[type]->Fill(mgg, m);
      h_dth_vs_mass_gg_epair[type]->Fill(m, dth);
    }
  }

  // now that we have closest pair in OA and M to what we want, fill
  if (i_btb_e >= 0 and i_btb_g >= 0 ) {
    calc_kin(org_gg[i_btb_g], org_epem[i_btb_e], m, mgg, dth);
    h_dth_vs_mass_gg_epair_btb[type]->Fill(m, dth);
    h_m_gg_btb[type]->Fill(org_gg[i_btb_g]->M());
    h_m_epem_btb[type]->Fill(org_gg[i_btb_g]->M());
    //pi0_btb.Append(gg[i_btb]);
  }

  if (i_cts_e >= 0 and i_cts_g >= 0) {
    calc_kin(gg[i_cts_g], org_epem[i_cts_e], m, mgg, dth);
    h_dth_vs_mass_gg_epair_cts[type]->Fill(m, dth);
    h_m_gg_cts[type]->Fill(gg[i_cts_g]->M());
    h_m_epem_cts[type]->Fill(gg[i_cts_g]->M());
      //pi0_cts.Append(gg[i_cts]);
  }

}

void AnaTda::write_manual_kin_fit_hists(const int& type){
  //h_dth_gg_epair[type]->Write();
  //h_mass_gg_epair[type]->Write();
  h_dth_vs_mass_gg_epair[type]->Write();
  h_dth_gg_epair_vs_mass_gg[type]->Write();
  h_mass_gg_epair_vs_mass_gg[type]->Write();
  h_dth_vs_mass_gg_epair_btb[type]->Write();      // btb = most-back-to-back
  h_dth_vs_mass_gg_epair_cts[type]->Write();      // cts = closest-to-s
  h_m_gg_btb[type]->Write();     // btb = most-back-to-back
  h_m_gg_cts[type]->Write();     // cts = closest-to-s
  h_m_epem_btb[type]->Write();     // btb = most-back-to-back
  h_m_epem_cts[type]->Write();     // cts = closest-to-s
}







// --- below are some old bad ideas
void AnaTda::def_kin_fit_hists(const int& type, const int& hyp) {
  const char *n[] = {"","_btb","_cts"};
  const char *nhyp[] = {"epem", "pippim"};
  const char *thyp[] = {"e^{+}e^{-}", "#pi^{+}#pi^{-}"};

  h_4c_chi2[type][hyp] = new TH1F(
    Form("h_4c_chi2%s_%spi0", n[type], nhyp[hyp]),
    Form("#chi^{2} of 4C fit on %s#pi^{0} system; #chi^{2}", thyp[hyp]),
    100, 0, 200);

  h_4c_prob[type][hyp] = new TH1F(
    Form("h_4c_prob%s_%spi0", n[type], nhyp[hyp]),
    Form("Probability of #chi^{2} of 4C fit on %s#pi^{0} system; prob", thyp[hyp]),
    100, -0.1, 1.1);

  h_4c_m[type][hyp] = new TH1F(
    Form("h_4c_m%s_%s_4c", n[type], nhyp[hyp]),
    Form("Invariant mass of %s After 4C Fit;M_{%s}[GeV/c^{2}]", thyp[hyp], thyp[hyp]),
    100, 0, 4.5);

  h_4c_prob_vs_m[type][hyp] = new TH2F(
    Form("h_4c_prob_vs_m%s_%s_4c", n[type], nhyp[hyp]),
    Form("Prob. vs. Invariant mass of %s After 4C Fit;M_{%s}[GeV/c^{2}]; prob", thyp[hyp], thyp[hyp]),
    100, 2.9, 3.7, 100, 0, 1);
}

void AnaTda::def_kin_fit_hists() {
  for (int it = 0; it < 3; ++it)
    for (int ih = 0; ih < 2; ++ih)
      def_kin_fit_hists(it, ih);
}

void AnaTda::write_kin_fit_hists() {
  for (int it = 0; it < 3; ++it) {
    for (int ih = 0; ih < 2; ++ih) {
      h_4c_chi2[it][ih]->Write();
      h_4c_prob[it][ih]->Write();
      h_4c_m[it][ih]->Write();
      h_4c_prob_vs_m[it][ih]->Write();
    }
  }
}

//void AnaTda::kin_fit_full_sys(RhoCandList& org, TH1F* h_chi2, TH1F* h_prob, TH2F* h_prob_m) {
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

void AnaTda::kin_fit_pi0_nearest(RhoCandList& org, const int& type) {
  for (int j = 0; j < org.GetLength(); ++j) {
    //if (type==0)
    //  h_m_epem_pi0n->Fill(org[j]->M());
    //else
    //  h_m_pippim_pi0n->Fill(org[j]->M());
    PndKinFitter fitter(org[j]);  // instantiate the kin fitter
    fitter.Add4MomConstraint(ini);  //  set 4 constraint
    fitter.Fit();  //  do fit
    double chi2_4c = fitter.GetChi2();  //  get chi2 of fit
    double prob_4c = fitter.GetProb();  //  access probability of fit
    h_4c_chi2[0][type]->Fill(chi2_4c);
    h_4c_prob[0][type]->Fill(prob_4c);
    if ( prob_4c > 0.01 ) {
      RhoCandidate *_fit = org[j]->Daughter(0)->GetFit();
      h_4c_m[0][type]->Fill(_fit->M());
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
//  jpsi.SetType(443);
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
//  pi0.Combine(g1, g2);
//  pi0.SetType(111);
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

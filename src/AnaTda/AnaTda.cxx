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
#include "FairRuntimeDb.h"
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

// other headers
//#include "PndEmcBump.h"
//#include "PndEmcCluster.h"
//#include "PndEmcBump.h"
//#include "PndEmcRecoHit.h"
//#include "PndEmcMapper.h"
#include "PndMCTrack.h"

using std::cout;
using std::endl;
using std::vector;

AnaTda::AnaTda(const int &brem = 0) :
  FairTask("Radiation Length Profiler"),nevt(0) {
  fBremCorr = (brem != 0);
}

AnaTda::~AnaTda() { }

void AnaTda::def_hists() {

  // *** create some histograms
  hjpsim_all = new TH1F("hjpsim_all","J/#psi mass (all)",100,0,4.5);
  hpi0m_all  = new TH1F("hpi0m_all","#pi^{0} mass (all)",100,0,0.2);
  hjpsipi0m_all  = new TH1F("hjpsipi0m_all","J/#psi-#pi^{0} mass (all)",100,0,5);

  hjpsim_ftm = new TH1F("hjpsim_ftm","J/#psi mass (full truth match)",100,0,4.5);
  hpi0m_ftm  = new TH1F("hpi0m_ftm","#pi^{0} mass (full truth match)",100,0,0.2);
  hjpsipi0m_ftm  = new TH1F("hjpsipi0m_ftm","J/#psi-#pi^{0} mass (full truth match)",100,0,5);

  hjpsim_nm = new TH1F("hjpsim_nm","J/#psi mass (no truth match)",100,0,4.5);
  hpi0m_nm  = new TH1F("hpi0m_nm","#pi^{0} mass (no truth match)",100,0,0.2);
  hjpsipi0m_nm  = new TH1F("hjpsipi0m_nm","J/#psi-#pi^{0} mass (no truth match)",100,0,5);

  hjpsim_diff = new TH1F("hjpsim_diff","J/#psi mass diff to truth",100,-2,2);
  hpi0m_diff  = new TH1F("hpi0m_diff","#pi^{0} mass diff to truth",100,-2,2);
  hjpsipi0m_diff  = new TH1F("hjpsipi0m_diff","J/#psi-#pi^{0} mass diff to truth",100,-2,2);

  h_m_epem = new TH1F("h_m_epem", "Invariant mass of e^{+}e^{-};M_{e^{+}e^{-}}[GeV/c^{2}]", 100, 0, 5.);
  h_m_epem_tr = new TH1F("h_m_epem_tr", "Invariant mass of e^{+}e^{-} (Truth Match);M_{e^{+}e^{-}}[GeV/c^{2}]", 100, 0, 5.);
  h_m_pippim = new TH1F("h_m_pippim", "Invariant mass of #pi^{+}#pi^{-};M_{#pi^{+}#pi^{-}}[GeV/c^{2}]", 100, 0, 5.);
  h_m_pippim_tr = new TH1F("h_m_pippim_tr", "Invariant mass of #pi^{+}#pi^{-} (Truth Match);M_{#pi^{+}#pi^{-}}[GeV/c^{2}]", 100, 0, 5.);
  h_m_gg = new TH1F("h_m_gg", "Invariant mass of #gamma-#gamma pairs;M_{#gamma_{1}#gamma_{2}}[GeV/c^{2}]", 100, 0, 0.2);
  h_m_gg_tr = new TH1F("h_m_gg_tr", "Invariant mass of #gamma-#gamma pairs (Truth Match);M_{#gamma_{1}#gamma_{2}}[GeV/c^{2}]", 100, 0, 0.2);

  h_m_pi0n = new TH1F("h_m_pi0n","#pi^{0} mass (nearest to p0 pdg mass)",100,0,0.2);
  h_m_epem_pi0n = new TH1F("h_m_epem_pi0n","e^{+}e^{-}#pi^{0} mass (nearest to p0 pdg mass)",100,0,5.0);
  h_m_pippim_pi0n = new TH1F("h_m_pippim_pi0n","#pi^{+}#pi^{-}#pi^{0} mass (nearest to p0 pdg mass)",100,0,5.0);

  h_num_g = new TH1F("h_num_g", "Number of #gamma per event;N_{#gamma}", 25, 0, 25);
  h_num_epm = new TH1F("h_num_epm", "Number of electron tracks per event;N_{e^{#pm}}", 10, 0, 10);
  h_num_pipm = new TH1F("h_num_pipm", "Number of charged pion tracks per event;N_{#pi^{#pm}}", 10, 0, 10);
  h_e_g = new TH1F("h_e_g", "Energy of #gamma;E[GeV]", 100, 0, 2);
  h_mom_the_epm = new TH2F("h_mom_the_epm", "#theta vs. momentum of electron tracks;p[GeV/c];#theta[rad]", 100, 0, 5, 100, 0, TMath::Pi());
  h_mom_the_pipm = new TH2F("h_mom_the_pipm", "#theta vs. momentum of charge pion tracks;p[GeV/c];#theta[rad]", 100, 0, 5, 100, 0, TMath::Pi());
  h_num_g_tr = new TH1F("h_num_g_tr", "Number of #gamma per event (Truth Match);N_{#gamma}", 25, 0, 25);
  h_num_epm_tr = new TH1F("h_num_epm_tr", "Number of electron tracks per event (Truth Match);N_{e^{#pm}}", 10, 0, 10);
  h_num_pipm_tr = new TH1F("h_num_pipm_tr", "Number of charged pion tracks per event (Truth Match);N_{#pi^{#pm}}", 10, 0, 10);
  h_e_g_tr = new TH1F("h_e_g_tr", "Energy of #gamma (Truth Match);E[GeV]", 100, 0, 2);
  h_mom_the_epm_tr = new TH2F("h_mom_the_epm_tr", "#theta vs. momentum of electron tracks (Truth Match);p[GeV/c];#theta[rad]", 100, 0, 5, 100, 0, TMath::Pi());
  h_mom_the_pipm_tr = new TH2F("h_mom_the_pipm_tr", "#theta vs. momentum of charge pion tracks (Truth Match);p[GeV/c];#theta[rad]", 100, 0, 5, 100, 0, TMath::Pi());

  h_4c_chi2_epempi0 = new TH1F("h_4c_chi2_epempi0","#chi^2 of 4C fit on e^{+}e^{-}#pi^{0} system; #chi^2", 100, 0, 200);
  h_4c_prob_epempi0 = new TH1F("h_4c_prob_epempi0","Probability of #chi^2 of 4C fit on e^{+}e^{-}#pi^{0} system; prob", 100, -0.1, 1.1);
  h_4c_m_epem = new TH1F("h_4c_m_epem_4c", "Invariant mass of e^{+}e^{-} After 4C Fit;M_{e^{+}e^{-}}[GeV/c^{2}]", 100, 0, 5);

  h_4c_chi2_pippimpi0 = new TH1F("h_4c_chi2_pippimpi0","#chi^2 of 4C fit on #pi^{+}#pi^{-}#pi^{0} system; #chi^2", 100, 0, 200);
  h_4c_prob_pippimpi0 = new TH1F("h_4c_prob_pippimpi0","Probability of #chi^2 of 4C fit on #pi^{+}#pi^{-}#pi^{0} system; prob", 100, -0.1, 1.1);
  h_4c_m_pippim = new TH1F("h_4c_m_pippim_4c", "Invariant mass of #pi^{+}#pi^{-} After 4C Fit;M_{#pi^{+}#pi^{-}}[GeV/c^{2}]", 100, 0, 5);



}

InitStatus AnaTda::Init()
{

  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- AnaTda::Init: " << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  def_hists();

  fAnalysis = new PndAnalysis();

  // *** Mass selector for the jpsi cands
  m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi
  m0_pi0 = TDatabasePDG::Instance()->GetParticle("pi0")->Mass();   // Get nominal PDG mass of the J/psi
  cout << "AnaTda::Init() m0_jpsi= " << m0_jpsi << " m0_pi0= " << m0_pi0 << endl;
  jpsiMassSel = new RhoMassParticleSelector("jpsi",m0_jpsi,1.0);
  pi0MassSel = new RhoMassParticleSelector("pi0",m0_pi0,0.06);

  // *** the lorentz vector of the initial jpsi-pi0 system
  ini = TLorentzVector(0, 0, 5.513, 6.53023); //(0, 0, 6.231552, 7.240065);

  return kSUCCESS;

}

int AnaTda::SelectTruePid(PndAnalysis *ana, RhoCandList &l)
{
  int removed = 0;
  for (int ii=l.GetLength()-1;ii>=0;--ii)
    {
      if ( !(ana->McTruthMatch(l[ii])) )
	{
	  l.Remove(l[ii]);
	  removed++;
	}
    }
  return removed;
}

void AnaTda::Exec(Option_t* opt)
{
  if (++nevt%100==0)
    cout << "===== AnaTda::Exec -- Event " << nevt << " ====="<< endl;

  fAnalysis->GetEvent();

  // *** Select with no PID info ('All'); type and mass are set
  RhoCandList pip, pim, ep, em, g1, g2;
  fAnalysis->FillList(ep, (fBremCorr?"BremElectronAllPlus":"ElectronAllPlus"));
  fAnalysis->FillList(em, (fBremCorr?"BremElectronAllMinus":"ElectronAllMinus"));
  fAnalysis->FillList(g1, "Neutral");
  fAnalysis->FillList(g2, "Neutral");
  fAnalysis->FillList(pip, (fBremCorr?"BremPionAllPlus":"PionAllPlus"));
  fAnalysis->FillList(pim, (fBremCorr?"BremPionAllMinus":"PionAllMinus"));

  // Single distributions
  h_num_g->Fill(g1.GetLength());
  h_num_epm->Fill(ep.GetLength()+em.GetLength());
  h_num_pipm->Fill(pip.GetLength()+pim.GetLength());
  for (int i=0; i<g1.GetLength(); ++i) h_e_g->Fill(g1[i]->Energy());
  for (int i=0; i<ep.GetLength(); ++i) h_mom_the_epm->Fill(ep[i]->P3().Mag(), ep[i]->P3().Theta() );
  for (int i=0; i<em.GetLength(); ++i) h_mom_the_epm->Fill(em[i]->P3().Mag(), em[i]->P3().Theta() );
  for (int i=0; i<pip.GetLength(); ++i) h_mom_the_pipm->Fill(pip[i]->P3().Mag(), pip[i]->P3().Theta() );
  for (int i=0; i<pim.GetLength(); ++i) h_mom_the_pipm->Fill(pim[i]->P3().Mag(), pim[i]->P3().Theta() );

  RhoCandList pip_tr, pim_tr, ep_tr, em_tr, g1_tr, g2_tr;
  ep.SetType(-11);
  em.SetType(11);
  pip.SetType(211);
  pim.SetType(-211);
  g1.SetType(22);
  g2.SetType(22);
  for (int i=0; i<ep.GetLength(); ++i) if (fAnalysis->McTruthMatch(ep[i])) ep_tr.Append(ep[i]);
  for (int i=0; i<em.GetLength(); ++i) if (fAnalysis->McTruthMatch(em[i])) em_tr.Append(em[i]);
  for (int i=0; i<pip.GetLength(); ++i) if (fAnalysis->McTruthMatch(pip[i])) pip_tr.Append(pip[i]);
  for (int i=0; i<pim.GetLength(); ++i) if (fAnalysis->McTruthMatch(pim[i])) pim_tr.Append(pim[i]);
  for (int i=0; i<g1.GetLength(); ++i) if (fAnalysis->McTruthMatch(g1[i])) g1_tr.Append(g1[i]);
  for (int i=0; i<g2.GetLength(); ++i) if (fAnalysis->McTruthMatch(g2[i])) g2_tr.Append(g2[i]);

  // Single distributions (After truth match)
  h_num_g_tr->Fill(g1_tr.GetLength());
  h_num_epm_tr->Fill(ep_tr.GetLength()+em_tr.GetLength());
  h_num_pipm_tr->Fill(pip_tr.GetLength()+pim_tr.GetLength());
  for (int i=0; i<g1_tr.GetLength(); ++i) h_e_g_tr->Fill(g1_tr[i]->Energy());
  for (int i=0; i<ep_tr.GetLength(); ++i) h_mom_the_epm_tr->Fill(ep_tr[i]->P3().Mag(), ep_tr[i]->P3().Theta() );
  for (int i=0; i<em_tr.GetLength(); ++i) h_mom_the_epm_tr->Fill(em_tr[i]->P3().Mag(), em_tr[i]->P3().Theta() );
  for (int i=0; i<pip_tr.GetLength(); ++i) h_mom_the_pipm_tr->Fill(pip_tr[i]->P3().Mag(), pip_tr[i]->P3().Theta() );
  for (int i=0; i<pim_tr.GetLength(); ++i) h_mom_the_pipm_tr->Fill(pim_tr[i]->P3().Mag(), pim_tr[i]->P3().Theta() );

  RhoCandList epem, pippim, gg;
  RhoCandList epem_tr, pippim_tr, gg_tr;
  epem.Combine(ep, em);
  epem.SetType(443);
  epem_tr.Combine(ep_tr, em_tr);
  pippim.Combine(pip, pim);
  pippim_tr.Combine(pip_tr, pim_tr);
  gg.Combine(g1, g2);
  gg_tr.Combine(g1_tr, g2_tr);

  for (int j=0; j<epem.GetLength(); ++j) h_m_epem->Fill(epem[j]->M());
  for (int j=0; j<epem_tr.GetLength(); ++j) h_m_epem_tr->Fill(epem_tr[j]->M());
  for (int j=0; j<pippim.GetLength(); ++j) h_m_pippim->Fill(pippim[j]->M());
  for (int j=0; j<pippim_tr.GetLength(); ++j) h_m_pippim_tr->Fill(pippim_tr[j]->M());
  for (int j=0; j<gg.GetLength(); ++j) h_m_gg->Fill(gg[j]->M());
  for (int j=0; j<gg_tr.GetLength(); ++j) h_m_gg_tr->Fill(gg_tr[j]->M());

  RhoCandList pi0;
  pi0.Combine(g1,g2);
  pi0.SetType(111);

  RhoCandList pi0_true;
  for (int j=0; j<pi0.GetLength(); ++j) {
    hpi0m_all->Fill( pi0[j]->M() );
    if ( fAnalysis->McTruthMatch(pi0[j]) ) {
  	hpi0m_ftm->Fill( pi0[j]->M() );
  	hpi0m_diff->Fill( pi0[j]->GetMcTruth()->M() - pi0[j]->M() );
  	pi0_true.Append(pi0[j]);
    } else {
  	hpi0m_nm->Fill( pi0[j]->M() );
    }
  }

  RhoCandList pi0nearest;
  double dm_min = 1e10;
  int min_j = -1;
  for (int j=0; j<pi0.GetLength(); ++j) {
    const double dm = fabs(pi0[j]->M() - m0_pi0);
    if ( dm < dm_min ) {
  	dm_min = dm;
  	min_j = j;
    }
  }
  if (min_j>=0) {
    pi0nearest.Append(pi0[min_j]);
    h_m_pi0n->Fill(pi0nearest[0]->M());
  }



  // following is relevant only for signal
  RhoCandList jpsi;
  jpsi.Combine(ep, em);
  jpsi.SetType(443);

  RhoCandList jpsi_true;
  for (int j=0;j<jpsi.GetLength();++j) {
    hjpsim_all->Fill( jpsi[j]->M() );
    if (fAnalysis->McTruthMatch(jpsi[j])) {
  	hjpsim_ftm->Fill( jpsi[j]->M() );
  	hjpsim_diff->Fill( jpsi[j]->GetMcTruth()->M() - jpsi[j]->M() );
  	jpsi_true.Append(jpsi[j]);
    } else {
  	hjpsim_nm->Fill( jpsi[j]->M() );
    }
  }



  RhoCandList epem_pi0nearest;
  epem_pi0nearest.Combine(epem, pi0nearest);
  for (int j=0; j<epem_pi0nearest.GetLength(); ++j) {
    h_m_epem_pi0n->Fill( epem_pi0nearest[j]->M() );
    PndKinFitter fitter(epem_pi0nearest[j]);	// instantiate the kin fitter in psi(2S)
    fitter.Add4MomConstraint(ini);	// set 4 constraint
    fitter.Fit();		            // do fit
    double chi2_4c = fitter.GetChi2();	// get chi2 of fit
    double prob_4c = fitter.GetProb();	// access probability of fit
    //cout << "ini.M()= " << ini.M() << " epem_pi0nearest[j].M()= " << epem_pi0nearest[j]->M() <<  "chi2 = " << chi2_4c << " prob= " << prob_4c << endl;
    h_4c_chi2_epempi0->Fill(chi2_4c);
    h_4c_prob_epempi0->Fill(prob_4c);
    if ( prob_4c > 0.01 ) {			// when good enough, fill some histo
      RhoCandidate *epem_fit = epem_pi0nearest[j]->Daughter(0)->GetFit();	// get fitted epem
      h_4c_m_epem->Fill(epem_fit->M());
    }
  }

  RhoCandList pippim_pi0nearest;
  pippim_pi0nearest.Combine(pippim, pi0nearest);
  for (int j=0; j<pippim_pi0nearest.GetLength(); ++j) {
    h_m_pippim_pi0n->Fill( pippim_pi0nearest[j]->M() );
    PndKinFitter fitter(pippim_pi0nearest[j]);	// instantiate the kin fitter in psi(2S)
    fitter.Add4MomConstraint(ini);	// set 4 constraint
    fitter.Fit();		            // do fit
    double chi2_4c = fitter.GetChi2();	// get chi2 of fit
    double prob_4c = fitter.GetProb();	// access probability of fit
    h_4c_chi2_pippimpi0->Fill(chi2_4c);
    h_4c_prob_pippimpi0->Fill(prob_4c);
    if ( prob_4c > 0.01 ) {			// when good enough, fill some histo
      RhoCandidate *pippim_fit = pippim_pi0nearest[j]->Daughter(0)->GetFit();	// get fitted epem
      h_4c_m_pippim->Fill(pippim_fit->M());
    }
  }


}

void AnaTda::Finish()
{

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

  h_m_epem->Write();
  h_m_epem_tr->Write();
  h_m_pippim->Write();
  h_m_pippim_tr->Write();
  h_m_gg->Write();
  h_m_gg_tr->Write();

  h_m_pi0n->Write();
  h_m_epem_pi0n->Write();
  h_m_pippim_pi0n->Write();

  h_4c_chi2_epempi0->Write();
  h_4c_prob_epempi0->Write();
  h_4c_m_epem->Write();
  h_4c_chi2_pippimpi0->Write();
  h_4c_prob_pippimpi0->Write();
  h_4c_m_pippim->Write();

  h_num_g->Write();
  h_num_epm->Write();
  h_num_pipm->Write();
  h_e_g->Write();
  h_mom_the_epm->Write();
  h_mom_the_pipm->Write();
  h_num_g_tr->Write();
  h_num_epm_tr->Write();
  h_num_pipm_tr->Write();
  h_e_g_tr->Write();
  h_mom_the_epm_tr->Write();
  h_mom_the_pipm_tr->Write();

}

ClassImp(AnaTda)

#include "EffHists.h"

#include "FairTask.h"

#include "PndPidProbability.h"
#include "PndPidCandidate.h"
#include "PndMCTrack.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

static const double dpx_min[20] = { -0.0235856, -0.016433, -0.022361, -0.0311639, -0.0408277, -0.049248, -0.0626825, -0.0709645, -0.0834883, -0.0926895, -0.0961089, -0.100969, -0.0978446, -0.111202, -0.121147, -0.100926, -0.116786, -0.114434, -0.12469, -0.126147};
static const double dpx_max[20] = { 0.0238193, 0.0164182, 0.0217275, 0.0310145, 0.0407391, 0.049224, 0.0617949, 0.0714845, 0.0838231, 0.090658, 0.0957098, 0.099138, 0.096258, 0.106792, 0.114417, 0.100926, 0.116978, 0.111493, 0.124899, 0.125351};
static const double dpy_min[20] = { -0.0262052, -0.0177316, -0.0204651, -0.0286626, -0.0384467, -0.0491804, -0.0617386, -0.0560984, -0.057134, -0.0802229, -0.092882, -0.0883073, -0.0934353, -0.0983313, -0.0955621, -0.096944, -0.108257, -0.108708, -0.11006, -0.11348};
static const double dpy_max[20] = { 0.0264953, 0.0180457, 0.0207731, 0.0278699, 0.0384214, 0.04838, 0.0624423, 0.0558169, 0.0569003, 0.0813983, 0.090161, 0.0862677, 0.0935711, 0.0943912, 0.0964051, 0.0971344, 0.109794, 0.107839, 0.110128, 0.116133};
static const double dpz_min[20] = { -0.0161801, -0.0176935, -0.0254173, -0.0366754, -0.052316, -0.0633923, -0.0802275, -0.0916246, -0.103904, -0.12348, -0.118053, -0.121892, -0.146836, -0.144744, -0.144244, -0.154281, -0.168258, -0.188327, -0.16868, -0.16266};
static const double dpz_max[20] = { 0.0141423, 0.0172002, 0.0260039, 0.0375642, 0.052133, 0.0642037, 0.0803328, 0.0917838, 0.106182, 0.121883, 0.115329, 0.117401, 0.144266, 0.145639, 0.140742, 0.152066, 0.17101, 0.176431, 0.164111, 0.161177};

const TString EffHists::s_spc[EffHists::nsp_max] = {"posit","muplus","piplus","kplus","prot","elec","muminus","piminus","kminus","antiprot"};
const TString EffHists::s_spc_tex[EffHists::nsp_max] = {"e^{+}","#mu^{+}","#pi^{+}","K^{+}","p","e^{-}","#mu^{-}","#pi^{-}","K^{-}","#bar{p}"};
const TString EffHists::s_pid[EffHists::npid_max] = {"e_id", "mu_id", "pi_id", "k_id", "prot_id"};
const TString EffHists::s_det[EffHists::ndet] = {"emc", "stt", "mvd", "dirc", "disc"};

EffHists::EffHists(int a_sp):
  m_sp(a_sp),
  verb(false),
  det_var_max{1.5, 20, 10, 90, 90},
  mom_max(10.0),
  the_max(180.0),
  fAna()
{

  for (int iprob_cut=0; iprob_cut < 9; ++iprob_cut) {
    prob_cut.push_back( 0.1*iprob_cut );
    //prob_cut[iprob_cut] = (1.0/double(nprob_cut-1))*double(iprob_cut);
    //cout << "prob_cut[" << iprob_cut << "] = " << prob_cut[iprob_cut] << endl;
  }
  for (int iprob_cut=0; iprob_cut < 9; ++iprob_cut) {
    prob_cut.push_back( 0.9 + (0.01*iprob_cut) );
  }

  for (int iprob_cut=0; iprob_cut < 10; ++iprob_cut) {
    prob_cut.push_back( 0.99 + (0.001*iprob_cut) );
  }

}

EffHists::~EffHists() {
}

TClonesArray* EffHists::init_tca(TString name) {
  TClonesArray *tca = dynamic_cast<TClonesArray *> (m_ioman->GetObject(name));
  if ( ! tca ) {
    cout << "-W- EffHists::Init: "  << "No " << name << " array!" << endl;
    return NULL;
  } else {
    cout << "-I- EffHists::Init: "  << " finished reading " << name << " array!" << endl;
    return tca;
  }
}

InitStatus EffHists::init_tcas() {
  m_ioman = FairRootManager::Instance();
  if ( ! m_ioman ){
    cout << "-E- EffHists::Init: "
	 << "RootManager not instantiated!" << endl;
    return kFATAL;
  }
  m_cand_array = init_tca( "PidChargedCand");
  m_mc_array = init_tca( "MCTrack");
  m_drc_array = init_tca( "PidAlgoDrc");
  m_disc_array = init_tca( "PidAlgoDisc");
  m_stt_array = init_tca( "PidAlgoStt");
  m_mvd_array = init_tca( "PidAlgoMvd");
  m_emcb_array = init_tca( "PidAlgoEmcBayes");
  return kSUCCESS;
}

void EffHists::init_hists() {
  int nbin = 200;
  h_dpx = new TH2F("h_dpx", "px_{MC}-px_{REC} vs p_{MC};px_{MC}-px_{REC};p_{MC}",2000,-1,1,200,0,mom_max);
  h_dpy = new TH2F("h_dpy", "py_{MC}-py_{REC} vs p_{MC};py_{MC}-py_{REC};p_{MC}",2000,-1,1,200,0,mom_max);
  h_dpz = new TH2F("h_dpz", "pz_{MC}-pz_{REC} vs p_{MC};pz_{MC}-pz_{REC};p_{MC}",2000,-1,1,200,0,mom_max);

  for (int ipid=0; ipid < npid_max; ++ipid) {
    const char *tt = s_pid[ipid];
    h_prob[ipid][0] = new TH2F(Form("h_prob_comb_%s_vs_mom",tt),Form("p_{COMB}(%s) vs. p_{REC};p_{REC}[GeV/c];p_{COMB}(%s)",tt,tt), 200, 0, mom_max, 200, 0, 1);
    h_prob[ipid][1] = new TH2F(Form("h_prob_comb_sd_%s_vs_mom",tt),Form("p_{COMB}(%s, Sans dE/dx) vs. p_{REC};p_{REC}[GeV/c];p^{SD}_{COMB}(%s)",tt,tt), 200, 0, mom_max, 200, 0, 1);
    h_prob[ipid][2] = new TH2F(Form("h_prob_emc_%s_vs_eop",tt),Form("p_{EMC}(%s) vs. E/p; E/p; p_{EMC}(%s)",tt,tt), 200, 0, det_var_max[0], 200, 0, 1);
    h_prob[ipid][3] = new TH2F(Form("h_prob_stt_%s_vs_dedx",tt),Form("p_{STT}(%s) vs. dE/dx(STT);dE/dx(STT); p_{STT}(%s)",tt,tt), 200, 0, det_var_max[1], 200, 0, 1);
    h_prob[ipid][4] = new TH2F(Form("h_prob_mvd_%s_vs_dedx",tt),Form("p_{MVD}(%s) vs. dE/dx(MVD);dE/dx(MVD); p_{MVD}(%s)",tt,tt), 200, 0, det_var_max[2], 200, 0, 1);
    h_prob[ipid][5] = new TH2F(Form("h_prob_drc_%s_vs_thec",tt),Form("p_{DRC}(%s) vs. #theta_{C}(DRC);#theta_{C}(DRC); p_{DRC}(%s)",tt,tt), 200, 0, det_var_max[3], 200, 0, 1);
    h_prob[ipid][6] = new TH2F(Form("h_prob_disc_%s_vs_thec",tt),Form("p_{DISC}(%s) vs. #theta_{C}(DISC);#theta_{C}(DISC); p_{DISC}(%s)",tt,tt), 200, 0, det_var_max[4], 200, 0, 1);
  }

  for (int ipid = 0; ipid < npid_max; ++ipid) {

    TString title = Form("%s",s_spc_tex[m_sp].Data());
    for (int ipc=0; ipc < prob_cut.size(); ++ipc) {
      if (ipc>=nprob_cut) break;

      title = Form("efficiency of %s to pass %s cuts at prob>%4.2f%%",s_spc_tex[m_sp].Data(),s_pid[ipid].Data(),prob_cut[ipc]*100.);
      eff1d_mom[ipc][ipid] = new TEfficiency(Form("eff1d_mom_%s",s_pid[ipid].Data()), title+" (MC);p_{MC}[GeV/c]", nbin, 0, mom_max);
      eff1d_the[ipc][ipid] = new TEfficiency(Form("eff1d_the_%s",s_pid[ipid].Data()), title+" (MC);#theta_{MC}[rad]", nbin, 0, the_max);
      eff2d[ipc][ipid] = new TEfficiency(Form("eff2d_%s",s_pid[ipid].Data()), title+" (MC);p_{MC}[GeV/c];#theta_{MC}[rad]", nbin, 0, mom_max, nbin, 0, the_max);

      eff1d_mom_indiv[ipc][ipid] = new TEfficiency(Form("eff1d_mom_%s_indiv",s_pid[ipid].Data()), title+" (TOP);p_{MC}[GeV/c]", nbin, 0, mom_max);
      eff1d_the_indiv[ipc][ipid] = new TEfficiency(Form("eff1d_the_%s_indiv",s_pid[ipid].Data()), title+" (TOP);#theta_{MC}[rad]", nbin, 0, the_max);
      eff2d_indiv[ipc][ipid] = new TEfficiency(Form("eff2d_%s_indiv",s_pid[ipid].Data()), title+" (TOP);p_{MC}[GeV/c];#theta_{MC}[rad]", nbin, 0, mom_max, nbin, 0, the_max);

      eff1d_mom_sd[ipc][ipid] = new TEfficiency(Form("eff1d_mom_%s_sd",s_pid[ipid].Data()), title+" Sans dE/dx;p_{MC}[GeV/c]", nbin, 0, mom_max);
      eff1d_the_sd[ipc][ipid] = new TEfficiency(Form("eff1d_the_%s_sd",s_pid[ipid].Data()), title+" Sans dE/dx;#theta_{MC}[rad]", nbin, 0, the_max);
      eff2d_sd[ipc][ipid] = new TEfficiency(Form("eff2d_%s_sd",s_pid[ipid].Data()), title+" Sans dE/dx;p_{MC}[GeV/c];#theta_{MC}[rad]", nbin, 0, mom_max, nbin, 0, the_max);

      eff1d_mom_sd_indiv[ipc][ipid] = new TEfficiency(Form("eff1d_mom_%s_sd_indiv",s_pid[ipid].Data()), title+" Sans dE/dx (TOP);p_{MC}[GeV/c]", nbin, 0, mom_max);
      eff1d_the_sd_indiv[ipc][ipid] = new TEfficiency(Form("eff1d_the_%s_sd_indiv",s_pid[ipid].Data()), title+" Sans dE/dx (TOP);#theta_{MC}[rad]", nbin, 0, the_max);
      eff2d_sd_indiv[ipc][ipid] = new TEfficiency(Form("eff2d_%s_sd_indiv",s_pid[ipid].Data()), title+" Sans dE/dx  (TOP);p_{MC}[GeV/c];#theta_{MC}[rad]", nbin, 0, mom_max, nbin, 0, the_max);

    }

  }

  h_emc_mom_mc = new TH2F("h_emc_mom_mc", "EMC: E/p vs mom_{MC}; p_{MC}[GeV/c]; E/p", 200, 0, mom_max, 200, 0, det_var_max[0]);
  h_stt_mom_mc = new TH2F("h_stt_mom_mc", "STT: dedx vs mom_{MC}; p_{MC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[1]);
  h_mvd_mom_mc = new TH2F("h_mvd_mom_mc", "MVD: dedx vs mom_{MC}; p_{MC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[2]);
  h_dirc_mom_mc = new TH2F("h_dirc_mom_mc", "DIRC: #tehta_{C} vs mom_{MC}; p_{MC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[3]);
  h_disc_mom_mc = new TH2F("h_disc_mom_mc", "DISC: #tehta_{C} vs mom_{MC}; p_{MC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[4]);

  h_emc_mom_rec = new TH2F("h_emc_mom_rec", "EMC: E/p vs mom_{REC}; p_{REC}[GeV/c]; E/p", 200, 0, mom_max, 200, 0, det_var_max[0]);
  h_stt_mom_rec = new TH2F("h_stt_mom_rec", "STT: dedx vs mom_{REC}; p_{REC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[1]);
  h_mvd_mom_rec = new TH2F("h_mvd_mom_rec", "MVD: dedx vs mom_{REC}; p_{REC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[2]);
  h_dirc_mom_rec = new TH2F("h_dirc_mom_rec", "DIRC: #tehta_{C} vs mom_{REC}; p_{REC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[3]);
  h_disc_mom_rec = new TH2F("h_disc_mom_rec", "DISC: #tehta_{C} vs mom_{REC}; p_{REC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[4]);

  h_emc_th_mc = new TH2F("h_emc_th_mc", "EMC: E/p vs #theta_{MC}; #theta_{MC}[GeV/c]; E/p", 200, 0, the_max, 200, 0, det_var_max[0]);
  h_stt_th_mc = new TH2F("h_stt_th_mc", "STT: dedx vs #theta_{MC}; #theta_{MC}[GeV/c]; dE/dx", 200, 0, the_max, 200, 0, det_var_max[1]);
  h_mvd_th_mc = new TH2F("h_mvd_th_mc", "MVD: dedx vs #theta_{MC}; #theta_{MC}[GeV/c]; dE/dx", 200, 0, the_max, 200, 0, det_var_max[2]);
  h_dirc_th_mc = new TH2F("h_dirc_th_mc", "DIRC: #tehta_{C} vs #theta_{MC}; #theta_{MC}[GeV/c]; #theta_{C}(deg)", 200, 0, the_max, 200, 0, det_var_max[3]);
  h_disc_th_mc = new TH2F("h_disc_th_mc", "DISC: #tehta_{C} vs #theta_{MC}; #theta_{MC}[GeV/c]; #theta_{C}(deg)", 200, 0, the_max, 200, 0, det_var_max[4]);

  h_emc_th_rec = new TH2F("h_emc_th_rec", "EMC: E/p vs #theta_{REC}; #theta_{REC}[GeV/c]; E/p", 200, 0, the_max, 200, 0, det_var_max[0]);
  h_stt_th_rec = new TH2F("h_stt_th_rec", "STT: dedx vs #theta_{REC}; #theta_{REC}[GeV/c]; dE/dx", 200, 0, the_max, 200, 0, det_var_max[1]);
  h_mvd_th_rec = new TH2F("h_mvd_th_rec", "MVD: dedx vs #theta_{REC}; #theta_{REC}[GeV/c]; dE/dx", 200, 0, the_max, 200, 0, det_var_max[2]);
  h_dirc_th_rec = new TH2F("h_dirc_th_rec", "DIRC: #tehta_{C} vs #theta_{REC}; #theta_{REC}[GeV/c]; #theta_{C}(deg)", 200, 0, the_max, 200, 0, det_var_max[3]);
  h_disc_th_rec = new TH2F("h_disc_th_rec", "DISC: #tehta_{C} vs #theta_{REC}; #theta_{REC}[GeV/c]; #theta_{C}(deg)", 200, 0, the_max, 200, 0, det_var_max[4]);

}

InitStatus EffHists::Init() {
  cout << "EffHists::Init" << endl;
  fAna = new PndAnalysis();
  init_tcas();
  init_hists();
  return kSUCCESS;
}

bool EffHists::check_prob_indiv(prob_func func, double cut, bool sans_dedx) {
  if ((m_prob_emcb->*func)(NULL)<cut) return false;
  if (!sans_dedx) {
    if ((m_prob_stt->*func)(NULL)<cut) return false;
    if ((m_prob_mvd->*func)(NULL)<cut) return false;
  }
  if ((m_prob_drc->*func)(NULL)<cut) return false;
  if ((m_prob_disc->*func)(NULL)<cut) return false;
  return true;
}

double EffHists::get_comb_prob(prob_func func, bool sans_dedx) {
  Double_t prob_emc = (m_prob_emcb->*func)(NULL);
  Double_t prob_stt = (m_prob_stt->*func)(NULL);
  Double_t prob_mvd = (m_prob_mvd->*func)(NULL);
  Double_t prob_drc = (m_prob_drc->*func)(NULL);
  Double_t prob_disc = (m_prob_disc->*func)(NULL);
  Double_t xx = 1.0;
  xx *= (prob_drc/(1.-prob_drc));
  xx *= (prob_disc/(1.-prob_disc));
  if (!sans_dedx) {
    xx *= (prob_mvd/(1.-prob_mvd));
    xx *= (prob_stt/(1.-prob_stt));
  }
  xx *= (prob_emc/(1.-prob_emc));
  return xx/(xx+1.);
}

void EffHists::Exec(Option_t* opt) {
  if (verb>1 or nevt%1000==0)
    cout << "EffHists::Exec evt " << nevt << " ======== " << endl;
  fAna->GetEvent(); // this may not be necessary
  nevt++;
  int ntrk = m_cand_array->GetEntriesFast();
  int ntrkmc = m_mc_array->GetEntriesFast();
  if (ntrkmc==0) {
    cout << "WARNING: number of mc track==0 " << endl;
    return;
  }

  PndMCTrack *truth = (PndMCTrack*) m_mc_array->At(0);
  Double_t mom_mc = truth->GetMomentum().Mag();
  Double_t the_mc = TMath::RadToDeg()*truth->GetMomentum().Theta();
  double mom_diff_min = 1e9;
  int itrk = -1;
  TVector3 t = truth->GetMomentum();
  for (int ii = 0; ii < ntrk; ++ii) {
    PndPidCandidate *cand = (PndPidCandidate*) m_cand_array->At(ii);
    TVector3 r = cand->GetMomentum();

    int imom = mom_mc<5?int(mom_mc/0.25):19;

    double dpx = r.X()-t.X();
    double dpy = r.Y()-t.Y();
    double dpz = r.Z()-t.Z();
    h_dpx->Fill(dpx, mom_mc);
    h_dpy->Fill(dpy, mom_mc);
    h_dpz->Fill(dpz, mom_mc);

    if ((dpx<dpx_min[imom]||dpx>dpx_max[imom])||
	(dpy<dpy_min[imom]||dpy>dpy_max[imom])||
	(dpz<dpz_min[imom]||dpz>dpz_max[imom])) {
      continue;
    }

    double mom_diff = hypot(hypot(r.X()-t.X(),r.Y()-t.Y()),r.Z()-t.Z());
    if (mom_diff<mom_diff_min) {
      mom_diff_min = mom_diff; itrk = ii;
    }

  }

  if (itrk == -1 ) {
    return;
  }

  PndPidCandidate *cand = (PndPidCandidate*) m_cand_array->At(itrk);
  mom_rec = cand->GetMomentum().Mag();
  the_rec = TMath::RadToDeg()*cand->GetMomentum().Theta();
  eoverp = cand->GetEmcCalEnergy()/mom_rec;
  stt_dedx = cand->GetSttMeanDEDX();
  disc_thetaC = TMath::RadToDeg()*cand->GetDiscThetaC();
  drc_thetaC = TMath::RadToDeg()*cand->GetDrcThetaC();
  muo_iron = cand->GetMuoIron();
  mvd_dedx = 1000*cand->GetMvdDEDX();

  h_emc_mom_mc->Fill(mom_mc,eoverp);
  h_stt_mom_mc->Fill(mom_mc,stt_dedx);
  h_mvd_mom_mc->Fill(mom_mc,mvd_dedx);
  h_dirc_mom_mc->Fill(mom_mc,drc_thetaC);
  h_disc_mom_mc->Fill(mom_mc,disc_thetaC);
  h_emc_mom_rec->Fill(mom_rec,eoverp);
  h_stt_mom_rec->Fill(mom_rec,stt_dedx);
  h_mvd_mom_rec->Fill(mom_rec,mvd_dedx);
  h_dirc_mom_rec->Fill(mom_rec,drc_thetaC);
  h_disc_mom_rec->Fill(mom_rec,disc_thetaC);
  h_emc_th_mc->Fill(the_mc,eoverp);
  h_stt_th_mc->Fill(the_mc,stt_dedx);
  h_mvd_th_mc->Fill(the_mc,mvd_dedx);
  h_dirc_th_mc->Fill(the_mc,drc_thetaC);
  h_disc_th_mc->Fill(the_mc,disc_thetaC);
  h_emc_th_rec->Fill(the_rec,eoverp);
  h_stt_th_rec->Fill(the_rec,stt_dedx);
  h_mvd_th_rec->Fill(the_rec,mvd_dedx);
  h_dirc_th_rec->Fill(the_rec,drc_thetaC);
  h_disc_th_rec->Fill(the_rec,disc_thetaC);

  m_prob_drc = (PndPidProbability*) m_drc_array->At(itrk);
  m_prob_disc = (PndPidProbability*) m_disc_array->At(itrk);
  m_prob_mvd = (PndPidProbability*) m_mvd_array->At(itrk);
  m_prob_stt = (PndPidProbability*) m_stt_array->At(itrk);
  m_prob_emcb = (PndPidProbability*) m_emcb_array->At(itrk);

  // combined probability of this track being a given type
  double prob_comb[npid_max]= {0.0};
  bool prob_indiv[npid_max]= {false};
  prob_comb[iel] = get_comb_prob(&PndPidProbability::GetElectronPidProb,false);
  prob_comb[imu] = get_comb_prob(&PndPidProbability::GetMuonPidProb,false);
  prob_comb[ipi] = get_comb_prob(&PndPidProbability::GetPionPidProb,false);
  prob_comb[ik]  = get_comb_prob(&PndPidProbability::GetKaonPidProb,false);
  prob_comb[iprot] = get_comb_prob(&PndPidProbability::GetProtonPidProb,false);

  prob_indiv[iel] = check_prob_indiv(&PndPidProbability::GetElectronPidProb,0.05,false);
  prob_indiv[imu] = check_prob_indiv(&PndPidProbability::GetMuonPidProb,0.05,false);
  prob_indiv[ipi] = check_prob_indiv(&PndPidProbability::GetPionPidProb,0.05,false);
  prob_indiv[ik]  = check_prob_indiv(&PndPidProbability::GetKaonPidProb,0.05,false);
  prob_indiv[iprot] = check_prob_indiv(&PndPidProbability::GetProtonPidProb,0.05,false);

  double prob_comb_sd[npid_max]= {0.0};
  double prob_indiv_sd[npid_max]= {false};
  prob_comb_sd[iel] = get_comb_prob(&PndPidProbability::GetElectronPidProb,true);
  prob_comb_sd[imu] = get_comb_prob(&PndPidProbability::GetMuonPidProb,true);
  prob_comb_sd[ipi] = get_comb_prob(&PndPidProbability::GetPionPidProb,true);
  prob_comb_sd[ik]  = get_comb_prob(&PndPidProbability::GetKaonPidProb,true);
  prob_comb_sd[iprot] = get_comb_prob(&PndPidProbability::GetProtonPidProb,true);

  prob_indiv_sd[iel] = check_prob_indiv(&PndPidProbability::GetElectronPidProb,0.05,false);
  prob_indiv_sd[imu] = check_prob_indiv(&PndPidProbability::GetMuonPidProb,0.05,false);
  prob_indiv_sd[ipi] = check_prob_indiv(&PndPidProbability::GetPionPidProb,0.05,false);
  prob_indiv_sd[ik]  = check_prob_indiv(&PndPidProbability::GetKaonPidProb,0.05,false);
  prob_indiv_sd[iprot] = check_prob_indiv(&PndPidProbability::GetProtonPidProb,0.05,false);

  fill_prob_hists(iel, &PndPidProbability::GetElectronPidProb, prob_comb[iel], prob_comb_sd[iel]);
  fill_prob_hists(ipi, &PndPidProbability::GetPionPidProb, prob_comb[ipi], prob_comb_sd[ipi]);
  fill_prob_hists(imu, &PndPidProbability::GetMuonPidProb, prob_comb[imu], prob_comb_sd[imu]);
  fill_prob_hists(ik, &PndPidProbability::GetKaonPidProb, prob_comb[ik], prob_comb_sd[ik]);
  fill_prob_hists(iprot, &PndPidProbability::GetProtonPidProb, prob_comb[iprot], prob_comb_sd[iprot]);

  for (int ipid=0; ipid<npid_max; ++ipid) {
    bool top = isTop(ipid, prob_comb);
    bool top_sd = isTop(ipid, prob_comb_sd);

    for (int ipc=0; ipc < prob_cut.size(); ++ipc) {
      if (ipc>=nprob_cut) break;

      eff1d_the[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],the_mc);
      eff1d_mom[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],mom_mc);
      eff2d[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],mom_mc,the_mc);

      eff1d_the_indiv[ipc][ipid]->Fill(prob_indiv[ipid]&&prob_comb[ipid]>prob_cut[ipc],the_mc);
      eff1d_mom_indiv[ipc][ipid]->Fill(prob_indiv[ipid]&&prob_comb[ipid]>prob_cut[ipc],mom_mc);
      eff2d_indiv[ipc][ipid]->Fill(prob_indiv[ipid]&&prob_comb[ipid]>prob_cut[ipc],mom_mc,the_mc);

      eff1d_the_sd[ipc][ipid]->Fill(prob_comb_sd[ipid]>prob_cut[ipc],the_mc);
      eff1d_mom_sd[ipc][ipid]->Fill(prob_comb_sd[ipid]>prob_cut[ipc],mom_mc);
      eff2d_sd[ipc][ipid]->Fill(prob_comb_sd[ipid]>prob_cut[ipc],mom_mc,the_mc);

      eff1d_the_sd_indiv[ipc][ipid]->Fill(prob_indiv_sd[ipid]&&prob_comb_sd[ipid]>prob_cut[ipc],the_mc);
      eff1d_mom_sd_indiv[ipc][ipid]->Fill(prob_indiv_sd[ipid]&&prob_comb_sd[ipid]>prob_cut[ipc],mom_mc);
      eff2d_sd_indiv[ipc][ipid]->Fill(prob_indiv_sd[ipid]&&prob_comb_sd[ipid]>prob_cut[ipc],mom_mc,the_mc);

    }
  }

}

void EffHists::fill_prob_hists(int ipid, prob_func func, double _prob_comb, double _prob_comb_sd) {
  h_prob[ipid][0]->Fill(_prob_comb,mom_rec);
  h_prob[ipid][1]->Fill(_prob_comb_sd,mom_rec);
  h_prob[ipid][2]->Fill(eoverp, (m_prob_emcb->*func)(NULL));
  h_prob[ipid][3]->Fill(stt_dedx, (m_prob_stt->*func)(NULL));
  h_prob[ipid][4]->Fill(mvd_dedx, (m_prob_mvd->*func)(NULL));
  h_prob[ipid][5]->Fill(drc_thetaC, (m_prob_drc->*func)(NULL));
  h_prob[ipid][6]->Fill(disc_thetaC, (m_prob_disc->*func)(NULL));
}

bool EffHists::isTop(int ipid, double prob[]) {
    for (int ii=0; ii < npid_max; ++ii) {
      if (ipid == ii) continue;
      if (prob[ii]>prob[ipid]) {
	return false;
      }
    }
    return true;
}

void EffHists::FinishTask() {
  cout << "EffHists::FinishTask" << endl;
  write_hists();
}

void EffHists::FinishEvent() {
  //cout << "EffHists::FinishEvent" << endl;
  //fAna->Reset();
}

void EffHists::set_prob_cut(int a_pid, double a_cut) {
  //if (a_pid<0 || a_pid >= npid_max) {
  //  cout << "EffHists::set_prob_cut a_pid= " << a_pid << " outside of allowed range [0," << npid_max << ")" << endl;
  //  return;
  //}
  //if (a_cut<0 || a_cut >= 1) {
  //  cout << "EffHists::set_prob_cut a_cut= " << a_cut << " outside of sensible range [0,1]. consider using a value between 0 and 1" << endl;
  //  return;
  //}
  //prob_cut[a_pid] = a_cut;
}

void EffHists::write_hists() {

  const char* root = gDirectory->GetPath();

  h_dpx->Write();
  h_dpy->Write();
  h_dpz->Write();

  gDirectory->mkdir("probs");
  gDirectory->cd("probs");
  for (int ipid=0; ipid < npid_max; ++ipid) {
    for (int ihist=0; ihist < 7; ++ihist) {
      h_prob[ipid][ihist]->Write();
    }

  }
  gDirectory->cd(root);

  gDirectory->mkdir("detvars");
  gDirectory->cd("detvars");
  h_emc_mom_mc->Write();
  h_stt_mom_mc->Write();
  h_mvd_mom_mc->Write();
  h_dirc_mom_mc->Write();
  h_disc_mom_mc->Write();
  h_emc_mom_rec->Write();
  h_stt_mom_rec->Write();
  h_mvd_mom_rec->Write();
  h_dirc_mom_rec->Write();
  h_disc_mom_rec->Write();
  h_emc_th_mc->Write();
  h_stt_th_mc->Write();
  h_mvd_th_mc->Write();
  h_dirc_th_mc->Write();
  h_disc_th_mc->Write();
  h_emc_th_rec->Write();
  h_stt_th_rec->Write();
  h_mvd_th_rec->Write();
  h_dirc_th_rec->Write();
  h_disc_th_rec->Write();
  gDirectory->cd(root);

  for (int ipc=0; ipc < prob_cut.size(); ++ipc) {
    if (ipc>=nprob_cut) break;
    const char* subdir = Form("prob_cut_%d",ipc);
    gDirectory->mkdir(subdir);
    gDirectory->cd(subdir);
    for (int ipid=0; ipid<npid_max; ++ipid) {
      eff1d_the[ipc][ipid]->Write();
      eff1d_mom[ipc][ipid]->Write();
      eff2d[ipc][ipid]->Write();

      eff1d_the_indiv[ipc][ipid]->Write();
      eff1d_mom_indiv[ipc][ipid]->Write();
      eff2d_indiv[ipc][ipid]->Write();

      eff1d_the_sd[ipc][ipid]->Write();
      eff1d_mom_sd[ipc][ipid]->Write();
      eff2d_sd[ipc][ipid]->Write();

      eff1d_the_sd_indiv[ipc][ipid]->Write();
      eff1d_mom_sd_indiv[ipc][ipid]->Write();
      eff2d_sd_indiv[ipc][ipid]->Write();
    }
    gDirectory->cd(root);
  }

}

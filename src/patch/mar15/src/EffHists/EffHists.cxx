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
  if ((m_prob_drc->*func)(NULL)>cut) return false;
  if ((m_prob_disc->*func)(NULL)>cut) return false;
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
  for (int ii = 0; ii < ntrk; ++ii) {
    PndPidCandidate *cand = (PndPidCandidate*) m_cand_array->At(ii);
    TVector3 r = cand->GetMomentum();
    TVector3 t = truth->GetMomentum();
    double mom_diff = hypot(hypot(r.X()-t.X(),r.Y()-t.Y()),r.Z()-t.Z());
    if (mom_diff<mom_diff_min) { mom_diff_min = mom_diff; itrk = ii;}
  }
  if (itrk == -1 ) {
    return;
  }

  PndPidCandidate *cand = (PndPidCandidate*) m_cand_array->At(itrk);
  Double_t mom_rec = cand->GetMomentum().Mag();
  Double_t the_rec = TMath::RadToDeg()*cand->GetMomentum().Theta();
  Double_t eoverp = cand->GetEmcCalEnergy()/mom_rec;
  Double_t stt_dedx = cand->GetSttMeanDEDX();
  Double_t disc_thetaC = TMath::RadToDeg()*cand->GetDiscThetaC();
  Double_t drc_thetaC = TMath::RadToDeg()*cand->GetDrcThetaC();
  Double_t muo_iron = cand->GetMuoIron();
  Double_t mvd_dedx = 1000*cand->GetMvdDEDX();

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
  double prob_indiv_sd[npid_max]= {0.0};
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

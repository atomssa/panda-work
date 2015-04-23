#include "EffHists.h"

#include "FairTask.h"

#include "PndPidProbability.h"
#include "PndPidCandidate.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

const TString s_spc[EffHists::nsp_max] = {"posit","muplus","piplus","kplus","prot","elec","muminus","piminus","kminus","antiprot"};
const TString s_spc_tex[EffHists::nsp_max] = {"e^{+}","#mu^{+}","#pi^{+}","K^{+}","p","e^{-}","#mu^{-}","#pi^{-}","K^{-}","#bar{p}"};
const TString s_pid[EffHists::npid_max] = {"e_id", "mu_id", "pi_id", "k_id", "prot_id"};
const TString s_det[EffHists::ndet] = {"emc", "stt", "mvd", "dirc", "disc"};

EffHists::EffHists(int a_sp):
  m_sp(a_sp),
  verb(false),
  out_file_name("eff_hists.root"),
  prob_cut{0.5, 0.5, 0.5, 0.5, 0.5},
  det_var_max{1.5, 4, 4, 90, 90},
  mom_max(2.0),
  the_max(180.0),
  fAna()
{
}

EffHists::~EffHists() {
}

InitStatus EffHists::init_tca(TClonesArray *tca, TString name) {
  tca = dynamic_cast<TClonesArray *> (m_ioman->GetObject(name));
  if ( ! tca ) {
    cout << "-W- EffHists::Init: "
	 << "No " << name << " array!" << endl;
    return kERROR;
  }
  return kSUCCESS;
}

InitStatus EffHists::init_tcas() {
  m_ioman = FairRootManager::Instance();
  if ( ! m_ioman ){
    cout << "-E- BremPidReader::Init: "
	 << "RootManager not instantiated!" << endl;
    return kFATAL;
  }
  if ( init_tca(m_cand_array, "PidChargedCand") != kSUCCESS) return kERROR;
  if ( init_tca(m_drc_array, "PidAlgoDrc") != kSUCCESS) return kERROR;
  if ( init_tca(m_disc_array, "PidAlgoDisc") != kSUCCESS) return kERROR;
  if ( init_tca(m_stt_array, "PidAlgoEmcStt") != kSUCCESS) return kERROR;
  if ( init_tca(m_mvd_array, "PidAlgoEmcMvd") != kSUCCESS) return kERROR;
  if ( init_tca(m_emcb_array, "PidAlgoEmcBayes") != kSUCCESS) return kERROR;
  return kSUCCESS;
}

void EffHists::init_hists() {
  int nbin =100;
  for (int ipid = 0; ipid < npid_max; ++ipid) {
    TString title = Form("Eff %s to pass %s cuts at prob>%4.2f",s_spc[m_sp].Data(),s_pid[ipid].Data(),prob_cut[ipid]);
    eff_den[ipid] = new TH2F(Form("eff_den_%s",s_pid[ipid].Data()),title+";mom[GeV/c];#theta[rad]",nbin,0,mom_max,nbin,0,the_max);
    eff_num[ipid] = new TH2F(Form("eff_num_%s",s_pid[ipid].Data()),title+";mom[GeV/c];#theta[rad]",nbin,0,mom_max,nbin,0,the_max);
    eff1d_mom[ipid] = new TEfficiency(Form("eff1d_mom_%s",s_pid[ipid].Data()), title+";mom[GeV/c]", nbin, 0, mom_max);
    eff1d_the[ipid] = new TEfficiency(Form("eff1d_the_%s",s_pid[ipid].Data()), title+";#theta[rad]", nbin, 0, the_max);
    eff2d[ipid] = new TEfficiency(Form("eff2d_%s",s_pid[ipid].Data()), title+";mom[GeV/c];#theta[rad]", nbin, 0, mom_max, nbin, 0, the_max);
  }
}

InitStatus EffHists::Init() {
  cout << "EffHists::Init" << endl;
  fAna = new PndAnalysis();
  init_hists();
  return kSUCCESS;
}

double EffHists::get_comb_prob(prob_func func) {
  Double_t prob_emc = (m_prob_emcb->*func)(NULL);
  Double_t prob_stt = (m_prob_stt->*func)(NULL);
  Double_t prob_mvd = (m_prob_mvd->*func)(NULL);
  Double_t prob_drc = (m_prob_drc->*func)(NULL);
  Double_t prob_disc = (m_prob_disc->*func)(NULL);
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double EffHists::get_comb_prob_elec() {
  Double_t prob_emc = m_prob_emcb->GetElectronPidProb();
  Double_t prob_stt = m_prob_stt->GetElectronPidProb();
  Double_t prob_mvd = m_prob_mvd->GetElectronPidProb();
  Double_t prob_drc = m_prob_drc->GetElectronPidProb();
  Double_t prob_disc = m_prob_disc->GetElectronPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double EffHists::get_comb_prob_pion() {
  Double_t prob_emc = m_prob_emcb->GetPionPidProb();
  Double_t prob_stt = m_prob_stt->GetPionPidProb();
  Double_t prob_mvd = m_prob_mvd->GetPionPidProb();
  Double_t prob_drc = m_prob_drc->GetPionPidProb();
  Double_t prob_disc = m_prob_disc->GetPionPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double EffHists::get_comb_prob_proton() {
  Double_t prob_emc = m_prob_emcb->GetProtonPidProb();
  Double_t prob_stt = m_prob_stt->GetProtonPidProb();
  Double_t prob_mvd = m_prob_mvd->GetProtonPidProb();
  Double_t prob_drc = m_prob_drc->GetProtonPidProb();
  Double_t prob_disc = m_prob_disc->GetProtonPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double EffHists::get_comb_prob_kaon() {
  Double_t prob_emc = m_prob_emcb->GetKaonPidProb();
  Double_t prob_stt = m_prob_stt->GetKaonPidProb();
  Double_t prob_mvd = m_prob_mvd->GetKaonPidProb();
  Double_t prob_drc = m_prob_drc->GetKaonPidProb();
  Double_t prob_disc = m_prob_disc->GetKaonPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double EffHists::get_comb_prob_muon() {
  Double_t prob_emc = m_prob_emcb->GetMuonPidProb();
  Double_t prob_stt = m_prob_stt->GetMuonPidProb();
  Double_t prob_mvd = m_prob_mvd->GetMuonPidProb();
  Double_t prob_drc = m_prob_drc->GetMuonPidProb();
  Double_t prob_disc = m_prob_disc->GetMuonPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

void EffHists::Exec(Option_t* opt) {
  if (verb>1 or nevt%100==0)
    cout << "======== EffHists::Exec evt " << nevt << " ======== " << endl;
  fAna->GetEvent();
  nevt++;
  int ntrk = m_cand_array->GetEntriesFast();
  for (int itrk = 0; itrk < ntrk; ++itrk) {
    PndPidCandidate *cand = (PndPidCandidate*) m_cand_array->At(itrk);
    Double_t mom = cand->GetMomentum().Mag();
    Double_t the = TMath::RadToDeg()*cand->GetMomentum().Theta();
    Double_t eoverp = cand->GetEmcCalEnergy()/mom;
    Double_t stt_dedx = cand->GetSttMeanDEDX();
    Double_t disc_thetaC = TMath::RadToDeg()*cand->GetDiscThetaC();
    Double_t drc_thetaC = TMath::RadToDeg()*cand->GetDrcThetaC();
    Double_t muo_iron = cand->GetMuoIron();
    Double_t mvd_dedx = 1000*cand->GetMvdDEDX();

    m_prob_drc = (PndPidProbability*) m_drc_array->At(itrk);
    m_prob_disc = (PndPidProbability*) m_disc_array->At(itrk);
    m_prob_mvd = (PndPidProbability*) m_mvd_array->At(itrk);
    m_prob_stt = (PndPidProbability*) m_stt_array->At(itrk);
    m_prob_emcb = (PndPidProbability*) m_emcb_array->At(itrk);

    // combined probability of this track being a given type
    double prob_comb[npid_max]= {0.0};
    //prob_comb[iel] = get_comb_prob_elec();
    //prob_comb[imu] = get_comb_prob_muon();
    //prob_comb[ipi] = get_comb_prob_pion();
    //prob_comb[ik]  = get_comb_prob_kaon();
    //prob_comb[iprot] = get_comb_prob_proton();

    prob_comb[iel] = get_comb_prob(&PndPidProbability::GetElectronPidProb);
    prob_comb[imu] = get_comb_prob(&PndPidProbability::GetMuonPidProb);
    prob_comb[ipi] = get_comb_prob(&PndPidProbability::GetPionPidProb);
    prob_comb[ik]  = get_comb_prob(&PndPidProbability::GetKaonPidProb);
    prob_comb[iprot] = get_comb_prob(&PndPidProbability::GetProtonPidProb);
    //prob_comb[iposit]  = get_comb_prob(&PndPidProbability::GetElectronPidProb);

    for (int ipid=0; ipid<npid_max; ++ipid) {
      eff_den[ipid]->Fill(mom,the);
      if (prob_comb[ipid]>prob_cut[ipid]) { eff_num[ipid]->Fill(mom, the); }
      eff1d_the[ipid]->Fill(prob_comb[ipid]>prob_cut[ipid],the);
      eff1d_mom[ipid]->Fill(prob_comb[ipid]>prob_cut[ipid],mom);
      eff2d[ipid]->Fill(prob_comb[ipid]>prob_cut[ipid],mom,the);
    }
  }

}

void EffHists::FinishTask() {
  cout << "EffHists::FinishTask" << endl;
  write_hists();
}

void EffHists::FinishEvent() {
  cout << "EffHists::Finish" << endl;
  fAna->Reset();
}

void EffHists::write_hists() {
  for (int ipid=0; ipid<npid_max; ++ipid) {
    eff_den[ipid]->Write();
    eff_num[ipid]->Write();
    eff1d_the[ipid]->Write();
    eff1d_mom[ipid]->Write();
    eff2d[ipid]->Write();
  }
}

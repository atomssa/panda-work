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
  det_var_max{1.5, 4, 4, 90, 90},
  mom_max(10.0),
  the_max(180.0),
  fAna()
{
  for (int iprob_cut=0; iprob_cut < nprob_cut; ++iprob_cut) {
    prob_cut[iprob_cut] = (1.0/double(nprob_cut-1))*double(iprob_cut);
    //cout << "prob_cut[" << iprob_cut << "] = " << prob_cut[iprob_cut] << endl;
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
    eff_den_mc[ipid] = new TH2F(Form("eff_den_%s",s_pid[ipid].Data()),title+" (MC);p_{MC}[GeV/c];#theta_{MC}[rad]",nbin,0,mom_max,nbin,0,the_max);
    eff_den_rec[ipid] = new TH2F(Form("rec_eff_den_%s",s_pid[ipid].Data()),title+" (REC);p_{REC}[GeV/c];#theta_{REC}[rad]",nbin,0,mom_max,nbin,0,the_max);

    for (int ipc=0; ipc < nprob_cut; ++ipc) {

      title = Form("%s passing %s cuts with prob>%4.2f",s_spc_tex[m_sp].Data(),s_pid[ipid].Data(),prob_cut[ipc]);
      eff_num_mc[ipc][ipid] = new TH2F(Form("eff_num_%s_ipc%d",s_pid[ipid].Data(),ipc),title+" (MC);p_{MC}[GeV/c];#theta_{MC}[rad]",nbin,0,mom_max,nbin,0,the_max);
      eff_num_rec[ipc][ipid] = new TH2F(Form("rec_eff_num_%s_ipc%d",s_pid[ipid].Data(),ipc),title+" (REC);p_{REC}[GeV/c];#theta_{REC}[rad]",nbin,0,mom_max,nbin,0,the_max);

      title = Form("efficiency of %s to pass %s cuts at prob>%4.2f",s_spc_tex[m_sp].Data(),s_pid[ipid].Data(),prob_cut[ipc]);
      eff1d_mom_mc[ipc][ipid] = new TEfficiency(Form("eff1d_mom_%s_ipc%d",s_pid[ipid].Data(),ipc), title+" (MC);p_{MC}[GeV/c]", nbin, 0, mom_max);
      eff1d_the_mc[ipc][ipid] = new TEfficiency(Form("eff1d_the_%s_ipc%d",s_pid[ipid].Data(),ipc), title+" (MC);#theta_{MC}[rad]", nbin, 0, the_max);
      eff2d_mc[ipc][ipid] = new TEfficiency(Form("eff2d_%s_ipc%d",s_pid[ipid].Data(),ipc), title+" (MC);p_{MC}[GeV/c];#theta_{MC}[rad]", nbin, 0, mom_max, nbin, 0, the_max);

      eff1d_mom_rec[ipc][ipid] = new TEfficiency(Form("rec_eff1d_mom_%s_ipc%d",s_pid[ipid].Data(),ipc), title+" (REC);p_{REC}[GeV/c]", nbin, 0, mom_max);
      eff1d_the_rec[ipc][ipid] = new TEfficiency(Form("rec_eff1d_the_%s_ipc%d",s_pid[ipid].Data(),ipc), title+" (REC);#theta_{REC}[rad]", nbin, 0, the_max);
      eff2d_rec[ipc][ipid] = new TEfficiency(Form("rec_eff2d_%s_ipc%d",s_pid[ipid].Data(),ipc), title+" (REC);p_{REC}[GeV/c];#theta_{REC}[rad]", nbin, 0, mom_max, nbin, 0, the_max);

    }
  }

  h_emc_mc = new TH2F("h_emc_mc", "EMC: E/p vs mom_{MC}; p_{MC}[GeV/c]; E/p", 200, 0, mom_max, 200, 0, det_var_max[0]);
  h_stt_mc = new TH2F("h_stt_mc", "STT: dedx vs mom_{MC}; p_{MC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[1]);
  h_mvd_mc = new TH2F("h_mvd_mc", "MVD: dedx vs mom_{MC}; p_{MC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[2]);
  h_dirc_mc = new TH2F("h_dirc_mc", "DIRC: #tehta_{C} vs mom_{MC}; p_{MC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[3]);
  h_disc_mc = new TH2F("h_disc_mc", "DISC: #tehta_{C} vs mom_{MC}; p_{MC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[4]);

  h_emc_rec = new TH2F("h_emc_rec", "EMC: E/p vs mom_{REC}; p_{REC}[GeV/c]; E/p", 200, 0, mom_max, 200, 0, det_var_max[0]);
  h_stt_rec = new TH2F("h_stt_rec", "STT: dedx vs mom_{REC}; p_{REC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[1]);
  h_mvd_rec = new TH2F("h_mvd_rec", "MVD: dedx vs mom_{REC}; p_{REC}[GeV/c]; dE/dx", 200, 0, mom_max, 200, 0, det_var_max[2]);
  h_dirc_rec = new TH2F("h_dirc_rec", "DIRC: #tehta_{C} vs mom_{REC}; p_{REC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[3]);
  h_disc_rec = new TH2F("h_disc_rec", "DISC: #tehta_{C} vs mom_{REC}; p_{REC}[GeV/c]; #theta_{C}(deg)", 200, 0, mom_max, 200, 0, det_var_max[4]);

}

InitStatus EffHists::Init() {
  cout << "EffHists::Init" << endl;
  fAna = new PndAnalysis();
  init_tcas();
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

void EffHists::Exec(Option_t* opt) {
  if (verb>1 or nevt%1000==0)
    cout << "EffHists::Exec evt " << nevt << " ======== " << endl;
  fAna->GetEvent(); // this may not be necessary
  nevt++;
  int ntrk = m_cand_array->GetEntriesFast();
  int ntrkmc = m_mc_array->GetEntriesFast();
  if (ntrkmc!=1) {
    cout << "WARNING: number of mc track!= 1" << endl;
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

  h_emc_mc->Fill(mom_mc,eoverp);
  h_stt_mc->Fill(mom_mc,stt_dedx);
  h_mvd_mc->Fill(mom_mc,mvd_dedx);
  h_dirc_mc->Fill(mom_mc,drc_thetaC);
  h_disc_mc->Fill(mom_mc,disc_thetaC);

  h_emc_rec->Fill(mom_rec,eoverp);
  h_stt_rec->Fill(mom_rec,stt_dedx);
  h_mvd_rec->Fill(mom_rec,mvd_dedx);
  h_dirc_rec->Fill(mom_rec,drc_thetaC);
  h_disc_rec->Fill(mom_rec,disc_thetaC);

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
    eff_den_mc[ipid]->Fill(mom_mc,the_mc);
    eff_den_rec[ipid]->Fill(mom_rec,the_rec);
    for (int ipc=0; ipc < nprob_cut; ++ipc) {
      if (prob_comb[ipid]>prob_cut[ipc]) {
	eff_num_mc[ipc][ipid]->Fill(mom_mc, the_mc);
	eff_num_rec[ipc][ipid]->Fill(mom_rec, the_rec);
      }
      eff1d_the_mc[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],the_mc);
      eff1d_mom_mc[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],mom_mc);
      eff2d_mc[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],mom_mc,the_mc);

      eff1d_the_rec[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],the_rec);
      eff1d_mom_rec[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],mom_rec);
      eff2d_rec[ipc][ipid]->Fill(prob_comb[ipid]>prob_cut[ipc],mom_rec,the_rec);

    }
  }

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

  h_emc_mc->Write();
  h_stt_mc->Write();
  h_mvd_mc->Write();
  h_dirc_mc->Write();
  h_disc_mc->Write();
  h_emc_rec->Write();
  h_stt_rec->Write();
  h_mvd_rec->Write();
  h_dirc_rec->Write();
  h_disc_rec->Write();


  for (int ipid=0; ipid<npid_max; ++ipid) {
    eff_den_mc[ipid]->Write();
    eff_den_rec[ipid]->Write();
  }

  const char* root = gDirectory->GetPath();
  for (int ipc=0; ipc < nprob_cut; ++ipc) {
    const char* subdir = Form("prob_cut_%d",ipc);
    gDirectory->mkdir(subdir);
    gDirectory->cd(subdir);
    for (int ipid=0; ipid<npid_max; ++ipid) {
      eff_num_mc[ipc][ipid]->Write();
      eff1d_the_mc[ipc][ipid]->Write();
      eff1d_mom_mc[ipc][ipid]->Write();
      eff2d_mc[ipc][ipid]->Write();

      eff_num_rec[ipc][ipid]->Write();
      eff1d_the_rec[ipc][ipid]->Write();
      eff1d_mom_rec[ipc][ipid]->Write();
      eff2d_rec[ipc][ipid]->Write();
    }
    gDirectory->cd(root);
  }

}

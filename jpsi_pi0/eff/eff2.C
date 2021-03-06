typedef PndPidProbability* ppp_p;

double get_comb_prob(ppp_p _prob_emc, ppp_p _prob_stt, ppp_p _prob_mvd,
		     ppp_p _prob_drc, ppp_p _prob_disc
		     Double_t(PndPidProbability::*)(PndPidProbability*) func) {
  Double_t prob_emc = (_prob_emc->*func)();
  Double_t prob_stt = (_prob_stt->*func)();
  Double_t prob_mvd = (_prob_mvd->*func)();
  Double_t prob_drc = (_prob_drc->*func)();
  Double_t prob_disc = (_prob_disc->*func)();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double get_comb_prob_elec(ppp_p _prob_emc, ppp_p _prob_stt, ppp_p _prob_mvd,
			  ppp_p _prob_drc, ppp_p _prob_disc) {
  Double_t prob_emc = _prob_emc->GetElectronPidProb();
  Double_t prob_stt = _prob_stt->GetElectronPidProb();
  Double_t prob_mvd = _prob_mvd->GetElectronPidProb();
  Double_t prob_drc = _prob_drc->GetElectronPidProb();
  Double_t prob_disc = _prob_disc->GetElectronPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double get_comb_prob_pion(ppp_p _prob_emc, ppp_p _prob_stt, ppp_p _prob_mvd,
			  ppp_p _prob_drc, ppp_p _prob_disc) {
  Double_t prob_emc = _prob_emc->GetPionPidProb();
  Double_t prob_stt = _prob_stt->GetPionPidProb();
  Double_t prob_mvd = _prob_mvd->GetPionPidProb();
  Double_t prob_drc = _prob_drc->GetPionPidProb();
  Double_t prob_disc = _prob_disc->GetPionPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double get_comb_prob_proton(ppp_p _prob_emc, ppp_p _prob_stt, ppp_p _prob_mvd,
			  ppp_p _prob_drc, ppp_p _prob_disc) {
  Double_t prob_emc = _prob_emc->GetProtonPidProb();
  Double_t prob_stt = _prob_stt->GetProtonPidProb();
  Double_t prob_mvd = _prob_mvd->GetProtonPidProb();
  Double_t prob_drc = _prob_drc->GetProtonPidProb();
  Double_t prob_disc = _prob_disc->GetProtonPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double get_comb_prob_kaon(ppp_p _prob_emc, ppp_p _prob_stt, ppp_p _prob_mvd,
			  ppp_p _prob_drc, ppp_p _prob_disc) {
  Double_t prob_emc = _prob_emc->GetKaonPidProb();
  Double_t prob_stt = _prob_stt->GetKaonPidProb();
  Double_t prob_mvd = _prob_mvd->GetKaonPidProb();
  Double_t prob_drc = _prob_drc->GetKaonPidProb();
  Double_t prob_disc = _prob_disc->GetKaonPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

double get_comb_prob_muon(ppp_p _prob_emc, ppp_p _prob_stt, ppp_p _prob_mvd,
			  ppp_p _prob_drc, ppp_p _prob_disc) {
  Double_t prob_emc = _prob_emc->GetMuonPidProb();
  Double_t prob_stt = _prob_stt->GetMuonPidProb();
  Double_t prob_mvd = _prob_mvd->GetMuonPidProb();
  Double_t prob_drc = _prob_drc->GetMuonPidProb();
  Double_t prob_disc = _prob_disc->GetMuonPidProb();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

void eff2(int sp /*species index*/, int ifile ) {

  const double prob_cut = 0.5;
  const double mom_max = 2.0;
  const double the_max = 180.0;

  const TString basedir="/vol0/panda/work/jpsi_pi0/grid.out/jacek/";
  const TString subdir = Form("/runall.%d/",ifile);


  enum{iposit=0, imuonp, ipionp, ikaonp, iproton, ielec, imuonm, ipionm, ikaonm, iantiproton, nsp_max};
  const TString species[nsp_max] = {"posit","muplus","piplus","kplus","prot","elec","muminus","piminus","kminus","antiprot"};

  enum{iel=0,imu,ipi,ik,iprot,npid_max};
  const TString s_pid[npid_max] = {"e_id", "mu_id", "pi_id", "k_id", "prot_id"};

  enum{iemc = 0, istt, imvd, idirc, idisc, ndet};
  const TString det[ndet] = {"emc", "stt", "mvd", "dirc", "disc"};
  const double det_var_max[ndet] = {1.5, 4, 4, 90, 90};

  TString fName = basedir+species[sp]+subdir+"pid_complete.root";
  cout << "opening " << fName << " in read mode... " << endl;
  TFile *inFile = new TFile(fName,"READ");
  TTree *cbmsim = (TTree *) inFile->Get("cbmsim");

  TClonesArray* cCand_array=new TClonesArray("PndPidCandidate");
  cbmsim->SetBranchAddress("PidChargedCand", &cCand_array);

  TClonesArray* drc_array=new TClonesArray("PndPidProbability");
  cbmsim->SetBranchAddress("PidAlgoDrc", &drc_array);

  TClonesArray* disc_array=new TClonesArray("PndPidProbability");
  cbmsim->SetBranchAddress("PidAlgoDisc", &disc_array);

  TClonesArray* mvd_array=new TClonesArray("PndPidProbability");
  cbmsim->SetBranchAddress("PidAlgoMvd", &mvd_array);

  TClonesArray* stt_array=new TClonesArray("PndPidProbability");
  cbmsim->SetBranchAddress("PidAlgoStt", &stt_array);

  TClonesArray* emcb_array=new TClonesArray("PndPidProbability");
  cbmsim->SetBranchAddress("PidAlgoEmcBayes", &emcb_array);

  int nevt = cbmsim->GetEntries();
  cout << "nevt = " << nevt << endl;

  TH2F* prob[nsp_max][ndet];
  for (int idet=0; idet<ndet; ++idet) {
    for (int isp=0; isp<nsp_max; ++isp) {
      TString name = Form("%s_prob_%s", det[idet].Data(), species[isp].Data());
      TString title = Form("%s_prob_%s", det[idet].Data(), species[isp].Data());
      cout << "hist name  = "  << name << endl;
      prob[isp][idet] = new TH2F(name, title, 200, 0, mom_max, 200, det_var_max[idet]);
    }
  }

  TH2F* eff_den[nsp_max], *eff_num[nsp_max];
  TEfficiency *eff1d_the[nsp_max], eff1d_mom[nsp_max];
  for (int isp = 0; isp < nsp_max; ++isp) {
    eff_den[isp] = new TH2F(Form("eff_den_%s",species[isp].Data()),Form("eff_den_%s",species[isp].Data()),200,0,mom_max,200,0,the_max);
    eff_num[isp] = new TH2F(Form("eff_den_%s",species[isp].Data()),Form("eff_den_%s",species[isp].Data()),200,0,mom_max,200,0,the_max);
    eff1d_mom[isp] = new TEfficiency(Form("eff1d_the_%s",species[isp].Data()), Form("eff1d_the_%s",species[isp].Data()), 200, 0, the_max);
    eff1d_the[isp] = new TEfficiency(Form("eff1d_mom_%s",species[isp].Data()), Form("eff1d_mom_%s",species[isp].Data()), 200, 0, the_max);
  }

  for (int ievt = 0; ievt < 1000; ++ievt) {
    cout << "ievt= " << ievt << " / " << nevt <<endl;

    int ntrk = cbmsim->GetEntriesFast();
    for (int itrk = 0; itrk < ntrk; ++itrk) {

      PndPidCandidate *cand = (PndPidCandidate*) cCand_array->At(itrk);

      Double_t mom = cand->GetMomentum().Mag();
      Double_t the = TMath::RadToDeg()*cand->GetMomentum().Theta();
      Double_t eoverp = cand->GetEmcCalEnergy()/mom;
      Double_t stt_dedx = cand->GetSttMeanDEDX();
      Double_t disc_thetaC = TMath::RadToDeg()*cand->GetDiscThetaC();
      Double_t drc_thetaC = TMath::RadToDeg()*cand->GetDrcThetaC();
      Double_t muo_iron = cand->GetMuoIron();
      Double_t mvd_dedx = 1000*cand->GetMvdDEDX();

      PndPidProbability *prob_drc = (PndPidProbability*) drc_array->At(itrk);
      PndPidProbability *prob_disc = (PndPidProbability*) disc_array->At(itrk);
      PndPidProbability *prob_mvd = (PndPidProbability*) mvd_array->At(itrk);
      PndPidProbability *prob_stt = (PndPidProbability*) stt_array->At(itrk);
      PndPidProbability *prob_emc = (PndPidProbability*) emcb_array->At(itrk);

      // combined probability of this track being a
      double prob_comb[nsp_max] = {0.0};

      double prob_comb[ielec] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetElectronPidProb);
      double prob_comb[imuonm] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetMuonPidProb);
      double prob_comb[ipionm] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetPionPidProb);
      double prob_comb[ikaonm] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetKaonPidProb);
      double prob_comb[iprotm] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetProtonPidProb);

      eff_den->Fill(mom,the);
      if (prob_comb[ielec]>0.99) { eff_num[sp]->Fill(mom, the); }
      eff1d_the->Fill(the, prob_comb[ielec]>0.99);
      eff1d_mom->Fill(the, prob_comb[ielec]>0.99);

      //Double_t prob_drc_e = prob_drc->GetElectronPidProb();
      //Double_t prob_disc_e = prob_disc->GetElectronPidProb();
      //Double_t prob_mvd_e = prob_mvd->GetElectronPidProb();
      //Double_t prob_stt_e = prob_stt->GetElectronPidProb();
      //Double_t prob_emc_e = prob_emc->GetElectronPidProb();
      //Double_t xx_e = (prob_drc_e/(1-prob_drc_e))*(prob_disc_e/(1-prob_disc_e))
      //	*(prob_mvd_e/(1-prob_mvd_e))*(prob_stt_e/(1-prob_stt_e))
      //	*(prob_emc_e/(1-prob_emc_e));
      //Double_t prob_comb_e = xx_e/(xx_e+1);

      //prob[sp][iemc]->Fill(mom, eoverp, prob_emc_e);
      //prob[sp][istt]->Fill(mom, stt_dedx);
      //prob[sp][imvd]->Fill(mom, mvd_dedx);
      //prob[sp][idisc]->Fill(mom, disc_thetaC);
      //prob[sp][idrc]->Fill(mom, drc_thetaC);
      //prob[sp][imuo]->Fill(mom, muo_iron);

      //Double_t prob_drc_pi = prob_drc->GetPionPidProb();
      //Double_t prob_disc_pi = prob_disc->GetPionPidProb();
      //Double_t prob_mvd_pi = prob_mvd->GetPionPidProb();
      //Double_t prob_stt_pi = prob_stt->GetPionPidProb();
      //Double_t prob_emcb_pi = prob_emc->GetPionPidProb();
      //Double_t xx_pi = (prob_drc_pi/(1-prob_drc_pi))*(prob_disc_pi/(1-prob_disc_pi))
      //	*(prob_mvd_pi/(1-prob_mvd_pi))*(prob_stt_pi/(1-prob_stt_pi))
      //	*(prob_pimcb_pi/(1-prob_pimcb_pi));
      //Double_t prob_comb_pi = xx_pi/(xx_pi+1);
      //
      //Double_t prob_drc_prot = prob_drc->GetProtonPidProb();
      //Double_t prob_disc_prot = prob_disc->GetProtonPidProb();
      //Double_t prob_mvd_prot = prob_mvd->GetProtonPidProb();
      //Double_t prob_stt_prot = prob_stt->GetProtonPidProb();
      //Double_t prob_emcb_prot = prob_emc->GetProtonPidProb();
      //Double_t xx_prot = (prob_drc_prot/(1-prob_drc_prot))*(prob_disc_prot/(1-prob_disc_prot))
      //	*(prob_mvd_prot/(1-prob_mvd_prot))*(prob_stt_prot/(1-prob_stt_prot))
      //	*(prob_protmcb_prot/(1-prob_protmcb_prot));
      //Double_t prob_comb_prot = xx_prot/(xx_prot+1);

    }

  }

  TFile *fout = new TFile("eff2_out.root","RECREATE");
  fout->cd();
  for (int ipid=0; ipid<npid_max; ++ipid) {
    eff_den[ipid]->Write();
    eff_num[ipid]->Write();
    eff1d_the[ipid]->Write();
    eff1d_mom[ipid]->Write();
    eff2d[ipid]->Write();
  }

  TCanvas *tc= new TCanvas("tc", "tc");
  tc->Divide(3,2);
  for (int ipid=0; ipid<npid_max; ++ipid) {
    //eff_den[ipid]->Draw();
    //eff_num[ipid]->Write();
    //eff1d_the[ipid]->Write();
    //eff1d_mom[ipid]->Write();
    tc->cd(1+ipid);
    eff2d[ipid]->Draw("colz");
  }
  tc->Write();

}

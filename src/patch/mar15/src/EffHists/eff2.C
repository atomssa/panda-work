typedef PndPidProbability* ppp_p;

double get_comb_prob(ppp_p _prob_emc, ppp_p _prob_stt, ppp_p _prob_mvd,
		     ppp_p _prob_drc, ppp_p _prob_disc,
		     auto func) {
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

  //TH2F* prob[nsp_max][ndet];
  //for (int idet=0; idet<ndet; ++idet) {
  //  for (int isp=0; isp<nsp_max; ++isp) {
  //    TString name = Form("%s_prob_%s", det[idet].Data(), species[isp].Data());
  //    TString title = Form("%s_prob_%s", det[idet].Data(), species[isp].Data());
  //    prob[isp][idet] = new TH2F(name, title, 200, 0, mom_max, 200, det_var_max[idet]);
  //  }
  //}

  int nbin =100;
  TH2F* eff_den[npid_max], *eff_num[npid_max];
  TEfficiency *eff2d[npid_max], *eff1d_the[npid_max], *eff1d_mom[npid_max];
  for (int ipid = 0; ipid < npid_max; ++ipid) {
    TString title = Form("Eff %s to pass %s cuts at prob>%4.2f",species[sp].Data(),s_pid[ipid].Data(),prob_cut);
    eff_den[ipid] = new TH2F(Form("eff_den_%s",s_pid[ipid].Data()),title+";mom[GeV/c];#theta[rad]",nbin,0,mom_max,nbin,0,the_max);
    eff_num[ipid] = new TH2F(Form("eff_num_%s",s_pid[ipid].Data()),title+";mom[GeV/c];#theta[rad]",nbin,0,mom_max,nbin,0,the_max);
    eff1d_mom[ipid] = new TEfficiency(Form("eff1d_mom_%s",s_pid[ipid].Data()), title+";mom[GeV/c]", nbin, 0, mom_max);
    eff1d_the[ipid] = new TEfficiency(Form("eff1d_the_%s",s_pid[ipid].Data()), title+";#theta[rad]", nbin, 0, the_max);
    eff2d[ipid] = new TEfficiency(Form("eff2d_%s",s_pid[ipid].Data()), title+";mom[GeV/c];#theta[rad]", nbin, 0, mom_max, nbin, 0, the_max);
  }

  for (int ievt = 0; ievt < nevt; ++ievt) {
    if (ievt%1000==0) cout << "ievt= " << ievt << " / " << nevt <<endl;
    //if (ievt>10) break;
    cbmsim->GetEntry(ievt);
    int ntrk = cCand_array->GetEntriesFast();
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

      // combined probability of this track being a given type
      double prob_comb[npid_max] = {0.0};
      double prob_comb[iel] = get_comb_prob_elec(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc);
      double prob_comb[imu] = get_comb_prob_muon(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc);
      double prob_comb[ipi] = get_comb_prob_pion(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc);
      double prob_comb[ik]  = get_comb_prob_kaon(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc);
      double prob_comb[iprot] = get_comb_prob_proton(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc);
      double prob_comb[iposit]  = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetElectronPidProb);

      for (int ipid=0; ipid<npid_max; ++ipid) {
      	eff_den[ipid]->Fill(mom,the);
      	if (prob_comb[ipid]>prob_cut) { eff_num[ipid]->Fill(mom, the); }
      	eff1d_the[ipid]->Fill(prob_comb[ipid]>prob_cut,the);
      	eff1d_mom[ipid]->Fill(prob_comb[ipid]>prob_cut,mom);
      	eff2d[ipid]->Fill(prob_comb[ipid]>prob_cut,mom,the);
      }

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

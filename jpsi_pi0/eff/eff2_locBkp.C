//typedef double(PndPidProbability::*prob_func)(PndPidProbability*) const;
typedef PndPidProbability* ppp_p;

//double get_comb_prob(ppp_p prob_emc, ppp_p prob_stt, ppp_p prob_mvd, ppp_p prob_drc, ppp_p prob_disc, prob_func func) {
double get_comb_prob(ppp_p prob_emc, ppp_p prob_stt, ppp_p prob_mvd, ppp_p prob_drc, ppp_p prob_disc, double(PndPidProbability::*prob_func)(PndPidProbability*) func) {
  Double_t prob_emc = (prob_emc->*func)();
  Double_t prob_stt = (prob_stt->*func)();
  Double_t prob_mvd = (prob_mvd->*func)();
  Double_t prob_drc = (prob_drc->*func)();
  Double_t prob_disc = (prob_disc->*func)();
  Double_t xx = (prob_drc/(1-prob_drc))*(prob_disc/(1-prob_disc))
    *(prob_mvd/(1-prob_mvd))*(prob_stt/(1-prob_stt))
    *(prob_emc/(1-prob_emc));
  return xx/(xx+1);
}

void eff2(int sp /*species index*/, int ifile ) {

  const double mom_max = 2.0;
  const double the_max = 180.0;

  const TString basedir="/Users/tujuba/panda/work/jpsi_pi0/grid.out/jacek/";
  const TString subdir = Form("/runall.%d/",ifile);

  //const int nsp_max = 10;
  const enum{ielec=0, imuonm, ipionm, ikaonm, iantiprot, iposit, imuonp, ipionp, ikaonp, iprot, nsp_max};
  const TString species[nsp_max] = {"posit","muplus","piplus","kplus","prot","elec","muminus","piminus","kminus","antiprot"};

  //const int ndet = 5;
  const enum{iemc = 0, istt, imvd, idirc, idisc, ndet};
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
  TEfficiency *eff2d[nsp_max], *eff1d_the[nsp_max], eff1d_mom[nsp_max];
  for (int isp = 0; isp < nsp_max; ++isp) {

    eff_den[isp] = new TH2F(Form("eff_den_%s",species[isp].Data()),Form("eff_den_%s",species[isp].Data()),200,0,mom_max,200,0,the_max);
    eff_num[isp] = new TH2F(Form("eff_den_%s",species[isp].Data()),Form("eff_den_%s",species[isp].Data()),200,0,mom_max,200,0,the_max);
    eff1d_mom[isp] = new TEfficiency(Form("eff1d_the_%s",species[isp].Data()), Form("eff1d_the_%s",species[isp].Data()), 200, 0, the_max);
    eff1d_the[isp] = new TEfficiency(Form("eff1d_mom_%s",species[isp].Data()), Form("eff1d_mom_%s",species[isp].Data()), 200, 0, the_max);
    eff2d[isp] = new TEfficiency(Form("eff2d_%s",species[isp].Data()), Form("eff2d_%s",species[isp].Data()), 200, 0, mom_max, 200, 0, the_max);
  }

  //  for (int ievt = 0; ievt < nevt; ++ievt) {
  for (int ievt = 0; ievt < 1000; ++ievt) {
    cout << "ievt= " << ievt << " / " << nevt <<endl;

    int ntrk = cbmsim->GetEntriesFast();
    for (int itrk = 0; itrk < ntrk; ++itrk) {

      PndPidCandidate *cand = (PndPidCandidate*) cCand_array->At(itrk);

      Double_t mom = cand->GetMomentum().Mag();
      Double_t the = TMath::RadToDeg()*cand->GetMomentum().Mag();
      //Double_t eoverp = cand->GetEmcCalEnergy()/mom;
      //Double_t stt_dedx = cand->GetSttMeanDEDX();
      //Double_t disc_thetaC = TMath::RadToDeg()*cand->GetDiscThetaC();
      //Double_t drc_thetaC = TMath::RadToDeg()*cand->GetDrcThetaC();
      //Double_t muo_iron = cand->GetMuoIron();
      //Double_t mvd_dedx = 1000*cand->GetMvdDEDX();

      ppp_p prob_drc = (ppp_p) drc_array->At(itrk);
      ppp_p prob_disc =(ppp_p) disc_array->At(itrk);
      ppp_p prob_mvd = (ppp_p) mvd_array->At(itrk);
      ppp_p prob_stt = (ppp_p) stt_array->At(itrk);
      ppp_p prob_emc = (ppp_p) emcb_array->At(itrk);

      // combined probability of this track being a
      double prob_comb[nsp_max] = {0.0};
      double prob_comb[ielec] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetElectronPidProb);
      double prob_comb[imuonm] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetMuonPidProb);
      double prob_comb[ipionm] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetPionPidProb);
      double prob_comb[ikaonm] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetKaonPidProb);
      double prob_comb[iantiprot] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetProtonPidProb);
      double prob_comb[iposit] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetElectronPidProb);
      double prob_comb[imuonp] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetMuonPidProb);
      double prob_comb[ipionp] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetPionPidProb);
      double prob_comb[ikaonp] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetKaonPidProb);
      double prob_comb[iprot] = get_comb_prob(prob_emc, prob_stt, prob_mvd, prob_drc, prob_disc, &PndPidProbability::GetProtonPidProb);

      for (int isp=0; isp<nsp_max; ++isp) {
	eff_den[isp]->Fill(mom,the);
	if (prob_comb[isp]>0.99) { eff_num[isp]->Fill(mom, the); }
	eff1d_the[isp]->Fill(the, prob_comb[isp]>0.99);
	eff1d_mom[isp]->Fill(the, prob_comb[isp]>0.99);
	eff1d_mom[isp]->Fill(the, prob_comb[isp]>0.99);
      }

    }

  }

}

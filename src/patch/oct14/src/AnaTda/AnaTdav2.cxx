#include "AnaTdav2.h"

#include "FairTask.h"
#include "RhoCandList.h"
#include "PndKinFitter.h"
#include "PndPidProbability.h"
#include "FairRootManager.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

using namespace std;

AnaTdav2::AnaTdav2(const int& _iplab, const int& itype, const int& brem, const int& eid_meth):
  verb(false),
  nevt(0),
  brem_corr(brem!=0),
  mc_type(itype),
  eid_param_method(eid_meth!=0),
  iplab((_iplab>=0&&_iplab<3)?_iplab:0),
  plab{5.513,8.,12.},
  mom_antip(plab[iplab]),
  boost_to_cm(),
  boost_to_lab(),
  p4pbar(),
  p4targ(),
  p4sys(),
  event_t(-9999.),
  event_u(-9999.),
  //tmin{-0.443789, -0.5, -0.5},
  //tmax{0.616486, 0.457248, 0.31538},
  //tmin{-0.443789, -2.76, -6.50},
  //tmax{0.616486, 0.457248, 0.31538},
  //tmin{-0.092, -1.3, -2.85},
  tmin{-0.092, -1.0, -1.0},
  tmid{-0.45, -2.76, -6.5},
  tmax{0.59, 0.43, 0.3},
  //nevt_sim{81874.0, 224120.0, 189015.0},
  // 1st index: 0->pi0pipm, 1->pi0jpsi, 2->pi02pipm, 3->pi0pipm2, 4->pi02jpsi
  nevt_sim{
    {814794.0,888292.0,898721.0},
      {32780.0,50142.0,51860.0},
	{214780.0,174864.0,160099.0},
	  {570751.0,609044.0,527506.0},
	    {200000.0,200000.0,200000.0},
	      {100000.0,100000.0,100000.0},
		{100000.0,100000.0,100000.0}},
  //nevt_xsect{4.0e11, 1e11, 1e9},
  // xsect={0.2mb, 0.05mb, 0.02mb}
  nevt_xsect{
    {4.0e11, 1e11, 2e10},
      {32780.0,50142.0,51860.0},
	{1.15e12, 3.15e11, 6.84e10},
	  {3.19e12, 1.14e12, 2.92e11},
	    {94243.0, 157947.3, 177361.2},
	      {100000.0,100000.0,100000.0},
		{100000.0,100000.0,100000.0}},
  eff_file_name("eff/effic_smooth.root"),
  eff_hist_name("eff_ep_em_rad"),
  eff_hist_rad(true),
  pi_eff_file_name("eff/hadd_out/eff.pi.root"),
  pi_eff_hist_name("prob_cut_5/eff2d_e_id"),
  pi_eff_hist_rad(false),
  eid_prob_min(0.5),
  mcList(),
  apply_pi0evsoa_cut(true),
  lw{{0.11,-0.05,+0.00},{0.12,-0.03,-0.03},{0.15,-0.06,-0.07}},
  up{{0.14,0.21,0.07},{0.12,0.15,0.15},{0.11,0.10,0.20}},
  apply_pi0m_cut(true),
  pi0m_cut_min(0.11),
  pi0m_cut_max(0.16),
  apply_mtot_cut(true),
  mtot_cut_min(3.0),
  mtot_cut_max(3.7),
  apply_dth_dph_cut(true),
  dth_sigma(0.4),
  dph_sigma(0.4),
  dth_dph_cm_cut_max(3.0),
  require_exclusivity(false),
  jpsi_m_3sig_min(2.8),
  jpsi_m_3sig_max(3.3),
  chi2_cut{20.0,50.0,100.0}
{
  assert(iplab>=0&&iplab<=2);
  rcl.resize(nrcl);
  gRandom->SetSeed();
  // bootstrap the pi0pi0jpsi (idx4) x-sect from pi0jpsi (idx1) xsect using the pi0pi0pipipm(idx2)/pi0pippim(idx0) ratio
  //nevt_xsect[4][0] = 60000.;
  //nevt_xsect[4][1] = 60000.*(nevt_xsect[1][1]/nevt_xsect[1][0]);
  //nevt_xsect[4][2] = 60000.*(nevt_xsect[1][2]/nevt_xsect[1][0]);
  for (int ii=0; ii < 3; ++ii) {
    nevt_xsect[4][ii] = nevt_xsect[1][ii]*(nevt_xsect[2][ii]/nevt_xsect[0][ii]);
    //nevt_xsect[4][ii] = nevt_xsect[1][ii]*(nevt_xsect[2][ii]/nevt_xsect[0][ii]);
    cout << " pi0pipm nevt_xsect of iplab("<< iplab<< ") = " << nevt_xsect[0][ii] << endl;
    cout << " pi0jpsi nevt_xsect of iplab("<< iplab<< ") = " << nevt_xsect[1][ii] << endl;
    cout << " pi02pipm nevt_xsect of iplab("<< iplab<< ") = " << nevt_xsect[2][ii] << endl;
    cout << " pi0pipm2 nevt_xsect of iplab("<< iplab<< ") = " << nevt_xsect[3][ii] << endl;
    cout << " pi0pi0jpsi nevt_xsect of iplab("<< iplab<< ") = " << nevt_xsect[4][ii] << endl;
  }
}

AnaTdav2::~AnaTdav2() {

}

TClonesArray* AnaTdav2::init_tca(TString name) {
  TClonesArray *tca = dynamic_cast<TClonesArray *> (m_ioman->GetObject(name));
  if ( ! tca ) {
    cout << "-W- EffHists::Init: "  << "No " << name << " array!" << endl;
    return NULL;
  } else {
    cout << "-I- EffHists::Init: "  << " finished reading " << name << " array!" << endl;
    return tca;
  }
}

void AnaTdav2::init_tcas() {
  m_ioman = FairRootManager::Instance();
  if ( ! m_ioman ){
    cout << "-E- EffHists::Init: "
	 << "RootManager not instantiated!" << endl;
    return;
  }
  m_cand_array = init_tca( "PidChargedCand");
  m_drc_array = init_tca( "PidAlgoDrc");
  m_disc_array = init_tca( "PidAlgoDisc");
  m_stt_array = init_tca( "PidAlgoStt");
  m_mvd_array = init_tca( "PidAlgoMvd");
  m_emcb_array = init_tca( "PidAlgoEmcBayes");
}

double _pi_eff_func(double *x, double *p) {
  double xx=x[0];
  double f1 = p[0]+  p[2]*TMath::Sin(xx*p[1]*1.0)+   p[3]*TMath::Sin(xx*p[1]*2.0)+  p[4]*TMath::Cos(xx*p[1]*1.0) + p[5]*TMath::Cos(xx*p[1]*2.0);
  double t1 = f1/(p[6]+(p[7]*TMath::Power(xx,p[8])));
  double p2 = p[9]+  p[10]*x[0]+  p[11]*x[0]*x[0];
  double t2 = p[12]*p2;
  return p[13]+t1+t2;
}

void AnaTdav2::init_hists() {
  eff_file = TFile::Open(eff_file_name.c_str());

  heff_epm = (TH2F*) eff_file->Get(eff_hist_name.c_str())->Clone("heff_epm");

  pi_eff_file = TFile::Open(pi_eff_file_name.c_str());

  TEfficiency *tmp = (TEfficiency*) pi_eff_file->Get(pi_eff_hist_name.c_str());
  if (tmp->GetDimension()==2) {
    //pi_eff = smooth_eff2d((TEfficiency*)pi_eff_file->Get(pi_eff_hist_name.c_str())->Clone("pi_eff"),1000);
    //pi_eff = smooth_eff2d(rebin2d((TEfficiency*)pi_eff_file->Get(pi_eff_hist_name.c_str())->Clone("pi_eff"), 5),1000);
    pi_eff = (TEfficiency*) pi_eff_file->Get(pi_eff_hist_name.c_str())->Clone("pi_eff2d");
  } else {
    //pi_eff = (TEfficiency*) pi_eff_file->Get(pi_eff_hist_name.c_str())->Clone("pi_eff");
    pi_eff = smooth_eff1d((TEfficiency*) pi_eff_file->Get(pi_eff_hist_name.c_str())->Clone("pi_eff1d"));
    double pars[14] = { 6.35356e+00, 1.0, 4.13113e+00, -4.43669e+00,
			0.1, 0.01, 0.1, 9.64513e+03, 1.22279e+00,
			4.66147e-04, 2.96494e-05, -6.21090e-06, 1.0, -3.23049e-06 };
    pi_eff_func = new TF1("pi_eff_func",_pi_eff_func,0.0001,10,14);
    for (int ii=0; ii < 14; ++ii) {
      if (ii==9||ii==10||ii==11)
	pi_eff_func->FixParameter(ii,pars[ii]);
      else
	pi_eff_func->SetParameter(ii,pars[ii]);
    }
    pi_eff->Fit(pi_eff_func,"+RME");
  }
  hwt = new TH1F("hwt", "event by event weight distribution", 2000,0,100);

  // These do not depend on steps
  hng = new TH1F(Form("hng"),Form("hng"),11,-0.5,10.5);
  hng20mev = new TH1F(Form("hng20mev"),Form("hng20mev"),11,-0.5,10.5);
  hnch = new TH1F(Form("hnch"),Form("hnch"),11,-0.5,10.5);

  for (int is = 0; is< nstep; ++is) {
    hmep[is] = new TH1F(Form("hmep_%d",is),Form("hmep_%d",is),400,0,5);
    hmep_valid[is] = new TH1F(Form("hmep_valid_%d",is),Form("hmep_valid_%d",is),400,0,5);
    hmgg[is] = new TH1F(Form("hmgg_%d",is),Form("hmgg_%d",is),400,0,0.2);
    hmgg_valid[is] = new TH1F(Form("hmgg_valid_%d",is),Form("hmgg_valid_%d",is),400,0,0.2);
    hnep[is] = new TH1F(Form("hnep_%d",is),Form("hnep_%d",is),6,-0.5,5.5);
    hngg[is] = new TH1F(Form("hngg_%d",is),Form("hngg_%d",is),26,-0.5,25.5);
    hnpi0jpsi[is] = new TH1F(Form("hnpi0jpsi_%d",is),Form("hnpi0jpsi_%d",is),6,-0.5,5.5);

    hmep_mct[is] = new TH1F(Form("hmep_mct_%d",is),Form("hmep_mct_%d",is),400,0,5);
    hmep_non_mct[is] = new TH1F(Form("hmep_non_mct_%d",is),Form("hmep_non_mct_%d",is),400,0,5);
    hthe_ep_mct[is] = new TH1F(Form("hthe_ep_mct_%d",is),Form("hthe_ep_mct_%d",is),400,0,TMath::Pi());
    hthe_ep_mct_fwd[is] = new TH1F(Form("hthe_ep_mct_fwd_%d",is),Form("hthe_ep_mct_fwd_%d",is),400,0,TMath::Pi());
    hthe_ep_mct_bwd[is] = new TH1F(Form("hthe_ep_mct_bwd_%d",is),Form("hthe_ep_mct_bwd_%d",is),400,0,TMath::Pi());
    hmgg_mct[is] = new TH1F(Form("hmgg_mct_%d",is),Form("hmgg_mct_%d",is),400,0.0,0.2);
    hmgg_non_mct[is] = new TH1F(Form("hmgg_non_mct_%d",is),Form("hmgg_non_mct_%d",is),400,0.0,0.2);
    hthe_gg_mct[is] = new TH1F(Form("hthe_gg_mct_%d",is),Form("hthe_gg_mct_%d",is),400,0,TMath::Pi());
    hthe_gg_mct_fwd[is] = new TH1F(Form("hthe_gg_mct_fwd_%d",is),Form("hthe_gg_mct_fwd_%d",is),400,0,TMath::Pi());
    hthe_gg_mct_bwd[is] = new TH1F(Form("hthe_gg_mct_bwd_%d",is),Form("hthe_gg_mct_bwd_%d",is),400,0,TMath::Pi());
    hoa_gg_mct[is] = new TH1F(Form("hoa_gg_mct_%d",is),Form("hoa_gg_mct_%d",is),400,0,2*TMath::Pi());
  }

  hpi0th = new TH1F("hpi0th", "hpi0th", 1000, 0, TMath::Pi());
  hpi0cost_cm = new TH1F("hpi0cost_cm", "hpi0cost_cm", 1100, -1.1, 1.1);
  hpi0th_mcut = new TH1F("hpi0th_mcut", "hpi0th_mcut", 1000, 0, TMath::Pi());
  hpi0cost_cm_mcut = new TH1F("hpi0cost_cm_mcut", "hpi0cost_cm_mcut", 1100, -1.1, 1.1);
  for (int ib=0; ib<tu_binning.size()-1; ++ib) {
    hmep_pi0cost_cm.push_back( new TH1F(Form("hmep_pi0cost_cm_%d", ib), Form("hmep_pi0cost_cm_%d", ib), 200, 0, 5));
    hmep_pi0th.push_back( new TH1F(Form("hmep_pi0th_%d", ib), Form("hmep_pi0th_%d", ib), 200, 0, 5));
    hmept.push_back(new TH1F(Form("hmept%d", ib), Form("%4.2f < t < %4.2f;M_{inv}", tu_binning[ib], tu_binning[ib+1]), 200, 0, 5));
    hmepu.push_back(new TH1F(Form("hmepu%d", ib), Form("%4.2f < u < %4.2f;M_{inv}", tu_binning[ib], tu_binning[ib+1]), 200, 0, 5));
  }

  for (int iby=0; iby < costh_binning_2d.size()-1; ++iby) {
    for (int ibx=0; ibx < tu_binning_2d.size()-1; ++ibx) {

      hmeptcth.push_back(new TH1F(Form("hmep_t%d_cth%d", ibx, iby), Form("%5.3f < t < %5.3f & %4.2f < cos(#theta) < %4.2f;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      hmeptcth[hmeptcth.size()-1]->Sumw2();
      hmeptcth0.push_back(new TH1F(Form("hmep_t%d_cth%d_wt0", ibx, iby), Form("%5.3f < t < %5.3f & %4.2f < cos(#theta) < %4.2f, A=1.0;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      hmeptcth0[hmeptcth0.size()-1]->Sumw2();
      hmeptcth1.push_back(new TH1F(Form("hmep_t%d_cth%d_wt1", ibx, iby), Form("%5.3f < t < %5.3f & %4.2f < cos(#theta) < %4.2f, A=0.4;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      hmeptcth1[hmeptcth1.size()-1]->Sumw2();

      f_hmeptcth.push_back(new TH1F(Form("f_hmep_t%d_cth%d", ibx, iby), Form("%5.3f < t < %5.3f & %4.2f < cos(#theta) < %4.2f;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      f_hmeptcth[f_hmeptcth.size()-1]->Sumw2();
      f_hmeptcth0.push_back(new TH1F(Form("f_hmep_t%d_cth%d_wt0", ibx, iby), Form("%5.3f < t < %5.3f & %4.2f < cos(#theta) < %4.2f, A=1.0;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      f_hmeptcth0[f_hmeptcth0.size()-1]->Sumw2();
      f_hmeptcth1.push_back(new TH1F(Form("f_hmep_t%d_cth%d_wt1", ibx, iby), Form("%5.3f < t < %5.3f & %4.2f < cos(#theta) < %4.2f, A=0.4;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      f_hmeptcth1[f_hmeptcth1.size()-1]->Sumw2();

    }
  }
  for (int iby=0; iby < costh_binning_2d.size()-1; ++iby) {
    for (int ibx=0; ibx < tu_binning_2d.size()-1; ++ibx) {

      hmepucth.push_back(new TH1F(Form("hmep_u%d_cth%d", ibx, iby), Form("%5.3f < u < %5.3f & %4.2f < cos(#theta) < %4.2f;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      hmepucth[hmepucth.size()-1]->Sumw2();
      hmepucth0.push_back(new TH1F(Form("hmep_u%d_cth%d_wt0", ibx, iby), Form("%5.3f < u < %5.3f & %4.2f < cos(#theta) < %4.2f, A=1.0;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      hmepucth0[hmepucth0.size()-1]->Sumw2();
      hmepucth1.push_back(new TH1F(Form("hmep_u%d_cth%d_wt1", ibx, iby), Form("%5.3f < u < %5.3f & %4.2f < cos(#theta) < %4.2f, A=0.4;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      hmepucth1[hmepucth1.size()-1]->Sumw2();

      f_hmepucth.push_back(new TH1F(Form("f_hmep_u%d_cth%d", ibx, iby), Form("%5.3f < u < %5.3f & %4.2f < cos(#theta) < %4.2f;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      f_hmepucth[f_hmepucth.size()-1]->Sumw2();
      f_hmepucth0.push_back(new TH1F(Form("f_hmep_u%d_cth%d_wt0", ibx, iby), Form("%5.3f < u < %5.3f & %4.2f < cos(#theta) < %4.2f, A=1.0;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      f_hmepucth0[f_hmepucth0.size()-1]->Sumw2();
      f_hmepucth1.push_back(new TH1F(Form("f_hmep_u%d_cth%d_wt1", ibx, iby), Form("%5.3f < u < %5.3f & %4.2f < cos(#theta) < %4.2f, A=0.4;M_{inv}", tu_binning_2d[ibx], tu_binning_2d[ibx+1], costh_binning_2d[iby], costh_binning_2d[iby+1]), 200, 0, 5));
      f_hmepucth1[f_hmepucth1.size()-1]->Sumw2();

    }
  }

  hmmiss = new TH1F("hmmiss", "hmmiss", 400, -0.6, 0.6);
  hmmiss2 = new TH1F("hmmiss2", "hmmiss2", 400, -0.2, 0.3);
  hmmiss_jpsi = new TH1F("hmmiss_jpsi", "hmmiss_jpsi", 400, -0.6, 0.9);
  hmmiss2_jpsi = new TH1F("hmmiss2_jpsi", "hmmiss2_jpsi", 400, -0.2, 0.5);

  hmtot = new TH1F("hmtot", "hmtot", 200, 0, 8);
  hcmoa = new TH2F("hcmoa", "hcmoa", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());

  hmep_mconst = new TH1F(Form("hmep_mconst"),Form("hmep_mconst"),200,0,5);
  hmtot_mconst = new TH1F("hmtot_mconst", "hmtot_mconst", 200, 0, 8);
  hcmoa_mconst = new TH2F("hcmoa_mconst", "hcmoa_mconst", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());
  hmtot_mconst_cut = new TH1F("hmtot_mconst_cut", "hmtot_mconst_cut", 200, 0, 8);
  hcmoa_mconst_cut = new TH2F("hcmoa_mconst_cut", "hcmoa_mconst_cut", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());

  double _tmin = iplab==0?-1.9:(iplab==1?-6.5:-14);
  if (is_dpm()) {
    _tmin = iplab==0?-11:(iplab==1?-16:-24);
  }

  double tbins_da[40] = {0.};
  for (int ib=0; ib<tu_binning.size(); ++ib) { tbins_da[ib] = tu_binning[ib]; }
  httrumc_vb = new TH1F("httrumc_vb", "httrumc_vb", tu_binning.size()-1, tbins_da);
  hutrumc_vb = new TH1F("hutrumc_vb", "hutrumc_vb", tu_binning.size()-1, tbins_da);

  httrumc = new TH1F("httrumc", "httrumc", 1000, _tmin, 1);
  hutrumc = new TH1F("hutrumc", "hutrumc", 1000, _tmin, 1);
  htrecgg = new TH1F("htrecgg", "htrecgg", 1000, _tmin, 1);
  hurecgg = new TH1F("hurecgg", "hurecgg", 1000, _tmin, 1);
  htrecep = new TH1F("htrecep", "htrecep", 1000, _tmin, 1);
  hurecep = new TH1F("hurecep", "hurecep", 1000, _tmin, 1);
  htrecgg_mcut = new TH1F("htrecgg_mcut", "htrecgg_mcut", 1000, _tmin, 1);
  hurecgg_mcut = new TH1F("hurecgg_mcut", "hurecgg_mcut", 1000, _tmin, 1);
  htrecep_mcut = new TH1F("htrecep_mcut", "htrecep_mcut", 1000, _tmin, 1);
  hurecep_mcut = new TH1F("hurecep_mcut", "hurecep_mcut", 1000, _tmin, 1);
  htresgg = new TH1F("htresgg", "htresgg", 1000, -3, 3);
  huresgg = new TH1F("huresgg", "huresgg", 1000, -3, 3);
  htresep = new TH1F("htresep", "htresep", 1000, -3, 3);
  huresep = new TH1F("huresep", "huresep", 1000, -3, 3);

  hepcosth_res = new TH1F("hepcosth_res", "hepcosth_res", 1000, -0.5, 0.5);
  f_hepcosth_res = new TH1F("f_hepcosth_res", "f_hepcosth_res", 1000, -0.5, 0.5);

  hepcosth_jpsi_mc_all = new TH1F("hepcosth_jpsi_mc_all", "hepcosth_jpsi_mc_all", 1000, -1.1, 1.1);
  hepcosth_jpsi_mc_all_wt0 = new TH1F("hepcosth_jpsi_mc_all_wt0", "hepcosth_jpsi_mc_all_wt0", 1000, -1.1, 1.1);
  hepcosth_jpsi_mc_all_wt1 = new TH1F("hepcosth_jpsi_mc_all_wt1", "hepcosth_jpsi_mc_all_wt1", 1000, -1.1, 1.1);
  hepcosth_jpsi_vs_epthlab_mc_all = new TH2F("hepcosth_jpsi_vs_epthlab_mc_all", "hepcosth_jpsi_vs_epthlab_mc_all", 200, -1.1, 1.1, 200, 0, 180);
  hepcosth_jpsi_vs_emthlab_mc_all = new TH2F("hepcosth_jpsi_vs_emthlab_mc_all", "hepcosth_jpsi_vs_emthlab_mc_all", 200, -1.1, 1.1, 200, 0, 180);

  hepcosth_jpsi_rec_all = new TH1F("hepcosth_jpsi_rec_all", "hepcosth_jpsi_rec_all", 1000, -1.1, 1.1);
  hepcosth_jpsi_vs_epthlab_rec_all = new TH2F("hepcosth_jpsi_vs_epthlab_rec_all", "hepcosth_jpsi_vs_epthlab_rec_all", 200, -1.1, 1.1, 200, 0, 180);
  hepcosth_jpsi_vs_emthlab_rec_all = new TH2F("hepcosth_jpsi_vs_emthlab_rec_all", "hepcosth_jpsi_vs_emthlab_rec_all", 200, -1.1, 1.1, 200, 0, 180);
  f_hepcosth_jpsi_rec_all = new TH1F("f_hepcosth_jpsi_rec_all", "f_hepcosth_jpsi_rec_all", 1000, -1.1, 1.1);
  f_hepcosth_jpsi_vs_epthlab_rec_all = new TH2F("f_hepcosth_jpsi_vs_epthlab_rec_all", "f_hepcosth_jpsi_vs_epthlab_rec_all", 200, -1.1, 1.1, 200, 0, 180);
  f_hepcosth_jpsi_vs_emthlab_rec_all = new TH2F("f_hepcosth_jpsi_vs_emthlab_rec_all", "f_hepcosth_jpsi_vs_emthlab_rec_all", 200, -1.1, 1.1, 200, 0, 180);

  for (int ii=0; ii < 4; ++ii) {
    hepcosth_jpsi_rec[ii] = new TH1F(Form("hepcosth_jpsi_rec_itu%d",ii), Form("hepcosth_jpsi_rec_itu%d",ii), 1000, -1.1, 1.1);
    hepcosth_jpsi_vs_epthlab_rec[ii] = new TH2F(Form("hepcosth_jpsi_vs_epthlab_rec_itu%d",ii), Form("hepcosth_jpsi_vs_epthlab_rec_itu%d",ii), 1000, -1.1, 1.1, 200, 0, 180);
    hepcosth_jpsi_vs_emthlab_rec[ii] = new TH2F(Form("hepcosth_jpsi_vs_emthlab_rec_itu%d",ii), Form("hepcosth_jpsi_vs_emthlab_rec_itu%d",ii), 1000, -1.1, 1.1, 200, 0, 180);

    f_hepcosth_jpsi_rec[ii] = new TH1F(Form("f_hepcosth_jpsi_rec_itu%d",ii), Form("hepcosth_jpsi_rec_itu%d",ii), 1000, -1.1, 1.1);
    f_hepcosth_jpsi_vs_epthlab_rec[ii] = new TH2F(Form("f_hepcosth_jpsi_vs_epthlab_rec_itu%d",ii), Form("hepcosth_jpsi_vs_epthlab_rec_itu%d",ii), 1000, -1.1, 1.1, 200, 0, 180);
    f_hepcosth_jpsi_vs_emthlab_rec[ii] = new TH2F(Form("f_hepcosth_jpsi_vs_emthlab_rec_itu%d",ii), Form("hepcosth_jpsi_vs_emthlab_rec_itu%d",ii), 1000, -1.1, 1.1, 200, 0, 180);

    hepcosth_jpsi_mc[ii] = new TH1F(Form("hepcosth_jpsi_mc_itu%d",ii), Form("hepcosth_jpsi_mc_itu%d",ii), 1000, -1.1, 1.1);
    hepcosth_jpsi_vs_epthlab_mc[ii] = new TH2F(Form("hepcosth_jpsi_vs_epthlab_mc_itu%d",ii), Form("hepcosth_jpsi_vs_epthlab_mc_itu%d",ii), 1000, -1.1, 1.1, 200, 0, 180);
    hepcosth_jpsi_vs_emthlab_mc[ii] = new TH2F(Form("hepcosth_jpsi_vs_emthlab_mc_itu%d",ii), Form("hepcosth_jpsi_vs_emthlab_mc_itu%d",ii), 1000, -1.1, 1.1, 200, 0, 180);
  }

  htrupi0thcm = new TH1F("htrupi0thcm", "htrupi0thch", 1000, 0., TMath::Pi());
  htrupi0costhcm = new TH1F("htrupi0costhcm", "htrupi0costhcm", 1100, -1.1, 1.1);
  htrupi0thlab = new TH1F("htrupi0thlab", "htrupi0thlab", 1000, 0., TMath::Pi());
  htrupi0thcm_mcut = new TH1F("htrupi0thcm_mcut", "htrupi0thch_mc", 1000, 0., TMath::Pi());
  htrupi0costhcm_mcut = new TH1F("htrupi0costhcm_mcut", "htrupi0costhcm_mcut", 1100, -1.1, 1.1);
  htrupi0thlab_mcut = new TH1F("htrupi0thlab_mcut", "htrupi0thlab_mcut", 1000, 0., TMath::Pi());
  htrupi0thcm_tcut = new TH1F("htrupi0thcm_tcut", "htrupi0thch_tc", 1000, 0., TMath::Pi());
  htrupi0costhcm_tcut = new TH1F("htrupi0costhcm_tcut", "htrupi0costhcm_tcut", 1100, -1.1, 1.1);
  htrupi0thlab_tcut = new TH1F("htrupi0thlab_tcut", "htrupi0thlab_tcut", 1000, 0., TMath::Pi());
  htrupi0thcm_tcut_mcut = new TH1F("htrupi0thcm_tcut_mcut", "htrupi0thch_tc_mc", 1000, 0., TMath::Pi());
  htrupi0costhcm_tcut_mcut = new TH1F("htrupi0costhcm_tcut_mcut", "htrupi0costhcm_tcut_mcut", 1100, -1.1, 1.1);
  htrupi0thlab_tcut_mcut = new TH1F("htrupi0thlab_tcut_mcut", "htrupi0thlab_tcut_mcut", 1000, 0., TMath::Pi());

  htrupi0thcm_vs_m = new TH2F("htrupi0thcm_vs_m", "htrupi0thch_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());
  htrupi0costhcm_vs_m = new TH2F("htrupi0costhcm_vs_m", "htrupi0costhcm_vs_m", 200, 0, 5, 1100, -1.1, 1.1);
  htrupi0thlab_vs_m = new TH2F("htrupi0thlab_vs_m", "htrupi0thlab_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());
  htrupi0thcm_mcut_vs_m = new TH2F("htrupi0thcm_mcut_vs_m", "htrupi0thch_mc_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());
  htrupi0costhcm_mcut_vs_m = new TH2F("htrupi0costhcm_mcut_vs_m", "htrupi0costhcm_mcut_vs_m", 200, 0, 5, 1100, -1.1, 1.1);
  htrupi0thlab_mcut_vs_m = new TH2F("htrupi0thlab_mcut_vs_m", "htrupi0thlab_mcut_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());

  hnevt =  new TH1F("hnevt","hnevt", 10,0,10);
  hnevt->SetBinContent(3, nevt_sim[mc_type][iplab]);
  hnevt->SetBinContent(4, nevt_xsect[mc_type][iplab]);

  hpi0jpsi_chi24c = new TH1F("hpi0jpsi_chi24c","hpi0jpsi_chi24c",2000,0,10000);
  hpi0jpsi_chi24c_c = new TH1F("hpi0jpsi_chi24c_c","hpi0jpsi_chi24c_c",1000,0,500);
  hpi0jpsi_prob4c  = new TH1F("hpi0jpsi_prob4c ","hpi0jpsi_prob4c ",1000,0,1.0);
  hpi0jpsi_pull4c  = new TH1F("hpi0jpsi_pull4c ","hpi0jpsi_pull4c ",1000,-6,6);
  hpi0jpsi_chi2diff4c  = new TH1F("hpi0jpsi_chi2diff4c ","hpi0jpsi_chi2diff4c ",1000,-100,100);

  double hist_mmin[3] = {2.6,2.8,3.8};
  double hist_mmax[3] = {3.8,5.,6.};
  hpi0jpsi_chi24c_vs_mtot_r = new TH2F("hpi0jpsi_chi24c_vs_mtot_r","hpi0jpsi_chi24c_vs_mtot_r",2000,0,10000,200,hist_mmin[iplab],hist_mmax[iplab]);
  hpi0jpsi_chi24c_vs_cm_dth_r = new TH2F("hpi0jpsi_chi24c_vs_cm_dth_r","hpi0jpsi_chi24c_vs_cm_dth_r",2000,0,10000,200,1.0,5.0);
  hpi0jpsi_chi24c_vs_cm_dph_r = new TH2F("hpi0jpsi_chi24c_vs_cm_dph_r","hpi0jpsi_chi24c_vs_cm_dph_r",2000,0,10000,200,1.5,5.0);
  hpi0jpsi_chi24c_vs_mtot_f = new TH2F("hpi0jpsi_chi24c_vs_mtot_f","hpi0jpsi_chi24c_vs_mtot_f",500,0,1000,200,hist_mmin[iplab],hist_mmax[iplab]);
  hpi0jpsi_chi24c_vs_cm_dth_f = new TH2F("hpi0jpsi_chi24c_vs_cm_dth_f","hpi0jpsi_chi24c_vs_cm_dth_f",500,0,1000,500,1.0,5.0);
  hpi0jpsi_chi24c_vs_cm_dph_f = new TH2F("hpi0jpsi_chi24c_vs_cm_dph_f","hpi0jpsi_chi24c_vs_cm_dph_f",500,0,1000,500,1.5,5.0);

  hpi0pi0jpsi_chi24c = new TH1F("hpi0pi0jpsi_chi24c","hpi0pi0jpsi_chi24c",2000,0,10000);
  hpi0pi0jpsi_chi24c_c = new TH1F("hpi0pi0jpsi_chi24c_c","hpi0pi0jpsi_chi24c_c",1000,0,500);
  hpi0vs2pi0_chi24c = new TH2F("hpi0vs2pi0_chi24c","hpi0vs2pi0_chi24c",2000,0,10000,2000,0,10000);
  hpi0vs2pi0_chi24c_c = new TH2F("hpi0vs2pi0_chi24c_c","hpi0vs2pi0_chi24c_c",1000,0,500,1000,0,500);

  // pull distributions
  double zoom = (mc_type == 1||mc_type == 5)? 1.0:2.0;
  hmom_pull_ep_r = new TH1F("hmom_pull_ep_r","mom pull e^{+} (reco)",150,-3.*zoom,3.*zoom);
  hmom_pull_ep_f = new TH1F("hmom_pull_ep_f","mom pull e^{+} (fit)",150,-3.*zoom,3.*zoom);
  hmom_pull_em_r = new TH1F("hmom_pull_em_r","mom pull e^{-} (reco)",150,-3.*zoom,3.*zoom);
  hmom_pull_em_f = new TH1F("hmom_pull_em_f","mom pull e^{-} (fit)",150,-3.*zoom,3.*zoom);
  hpx_pull_ep_r = new TH1F("hpx_pull_ep_r","px pull e^{+} (reco)",150,-3.*zoom,3.*zoom);
  hpx_pull_ep_f = new TH1F("hpx_pull_ep_f","px pull e^{+} (fit)",150,-3.*zoom,3.*zoom);
  hpx_pull_em_r = new TH1F("hpx_pull_em_r","px pull e^{-} (reco)",150,-3.*zoom,3.*zoom);
  hpx_pull_em_f = new TH1F("hpx_pull_em_f","px pull e^{-} (fit)",150,-3.*zoom,3.*zoom);
  hpy_pull_ep_r = new TH1F("hpy_pull_ep_r","py pull e^{+} (reco)",150,-3.*zoom,3.*zoom);
  hpy_pull_ep_f = new TH1F("hpy_pull_ep_f","py pull e^{+} (fit)",150,-3.*zoom,3.*zoom);
  hpy_pull_em_r = new TH1F("hpy_pull_em_r","py pull e^{-} (reco)",150,-3.*zoom,3.*zoom);
  hpy_pull_em_f = new TH1F("hpy_pull_em_f","py pull e^{-} (fit)",150,-3.*zoom,3.*zoom);
  hpz_pull_ep_r = new TH1F("hpz_pull_ep_r","pz pull e^{+} (reco)",150,-3.*zoom,3.*zoom);
  hpz_pull_ep_f = new TH1F("hpz_pull_ep_f","pz pull e^{+} (fit)",150,-3.*zoom,3.*zoom);
  hpz_pull_em_r = new TH1F("hpz_pull_em_r","pz pull e^{-} (reco)",150,-3.*zoom,3.*zoom);
  hpz_pull_em_f = new TH1F("hpz_pull_em_f","pz pull e^{-} (fit)",150,-3.*zoom,3.*zoom);

}

void AnaTdav2::beam_cond(){
  //TLorentzVector ini = TLorentzVector(0, 0, 5.513, 6.53023); // need some info about
  if (mom_antip<0) {
    cout << "Antiproton momentum has to be given before intialization" << endl;
    exit(1);
  }
  const double mass_prot= 0.938;
  const double E_antip = TMath::Hypot(mass_prot, mom_antip);
  const double beta_cm = mom_antip/(E_antip + mass_prot);
  cout << "betac_cm = " << beta_cm << endl;
  boost_to_cm.SetZ(-beta_cm);
  boost_to_lab.SetZ(beta_cm);
  boost_to_cm.Print();
  boost_to_lab.Print();

  p4pbar.SetPxPyPzE(0,0,mom_antip,TMath::Hypot(mass_prot,mom_antip));
  p4targ.SetPxPyPzE(0,0,0,mass_prot);
  p4sys = p4pbar + p4targ;

  double ftmin[3] = {-0.45, -2.0, -5.0};
  double delta[3] = {0.1, 0.2, 0.2};
  //double ftmin[3] = {-0.1, -1.3, -2.8};
  //double delta[3] = {0.08, 0.2, 0.4};

  cout << "double range[iplab=" << iplab << "][]= {";
  for (int i=0; i<30; ++i) {
    double xx = ftmin[iplab] + delta[iplab]*i;
    //double xx = i==0 ? tmin[iplab]: (ftmin[iplab] + delta[iplab]*i);
    if (xx>tmax[iplab]) {
      cout << tmax[iplab] << " }; // n=" << i << endl;;
      tu_binning.push_back(tmax[iplab]);
      break;
    } else {
      tu_binning.push_back(xx);
      cout << xx << ",";
    }
  }
  print_binning(tu_binning, "tu_binning");

  tu_binning_2d.push_back(tmid[iplab]);
  tu_binning_2d.push_back(tmin[iplab]);
  //tu_binning_2d.push_back((tmin[iplab]+tmax[iplab])/2.0);
  tu_binning_2d.push_back(tmax[iplab]);
  //for (int icth=0; icth < 9; ++icth) { costh_binning_2d.push_back(-1.0 + (2.0*icth/8.0)); }
  for (int icth=0; icth < 9; ++icth) { costh_binning_2d.push_back(-0.8 + (2.0*icth/10.0)); }

  // Equal subdivisions in costh_cm, boost to lab for th bins
  TLorentzVector pi0;
  TVector3 u;
  for (int ii = 0; ii < tu_binning.size(); ++ii) {
    pi0cost_cm_binning.push_back(-1.0 + ( 2.0*ii/(tu_binning.size()-1) ));
  }

  for (int ii = pi0cost_cm_binning.size()-1; ii >=0; --ii) {
    u.SetMagThetaPhi(1,TMath::ACos(pi0cost_cm_binning[ii]),0);
    pi0.SetPxPyPzE(u.Px(),u.Py(),u.Pz(),TMath::Hypot(0.1349766,1.0));
    pi0.Boost(boost_to_lab);
    pi0th_binning.push_back( pi0.Theta());
  }
  pi0th_binning[0]-=0.0001; // safety for numerical error
  pi0th_binning[tu_binning.size()-1]+=0.0001; // safety for numerical error
  print_binning(pi0cost_cm_binning,"pi0cost_cm_binning");
  print_binning(pi0th_binning,"pi0th_binning");

}

void AnaTdav2::print_binning(const vector<double> &b, const char* n) {
  cout << n << ": {";
  for (int ii = 0; ii < b.size(); ++ii) cout << b[ii] << (ii==b.size()-1?"":", ");
  cout << "}; nbins= //" << b.size()-1 << endl;
}

InitStatus AnaTdav2::Init() {
  cout << "AnaTdav2::Init" << endl;
  fAna = new PndAnalysis();
  beam_cond();
  init_hists();
  init_tcas();
  return kSUCCESS;
}

void AnaTdav2::fill_mtot(RhoCandList& _ep, RhoCandList& _gg, TH1F* dest) {
  for (int jep = 0; jep < _ep.GetLength(); ++jep)
    for (int jgg = 0; jgg < _gg.GetLength(); ++jgg)
      dest->Fill(m(_gg[jgg], _ep[jep]), m_evt_wt);
}

void AnaTdav2::fill_mmiss(RhoCandList& _ep, RhoCandList& _gg, TH1F* dest, TH1F* dest2) {
  assert(_ep.GetLength()<=1&&_gg.GetLength()<=1);
  if (_ep.GetLength()==1&&_gg.GetLength()==1) {
    double _mmiss = (p4sys - (_gg[0]->P4() + _ep[0]->P4())).M();
    double _mmiss2 = (p4sys - (_gg[0]->P4() + _ep[0]->P4())).M2();
    dest->Fill(_mmiss, m_evt_wt);
    dest2->Fill(_mmiss2, m_evt_wt);
  }
}

void AnaTdav2::fill_mmiss_jpsi(RhoCandList& _ep, TH1F* dest, TH1F* dest2) {
  assert(_ep.GetLength()<=1);
  if (_ep.GetLength()==1) {
    double _mmiss = (p4sys -  _ep[0]->P4()).M();
    double _mmiss2 = (p4sys -  _ep[0]->P4()).M2();
    dest->Fill(_mmiss, m_evt_wt);
    dest2->Fill(_mmiss2, m_evt_wt);
  }
}

void AnaTdav2::fill_dth_dph_cm(RhoCandList& _ep, RhoCandList& _gg, TH2F* dest) {
  double _dth_cm=0, _dph_cm=0;
  for (int jep = 0; jep < _ep.GetLength(); ++jep)
    for (int jgg = 0; jgg < _gg.GetLength(); ++jgg) {
      dth_dph_cm(_gg[jgg], _ep[jep], _dth_cm, _dph_cm);
      dest->Fill(_dth_cm, _dph_cm);
    }
}

// Here, weight that is applied is the "event" weight. For mutli-pair events all pairs get the same weight.
// It is possible to do better, by calculating the weight for each pair based on the MC truth info
// but this will require big changes in structure
void AnaTdav2::fill_pair_mass(RhoCandList& org, TH1F* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->M(),m_evt_wt);
}

// Here, weight that is applied is the "event" weight. For mutli-pair events all pairs get the same weight.
// It is possible to do better, by calculating the weight for each pair based on the MC truth info
// but this will require big changes in structure
void AnaTdav2::fill_pair_mass(RhoCandList& org, TH1F* dest, double more_weight) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->M(),m_evt_wt*more_weight);
}

bool AnaTdav2::check_mct_jpsi(RhoCandidate* _cand) {
  if (mct_uid_e<0||mct_uid_p<0) return false;
  if (_cand->NDaughters()!=2) return false;
  int i0 = _cand->Daughter(0)->Uid();
  int i1 = _cand->Daughter(1)->Uid();
  return (i0==mct_uid_e&&i1==mct_uid_p)||(i0==mct_uid_p&&i1==mct_uid_e);
}

bool AnaTdav2::check_mct_pi0(RhoCandidate* _cand) {
  if (mct_uid_g1<0||mct_uid_g2<0) return false;
  if (_cand->NDaughters()!=2) return false;
  int i0 = _cand->Daughter(0)->Uid();
  int i1 = _cand->Daughter(1)->Uid();
  return (i0==mct_uid_g1&&i1==mct_uid_g2)||(i0==mct_uid_g2&&i1==mct_uid_g1);
}

double AnaTdav2::the_bwd(RhoCandidate* _cand) {
  if (_cand->NDaughters()!=2) return -9999.;
  double the0 = _cand->Daughter(0)->P3().Theta();
  double the1 = _cand->Daughter(1)->P3().Theta();
  return TMath::RadToDeg()*(the0>the1?the0:the1);
}

double AnaTdav2::the_fwd(RhoCandidate* _cand) {
  if (_cand->NDaughters()!=2) return -9999.;
  double the0 = _cand->Daughter(0)->P3().Theta();
  double the1 = _cand->Daughter(1)->P3().Theta();
  return TMath::RadToDeg()*(the0<the1?the0:the1);
}

double AnaTdav2::oa(RhoCandidate* _cand) {
  if (_cand->NDaughters()!=2) return -9999.;
  RhoCandidate *d1 = _cand->Daughter(0);
  RhoCandidate *d2 = _cand->Daughter(1);
  return oa(d1,d2);
}

void AnaTdav2::fill_mctruth(RhoCandList& _ep, RhoCandList& _gg, int step) {
  for (int ii=0; ii < _ep.GetLength(); ++ii) {
    for (int jj=0; jj < _gg.GetLength(); ++jj) {
      if (check_mct_jpsi(_ep[ii]) and
	  check_mct_pi0(_gg[jj])) {
	hmep_mct[step]->Fill(_ep[ii]->M(), m_evt_wt);
	hthe_ep_mct[step]->Fill(_ep[ii]->P3().Theta()*TMath::RadToDeg(), m_evt_wt);
	hthe_ep_mct_fwd[step]->Fill(the_fwd(_ep[ii]), m_evt_wt);
	hthe_ep_mct_bwd[step]->Fill(the_bwd(_ep[ii]), m_evt_wt);
	hmgg_mct[step]->Fill(_gg[jj]->M(),m_evt_wt);
	hthe_gg_mct[step]->Fill(_gg[jj]->P3().Theta()*TMath::RadToDeg(), m_evt_wt);
	hthe_gg_mct_fwd[step]->Fill(the_fwd(_gg[jj]), m_evt_wt);
	hthe_gg_mct_bwd[step]->Fill(the_bwd(_gg[jj]), m_evt_wt);
	hoa_gg_mct[step]->Fill(oa(_gg[jj]), m_evt_wt);
      } else {
	hmep_non_mct[step]->Fill(_ep[ii]->M(), m_evt_wt);
	hmgg_non_mct[step]->Fill(_gg[jj]->M(), m_evt_wt);
      }
    }
  }
}

void AnaTdav2::fill_count_hists(int _gg, int _ep, int ihist) {
  hnep[ihist]->Fill(rcl[_ep].GetLength());
  hngg[ihist]->Fill(rcl[_gg].GetLength());
  hnpi0jpsi[ihist]->Fill(rcl[_gg].GetLength()*rcl[_ep].GetLength());
}

void AnaTdav2::print_indices() {
  cout << "e= " << e << " p= " <<  p << " g= " <<  g
       << "ie = " << ie<< " ip = " <<  ip<< " gg = " <<  gg<< " gg_sel = "
       <<  gg_sel<< " ep = " <<  ep<< " iep = " << iep<< " iep_uniq = "
       << iep_uniq<< " iep_asso = " <<  iep_asso<< " iep_excl = "
       << iep_excl<< " gg_excl = " <<  gg_excl << " nrcl= " << nrcl << endl;
}

bool AnaTdav2::check_eid(RhoCandidate* cand) {
  if (eid_param_method) { // Eid Paramterisation
    cand->SetType(-11*cand->Charge());
    double mom = cand->P3().Mag();
    double theta = cand->P3().Theta();
    if (fAna->McTruthMatch(cand)) {
      mom = cand->GetMcTruth()->P3().Mag();
      theta = cand->GetMcTruth()->P3().Theta();
    }
    if (!eff_hist_rad) theta *= TMath::RadToDeg();
    const int ix = heff_epm->GetXaxis()->FindBin(mom);
    const int iy = heff_epm->GetYaxis()->FindBin(theta);
    const double eff = heff_epm->GetBinContent(ix,iy);
    const double prob = gRandom->Uniform();
    //cout << "eff = " << eff << endl;
    return (prob<eff);
  } else { // Ronald's Method
    return bayes_pid(cand);
  }
}

inline
double AnaTdav2::dist_photon_match(RhoCandidate *rec, RhoCandidate *mc) {
  const double _dth = (rec->P3().Theta()-mc->P3().Theta())/8.0e-3;
  const double _dph = (rec->P3().Phi()-mc->P3().Phi())/8.22e-3;
  return TMath::Hypot(_dph, _dth) / 3.0;
}

inline
double AnaTdav2::dist_chpi_match(RhoCandidate *rec, RhoCandidate *mc) {
  const double _dth = (rec->P3().Theta()-mc->P3().Theta())/1.54e-3;
  const double _dph = (rec->P3().Phi()-mc->P3().Phi())/3.95e-3;
  return TMath::Hypot(_dth, _dph);
}

void AnaTdav2::charged_pion_filter(RhoCandList& outp, RhoCandList& outm, RhoCandList& inp_pi, RhoCandList& inm_pi, RhoCandList& inp_e, RhoCandList& inm_e) {
  // Assume tracks returned by pion and electrn hypothesis come in the same order
  // The least requirement is that they have the same number of tracks
  assert(inp_pi.GetLength()==inp_e.GetLength());
  assert(inm_pi.GetLength()==inm_e.GetLength());
  // Quck test of "identity"-- seems like everythign is legit here
  //for (int ii = 0; ii < inp_pi.GetLength(); ++ii) assert(inp_pi[ii]->GetTrackNumber()==inp_e[ii]->GetTrackNumber());
  //for (int ii = 0; ii < inm_pi.GetLength(); ++ii) assert(inm_pi[ii]->GetTrackNumber()==inm_e[ii]->GetTrackNumber());
  //for (int ii = 0; ii < inm_pi.GetLength(); ++ii) {
  //  cout << "p(pihyp)= " << inm_pi[ii]->P3().Mag() << " p(elhyp)= " << inm_e[ii]->P3().Mag() << endl;;
  //  cout << "E(pihyp)= " << inm_pi[ii]->Energy() << " E(elhyp)= " << inm_e[ii]->Energy() << endl;;
  //}

  int pdg_pip = 211, pdg_pim = -211;
  int mc_p = -1, mc_m = -1;
  for (int j=0;j<mcList.GetLength();++j) {
    if (mcList[j]->PdgCode()!=pdg_pip and mcList[j]->PdgCode()!=pdg_pim) continue;
    if (!mcList[j]->TheMother()) {
      if (mcList[j]->PdgCode()==pdg_pip) mc_p = j;
      if (mcList[j]->PdgCode()==pdg_pim) mc_m = j;
    }
    if (mc_p>=0&&mc_m>=0) break;
  }
  // These shouldn't really happen!
  assert(0 <= mc_p and mc_p < mcList.GetLength());
  assert(0 <= mc_m and mc_m < mcList.GetLength());

  if (inp_pi.GetLength()==0 or inm_pi.GetLength()==0) return; // not interested in such events
  if (inp_pi.GetLength()==1 and inm_pi.GetLength()==1) { // accept such events without condition...
    outp.Append(inp_e[0]);
    outm.Append(inm_e[0]);
    return;
  }

  // Here keep the cosest matching pair
  double dmin = 1e9;
  int match_p = -1, match_m = -1;
  for (int iip = 0; iip < inp_pi.GetLength(); ++iip) {
    for (int iim = 0; iim < inm_pi.GetLength(); ++iim) {
      const double dp = dist_chpi_match(inp_pi[iip],mcList[mc_p]);
      const double dm = dist_chpi_match(inm_pi[iim],mcList[mc_m]);
      const double dtot = TMath::Hypot(dp,dm);
      if (dtot < dmin) {
	match_p = iip;
	match_m = iim;
	dmin = dtot;
      }
    }
  }
  if (dmin<10000) {
    outp.Append(inp_e[match_p]);
    outm.Append(inm_e[match_m]);
  } else {
    cout << "pi+ pi- match couldn't be found because reco tracks are too far away from MC tracks dist= " << dmin << endl;
  }
}

void AnaTdav2::mctruth_match_jpsi(RhoCandList& elec, RhoCandList& posit) {
  mct_uid_e = -1;
  mct_uid_p = -1;
  if (elec.GetLength()==0 || posit.GetLength()==0) return;

  int pdg_jpsi = 443, pdg_elec = 11, pdg_posit = -11;
  //int mc_elec = -1, mc_posit = -1;
  mc_elec = -1;
  mc_posit = -1;
  for (int j=0;j<mcList.GetLength();++j) {
    if (mcList[j]->PdgCode()!=pdg_elec and mcList[j]->PdgCode()!=pdg_posit) continue;
    if (!mcList[j]->TheMother()) continue;
    if (mcList[j]->TheMother()->PdgCode()==pdg_jpsi) {
      if (mcList[j]->PdgCode()==pdg_elec) mc_elec = j;
      if (mcList[j]->PdgCode()==pdg_posit) mc_posit = j;
    }
    if (mc_elec>=0&&mc_posit>=0) break;
  }
  if (mc_elec<0||mc_posit<0) {
    cout << "McTrue e+e- from jpsi not found!!" << endl;
    return;
  }

  // These shouldn't really happen!
  assert(0 <= mc_elec and mc_elec < mcList.GetLength());
  assert(0 <= mc_posit and mc_posit < mcList.GetLength());

  double dmin = 1e9;
  int match_elec = -1, match_posit = -1;
  for (int iielec = 0; iielec < elec.GetLength(); ++iielec) {
    for (int iiposit = 0; iiposit < posit.GetLength(); ++iiposit) {
      const double delec = dist_chpi_match(elec[iielec],mcList[mc_elec]);
      const double dposit = dist_chpi_match(posit[iiposit],mcList[mc_posit]);
      const double dtot = TMath::Hypot(delec,dposit);
      if (dtot < dmin) {
	match_elec = iielec;
	match_posit = iiposit;
	dmin = dtot;
      }
    }
  }
  if (dmin<10000) {
    mct_uid_e = elec[match_elec]->Uid();
    mct_uid_p = posit[match_posit]->Uid();
  } else {
    cout << "e+ e- match couldn't be found because reco tracks are too far away from MC tracks dist= " << dmin << endl;
  }

}

void AnaTdav2::mctruth_match_pi0(RhoCandList& gamma) {
  mct_uid_g1 = -1;
  mct_uid_g2 = -1;
  if (gamma.GetLength()<1) return;

  int pdg_pi0 = 111, pdg_gamma = 22;
  //int mc_g1 = -1, mc_g2 = -1;
  mc_g1 = -1;
  mc_g2 = -1;
  for (int j=0;j<mcList.GetLength();++j) {
    if (mcList[j]->PdgCode()!=pdg_gamma) continue;
    if (!mcList[j]->TheMother()) continue;
    if (mcList[j]->TheMother()->PdgCode()==pdg_pi0) {
      if (mcList[j]->PdgCode()==pdg_gamma) {
	if (mc_g1 == -1) mc_g1 = j;
	else mc_g2 = j;
      }
    }
    if (mc_g1>=0&&mc_g2>=0) break;
  }
  if (mc_g1<0||mc_g1<0) {
    cout << "McTrue gg from pi0 not found!!" << endl;
    return;
  }

  // These shouldn't really happen!
  assert(0 <= mc_g1 and mc_g1 < mcList.GetLength());
  assert(0 <= mc_g2 and mc_g2 < mcList.GetLength());

  double dmin = 1e9;
  int match_g1 = -1, match_g2 = -1;
  for (int iig1 = 0; iig1 < gamma.GetLength(); ++iig1) {
    for (int iig2 = 0; iig2 < gamma.GetLength(); ++iig2) {
      const double dg1 = dist_photon_match(gamma[iig1],mcList[mc_g1]);
      const double dg2 = dist_photon_match(gamma[iig2],mcList[mc_g2]);
      const double dtot = TMath::Hypot(dg1,dg2);
      if (dtot < dmin) {
	match_g1 = iig1;
	match_g2 = iig2;
	dmin = dtot;
      }
    }
  }
  if (dmin<10000) {
    mct_uid_g1 = gamma[match_g1]->Uid();
    mct_uid_g2 = gamma[match_g2]->Uid();
  } else {
    cout << "pi0->gg match couldn't be found because reco tracks are too far away from MC tracks dist= " << dmin << endl;
  }

}

void AnaTdav2::eid_filter(RhoCandList&out, RhoCandList&in) {
  for (int i=0; i< in.GetLength(); ++i) {
    if ( check_eid(in[i]) )
      out.Append(in[i]);
  }
}

double AnaTdav2::eff_weight(const TVector3 &mom) {
  int gb = 0;
  if (pi_eff->GetDimension()==2) {
    int bx = pi_eff->GetTotalHistogram()->GetXaxis()->FindBin(mom.Mag());
    int by = pi_eff_hist_rad?
      pi_eff->GetTotalHistogram()->GetXaxis()->FindBin(mom.Theta()*TMath::RadToDeg()):
      pi_eff->GetTotalHistogram()->GetXaxis()->FindBin(mom.Theta());
    gb = pi_eff->GetGlobalBin(bx,by);
  } else {
    if (pi_eff_func) {
      if (mom.Mag()<0.1)
	return pi_eff_func->Eval(0.1);
      if (mom.Mag()>10)
	return pi_eff_func->Eval(10);
      else
	return pi_eff_func->Eval(mom.Mag());
    }
    int bx = pi_eff->GetTotalHistogram()->GetXaxis()->FindBin(mom.Mag());
    gb = pi_eff->GetGlobalBin(bx);
  }
  return pi_eff->GetEfficiency(gb);
}

//Signal
//Track 0 (PDG:88888) has mother -1 and daughter(s) 1  2
//Track 1 (PDG:443) has mother 0 and daughter(s) 3  4
//Track 2 (PDG:111) has mother 0 and daughter(s) 5  6
//Track 3 (PDG:-11) has mother 1 and daughter(s) 909  910  ...
//Track 4 (PDG:11) has mother 1 and daughter(s) 178  179  ...
//Track 5 (PDG:22) has mother 2 and daughter(s) 23  177
//Track 6 (PDG:22) has mother 2 and daughter(s) 7  22
//Background
//Track 0 (PDG:211) has mother -1 and daughter(s) 857  961
//Track 1 (PDG:111) has mother -1 and daughter(s) 2  3
//Track 2 (PDG:22) has mother 1 and daughter(s) 157  856
//Track 3 (PDG:22) has mother 1 and daughter(s) 6  156
bool AnaTdav2::calc_true_tu() {
  for (int j=0;j<mcList.GetLength();++j) {
    if (mcList[j]->PdgCode()==111) { // in case of multi pi0 events just use the first pi0 to define "true" t and u. Anyways, these don't make much sense
      RhoCandidate *mcmother = mcList[j]->TheMother();
      int muid = mcmother? mcmother->GetTrackNumber(): -1;
      //if (( (mc_type!=1&&mc_type!=4&&mc_type!=5&&mc_type!=6) and muid == -1) or  // DPM simulation doesn't have pbar-p system at index 0
      if (( is_dpm() and muid == -1) or // DPM simulation doesn't have pbar-p system at index 0
	  ( is_evt_gen() and muid == 0)) {  // EvtGen simulation has pbar-p system at index 0
	event_t = t_gg(mcList[j]);
	event_u = u_gg(mcList[j]);
	event_pi0costh_cm = pi0cost_cm(mcList[j]);
	event_pi0theta_cm = pi0theta_cm(mcList[j]);
	event_pi0theta_lab = mcList[j]->P4().Theta();

	int itbin_2d = find_bin(event_t, tu_binning_2d);
	int iubin_2d = find_bin(event_u, tu_binning_2d);
	int itu2d = itbin_2d>=0?itbin_2d:(iubin_2d>=0?2+iubin_2d:-1);

	m_epcth_wt0 = 1.0;
	m_epcth_wt1 = 1.0;
	event_epthe_jpsi = -9999.0;
	event_epcosth_jpsi = -9999.0;
	if (is_evt_gen()) {
	  for (int ik=0; ik < mcList.GetLength(); ++ik) {
	    if (mcList[ik]->PdgCode()==443) {
	      event_epthe_jpsi = the_b(get_p4ep(mcList[ik]), -mcList[ik]->P4().BoostVector());
	      event_epcosth_jpsi = cost_b(get_p4ep(mcList[ik]), -mcList[ik]->P4().BoostVector());
	      double _epth_lab = get_p4ep(mcList[ik]).Vect().Theta() *TMath::RadToDeg();
	      double _emth_lab = get_p4em(mcList[ik]).Vect().Theta() *TMath::RadToDeg();

	      m_epcth_wt0 = (1 + event_epcosth_jpsi*event_epcosth_jpsi)/1.5;
	      m_epcth_wt1 = (1 + 0.4*event_epcosth_jpsi*event_epcosth_jpsi)/1.2;

	      hepcosth_jpsi_mc_all->Fill(event_epcosth_jpsi);
	      hepcosth_jpsi_mc_all_wt0->Fill(event_epcosth_jpsi, m_epcth_wt0);
	      hepcosth_jpsi_mc_all_wt1->Fill(event_epcosth_jpsi, m_epcth_wt1);
	      hepcosth_jpsi_vs_epthlab_mc_all->Fill(event_epcosth_jpsi, _epth_lab  );
	      hepcosth_jpsi_vs_emthlab_mc_all->Fill(event_epcosth_jpsi, _emth_lab  );

	      if (itu2d>=0){
		hepcosth_jpsi_mc[itu2d]->Fill(event_epcosth_jpsi);
		hepcosth_jpsi_vs_epthlab_mc[itu2d]->Fill(event_epcosth_jpsi, _epth_lab);
		hepcosth_jpsi_vs_emthlab_mc[itu2d]->Fill(event_epcosth_jpsi, _emth_lab);
	      }

	    }
	  }
	}

	httrumc->Fill(event_t);
	hutrumc->Fill(event_u);
	httrumc_vb->Fill(event_t);
	hutrumc_vb->Fill(event_u);
	htrupi0thcm->Fill(event_pi0costh_cm);
	htrupi0costhcm->Fill(event_pi0theta_cm);
	htrupi0thlab->Fill(event_pi0theta_lab);

	bool t_ok = (tmin[iplab] < event_t && event_t < tmax[iplab]);
	bool u_ok = (tmin[iplab] < event_u && event_u < tmax[iplab]);

	if (t_ok || u_ok) {
	  htrupi0thcm_tcut->Fill(event_pi0theta_cm);
	  htrupi0costhcm_tcut->Fill(event_pi0costh_cm);
	  htrupi0thlab_tcut->Fill(event_pi0theta_lab);
	}

	// MC pi0 angluar distributions with jpsi mass cut on the pippim pair
	if (is_dpm()) {
	  int kpip = -1, kpim=-1;
	  for (int k=0; k<mcList.GetLength(); ++k) {
	    if (mcList[k]->PdgCode()==211) { kpip = k; }
	    if (mcList[k]->PdgCode()==-211) { kpim = k; }
	    if (kpip>0&&kpim>0) { // in case of multi pip+pi- events, use the first found pi+pi- event to define mass cuts
	      double mpippim = (mcList[kpip]->P4()+mcList[kpim]->P4()).M();
	      htrupi0thcm_vs_m->Fill(mpippim, event_pi0theta_cm);
	      htrupi0costhcm_vs_m->Fill(mpippim, event_pi0costh_cm);
	      htrupi0thlab_vs_m->Fill(mpippim, event_pi0theta_lab);
	      if ( mpippim>jpsi_m_3sig_min && mpippim<jpsi_m_3sig_max) {
		htrupi0thcm_mcut->Fill(event_pi0theta_cm);
		htrupi0costhcm_mcut->Fill(event_pi0costh_cm);
		htrupi0thlab_mcut->Fill(event_pi0theta_lab);

		if (t_ok || u_ok) {
		  htrupi0thcm_tcut_mcut->Fill(event_pi0theta_cm);
		  htrupi0costhcm_tcut_mcut->Fill(event_pi0costh_cm);
		  htrupi0thlab_tcut_mcut->Fill(event_pi0theta_lab);
		}

		htrupi0thcm_mcut_vs_m->Fill(mpippim, event_pi0theta_cm);
		htrupi0costhcm_mcut_vs_m->Fill(mpippim, event_pi0theta_cm);
		htrupi0thlab_mcut_vs_m->Fill(mpippim, event_pi0theta_lab);
	      }
	      break;
	    }
	  }
	} else {
	  if (t_ok || u_ok) {
	    htrupi0thcm_tcut_mcut->Fill(pi0theta_cm(mcList[j]));
	    htrupi0costhcm_tcut_mcut->Fill(pi0cost_cm(mcList[j]));
	    htrupi0thlab_tcut_mcut->Fill(mcList[j]->P4().Theta());
	  }
	}

	return true;

      }
    }
  }
  cout << "MCtrue Pi0 not found. Not normal This shouldn't happen" << endl;
  return false;
}

inline
bool comp_dth(std::pair<pair<int,int>,double> lhs, std::pair<pair<int,int>,double> rhs) {
  return lhs.second<rhs.second;
}

void AnaTdav2::calc_evt_wt() {
  if (mc_type==0 /*pi0pipm*/||mc_type==2 /*pi0pi0pipm*/) {
    bool pip_found = false, pim_found = false;
    for (int j=0;j<mcList.GetLength();++j) {
      if (mcList[j]->PdgCode()==211&&!pip_found) {
	m_pip_wt = eff_weight(mcList[j]->GetMomentum());
	pip_found = true;
      }
      if (mcList[j]->PdgCode()==-211&&!pim_found) {
	m_pim_wt = eff_weight(mcList[j]->GetMomentum());
	pim_found = true;
      }
      if (pip_found&&pim_found) break;
    }
    m_evt_wt = nevt_xsect[mc_type][iplab]*m_pip_wt*m_pim_wt/nevt_sim[mc_type][iplab];
  } else if (mc_type==1 /*pi0jpsi->epm*/ || mc_type==5 /*pi0jpsi_10cm*/ || mc_type==6 /*pi0jpsi_5cm*/ ){
    m_evt_wt = 1.0;
  } else if (mc_type==4 /*pi0pi0jpsi->epm*/) {
    m_evt_wt = nevt_xsect[mc_type][iplab]/nevt_sim[mc_type][iplab];
  } else if (mc_type==3 /*pi0pipmpippim*/) {
    // find the most back to back pip-pim pair to the pi0, and use that to set the event weight
    std::vector<TLorentzVector> p4pip, p4pim, p4pi0;
    for (int j=0;j<mcList.GetLength();++j) {
      if (mcList[j]->PdgCode()==211 && p4pip.size()<2) p4pip.push_back(mcList[j]->P4());
      if (mcList[j]->PdgCode()==-211 && p4pim.size()<2) p4pim.push_back(mcList[j]->P4());
      if (mcList[j]->PdgCode()==111 && p4pi0.size()<1) p4pi0.push_back(mcList[j]->P4());
      if (p4pip.size()==2 && p4pim.size()==2 && p4pi0.size()==1) break;
    }
    assert(p4pip.size()==2 && p4pim.size()==2 && p4pi0.size()==1); // if not, wrong type of event
    std::vector<std::pair<std::pair<int,int>,double> > combs;
    for (int iip=0; iip<2; iip++)
      for (int iim=0; iim<2; iim++)
	combs.push_back(make_pair(make_pair(iip,iim),dth_cm(p4pi0[0],p4pip[iip]+p4pim[iim])));
    sort(combs.begin(),combs.end(),comp_dth);
    m_pip_wt = eff_weight(p4pip[combs[0].first.first].Vect());
    m_pim_wt = eff_weight(p4pim[combs[0].first.second].Vect());
    m_evt_wt = nevt_xsect[mc_type][iplab]*m_pip_wt*m_pim_wt/nevt_sim[mc_type][iplab];
  } else {
    m_evt_wt = 1.0;
  }
  hwt->Fill(m_evt_wt);
}

void AnaTdav2::print_mc_list() {
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

/**
 * Fills reconstructed single particle lists
 * returns true
 */
void AnaTdav2::fill_lists() {
  //print_indices();
  cleanup_lists();
  mcList.Cleanup();

  // *** Select with no PID info ('All'); type and mass are set
  fAna->FillList(mcList, "McTruth");
  fAna->FillList(rcl[p], (brem_corr?"BremElectronAllPlus":"ElectronAllPlus"));
  fAna->FillList(rcl[e], (brem_corr?"BremElectronAllMinus":"ElectronAllMinus"));
  fAna->FillList(rcl[g], "Neutral");

  if (is_dpm()) { // if final state is not electrons, use efficiency weight
    fAna->FillList(rcl[pip], "PionAllPlus");
    fAna->FillList(rcl[pim], "PionAllMinus");
    // all pions accepted as identified, but they will be wieghted when filling with the corresponding
    // product of efficiency of the two tracks
    fAna->FillList(rcl[ip], (brem_corr?"BremElectronAllPlus":"ElectronAllPlus"));
    fAna->FillList(rcl[ie], (brem_corr?"BremElectronAllMinus":"ElectronAllMinus"));

    //fAna->FillList(rcl[ip], (brem_corr?"BremElectronVeryTightPlus":"ElectronAllPlus"), "PidAlgoEmcBayes");
    //fAna->FillList(rcl[ie], (brem_corr?"BremElectronVeryTightMinus":"ElectronAllMinus"), "PidAlgoEmcBayes");
    //fAna->FillList(rcl[ip], (brem_corr?"BremElectronAllPlus":"ElectronAllPlus"));
    //fAna->FillList(rcl[ie], (brem_corr?"BremElectronAllMinus":"ElectronAllMinus"));
    //charged_pion_filter(rcl[ip], rcl[ie], rcl[pip], rcl[pim], rcl[p], rcl[e]);
  } else {
    eid_filter(rcl[ip],rcl[p]);
    eid_filter(rcl[ie],rcl[e]);
  }

  // If signal sim, tag tracks from true jpsi and pi0
  if (mc_type==1||mc_type==4||mc_type==5) {
    mctruth_match_pi0(rcl[g]);
    mctruth_match_jpsi(rcl[e], rcl[p]);
  }

  hng->Fill(rcl[g].GetLength());
  ng20mev = 0;
  for (int ig=0; ig < rcl[g].GetLength(); ++ig) {
    if (rcl[g][ig]->Energy()>0.02) ng20mev++;
  }
  hng20mev->Fill(ng20mev);
  hnch->Fill(rcl[e].GetLength()+rcl[p].GetLength());

}

bool AnaTdav2::is_dpm() {
  return mc_type==0||mc_type==2||mc_type==3;
}

bool AnaTdav2::is_evt_gen() {
  return mc_type==1||mc_type==4||mc_type==5||mc_type==6;
}

void AnaTdav2::nocut_ref() {
  rcl[gg].Combine(rcl[g],rcl[g]);
  rcl[ep].Combine(rcl[e],rcl[p]);
  double tmp_eid_wt = m_evt_wt*nevt_sim[mc_type][iplab]/nevt_xsect[mc_type][iplab];
  if (is_dpm()) m_evt_wt = nevt_xsect[mc_type][iplab]/nevt_sim[mc_type][iplab];
  fill_pair_mass(rcl[ep], hmep[0]);
  fill_pair_mass(rcl[gg], hmgg[0]);
  fill_mctruth(rcl[ep], rcl[gg], 0);
  fill_count_hists(gg,ep,0);
  if (is_dpm()) m_evt_wt *= tmp_eid_wt;

  rcl[iep].Combine(rcl[ie],rcl[ip]);

  fill_pair_mass(rcl[iep], hmep[1]);
  fill_pair_mass(rcl[gg], hmgg[1]);
  fill_mctruth(rcl[iep], rcl[gg], 1);
  fill_count_hists(gg,iep,1);
}

/**
 * Makes a pid selected e-p pair list on the condition that there is
 * an exclusive e-p pair in a reasonable mass range, events with more than one
 * pair in a "reasonable" mass range are rejected. The pair found within that
 * mass range is selected for further analysis
 */
void AnaTdav2::ep_uniq() {
  if (rcl[iep].GetLength() == 1) {
    assert( rcl[ie].GetLength()<=1 and rcl[ip].GetLength()<=1 );
    rcl[iep_uniq].Combine(rcl[ie],rcl[ip]);
  }
  fill_pair_mass(rcl[iep_uniq], hmep[2]);
  fill_pair_mass(rcl[gg], hmgg[2]);
  fill_mctruth(rcl[iep_uniq], rcl[gg], 2);
  fill_count_hists(gg,iep_uniq,2);
}

/**
 * This alternative analysis combines all reconstructed pairs and leaves
 * it to later stages to select the best matching pi0-jpsi pair or fill all
 * of them with a weight
 */
void AnaTdav2::ep_all() {
  rcl[iep_all].Combine(rcl[ie],rcl[ip]);
  fill_pair_mass(rcl[iep_all], hmep[2]);
  fill_pair_mass(rcl[gg], hmgg[2]);
  fill_mctruth(rcl[iep_all],rcl[gg], 2);
  fill_count_hists(gg,iep_all,2);
}

bool AnaTdav2::oa_vs_avg_cut(const double& _oa, const double &_avg){
  bool acc = (_avg > (lw[iplab][2] + (lw[iplab][0]/(_oa-lw[iplab][1]))));
  if (_oa>up[iplab][1])
    acc = acc && (_avg < (up[iplab][2] + (up[iplab][0] / (_oa-up[iplab][1]))));
  return acc;
}

/**
 * Makes a sublist of gg pairs that satisfies minimal pi0 selection cuts
 * that reduce combinatorial
 */
void AnaTdav2::pi0_sel() {
  for (int i = 0; i < rcl[gg].GetLength(); ++i) {
    RhoCandidate *_g1 = rcl[gg][i]->Daughter(0);
    RhoCandidate *_g2 = rcl[gg][i]->Daughter(1);
    const double _oa = oa(_g1,_g2);
    const double _e_avg = 0.5*(_g1->Energy()+_g2->Energy());
    if (apply_pi0evsoa_cut and
	!oa_vs_avg_cut(_oa, _e_avg) ) continue;
    if (apply_pi0m_cut and
	(rcl[gg][i]->M() < pi0m_cut_min or rcl[gg][i]->M() > pi0m_cut_max) ) continue;
    rcl[gg_sel].Append(rcl[gg][i]);
  }
  fill_count_hists(gg_sel,iep_uniq,3);
}

/**
 * Looks in the event if there is an associated gg pair which
 * passes the pi0 selection cuts (loose cuts at this point)
 */
void AnaTdav2::ep_pi0_asso() {
  assert(rcl[iep_uniq].GetLength()<=1);
  for (int i = 0; i < rcl[iep_uniq].GetLength(); ++i) {
    if (rcl[gg_sel].GetLength()>0) {
      rcl[iep_asso].Append(rcl[iep]);
    }
  }
  fill_pair_mass(rcl[iep_asso], hmep[3]);
  fill_pair_mass(rcl[gg_sel], hmgg[3]);
  fill_mctruth(rcl[iep_asso], rcl[gg_sel], 3);
  fill_count_hists(gg_sel,iep_asso,4);
}

/**
 * Looks in the event if there is an associated gg pair which
 * passes the pi0 selection cuts (loose cuts at this point). Here
 * the assumption of uniquness of ep pair is not made. All possible
 * pairs are kept
 */
void AnaTdav2::ep_pi0_asso_all() {
  if (rcl[gg_sel].GetLength()>0) {
    rcl[iep_asso_all].Append(rcl[iep_all]);
  }
  fill_pair_mass(rcl[iep_asso_all], hmep[3]);
  fill_pair_mass(rcl[gg_sel], hmgg[3]);
  fill_mctruth(rcl[iep_asso_all], rcl[gg_sel], 3);
  fill_count_hists(gg_sel,iep_asso_all,4);
}

/**
 * After applying kinematic cuts, (dPhi, dTh and mTot), checks that there
 * is only one pi0 left in the event. This is to insure exclusivity. The
 * unique pi0 is appended to the exgg, and the e-p pair is filled to exep
 */
void AnaTdav2::kin_excl() {
  assert(rcl[iep_asso].GetLength()<=1);
  double _dth_dph_cm_min = 1e9;
  int jgg_most_btb = -1;
  for (int jep = 0; jep < rcl[iep_asso].GetLength(); ++jep) {
    int ngg_ok=0, jgg_ok;
    // count gg_sel's that satisfy the kinematic and exclusivity cuts
    for (int jgg = 0; jgg < rcl[gg_sel].GetLength(); ++jgg) {
      double _mtot = m(rcl[gg_sel][jgg], rcl[iep_asso][jep]);
      double _dth_cm=0, _dph_cm=0;
      dth_dph_cm(rcl[gg_sel][jgg], rcl[iep_asso][jep], _dth_cm, _dph_cm);
      double _dth_dph_cm = hypot((_dth_cm-TMath::Pi())/dth_sigma, (_dph_cm-TMath::Pi())/dph_sigma);
      if (_dth_dph_cm<_dth_dph_cm_min) {
	_dth_dph_cm_min = _dth_dph_cm;
	jgg_most_btb = jgg;
      }
      if (apply_dth_dph_cut and
	  (_dth_dph_cm > dth_dph_cm_cut_max) ) continue;
      if (apply_mtot_cut and
	  ( (_mtot < mtot_cut_min) or (_mtot > mtot_cut_max) ) ) continue;
      ngg_ok++;
      jgg_ok = jgg;
    }
    if (require_exclusivity) {
      if (ngg_ok==1) {
	rcl[iep_excl].Append(rcl[iep_asso][jep]);
	rcl[gg_excl].Append(rcl[gg_sel][jgg_ok]);
      }
    } else { // just pick the most back to back
      rcl[iep_excl].Append(rcl[iep_asso][jep]);
      rcl[gg_excl].Append(rcl[gg_sel][jgg_most_btb]);
    }
  }
  fill_dth_dph_cm(rcl[iep_excl],rcl[gg_excl], hcmoa);
  fill_mtot(rcl[iep_excl],rcl[gg_excl], hmtot);
  fill_mmiss(rcl[iep_excl],rcl[gg_excl], hmmiss, hmmiss2);
  fill_mmiss_jpsi(rcl[iep_excl], hmmiss_jpsi, hmmiss2_jpsi);
  fill_pair_mass(rcl[iep_excl], hmep[4]);
  fill_pair_mass(rcl[gg_excl], hmgg[4]);
  fill_mctruth(rcl[iep_excl], rcl[gg_excl], 4);
  fill_count_hists(gg_excl,iep_excl,5);
}

/**
 * Select the most back to back and the mclosest to sqrt(s) pair
 * Check that it passes the cut windows on dTh and mtot
 */
void AnaTdav2::kin_excl_all() {
  double _dth_dph_cm_min = 1e9;
  int jgg_most_btb = -1;
  int jep_most_btb = -1;

  double _dmtot_min = 1e9;
  int jgg_dmtot_best = -1;
  int jep_dmtot_best = -1;

  for (int jep = 0; jep < rcl[iep_asso_all].GetLength(); ++jep) {
    int ngg_ok=0, jgg_ok;
    // count gg_sel's that satisfy the kinematic and exclusivity cuts
    for (int jgg = 0; jgg < rcl[gg_sel].GetLength(); ++jgg) {
      double _mtot = m(rcl[gg_sel][jgg], rcl[iep_asso_all][jep]);
      double _dmtot = _mtot - p4sys.M();
      double _dth_cm=0, _dph_cm=0;
      dth_dph_cm(rcl[gg_sel][jgg], rcl[iep_asso_all][jep], _dth_cm, _dph_cm);
      double _dth_dph_cm = hypot((_dth_cm-TMath::Pi())/dth_sigma, (_dph_cm-TMath::Pi())/dph_sigma);
      if (_dth_dph_cm<_dth_dph_cm_min) {
	_dth_dph_cm_min = _dth_dph_cm;
	jgg_most_btb = jgg;
	jep_most_btb = jep;
      }
      if (_dmtot<_dmtot_min) {
	_dmtot_min = _dmtot;
	jgg_dmtot_best = jgg;
	jep_dmtot_best = jep;
      }
    }
  }

  // if no pairs, nothing to do
  if (jgg_most_btb<0||jep_most_btb<0) return;

  rcl[iep_excl].Append(rcl[iep_asso_all][jep_most_btb]);
  rcl[gg_excl].Append(rcl[gg_sel][jgg_most_btb]);

  fill_dth_dph_cm(rcl[iep_excl],rcl[gg_excl], hcmoa);
  fill_mtot(rcl[iep_excl],rcl[gg_excl], hmtot);
  fill_mmiss(rcl[iep_excl],rcl[gg_excl], hmmiss, hmmiss2);
  fill_mmiss_jpsi(rcl[iep_excl], hmmiss_jpsi, hmmiss2_jpsi);
  fill_pair_mass(rcl[iep_excl], hmep[4]);
  fill_pair_mass(rcl[gg_excl], hmgg[4]);
  fill_mctruth(rcl[iep_excl], rcl[gg_excl], 4);
  fill_count_hists(gg_excl,iep_excl,5);

}

void AnaTdav2::kin_fit() {
  assert(rcl[iep_excl].GetLength()<=1 and
	 rcl[gg_excl].GetLength()<=1);
  if (rcl[iep_excl].GetLength()==1 and
      rcl[gg_excl].GetLength()==1) {

    PndKinFitter fitter(rcl[iep_excl][0]);
    fitter.AddMassConstraint(3.096);
    fitter.Fit();

    hmep_mconst->Fill(rcl[iep_excl][0]->GetFit()->M());
    double _mtot = (rcl[iep_excl][0]->GetFit()->P4()+rcl[gg_excl][0]->P4()).M();

    hmtot_mconst->Fill( _mtot);
    double _dph, _dth;
    dth_dph_cm(rcl[gg_excl][0],rcl[iep_excl][0],_dth,_dph);
    hcmoa_mconst->Fill(_dth,_dph);

    if (fabs(_dth-TMath::Pi())<TMath::Pi()*20./180.) {
      fill_pair_mass(rcl[iep_excl], hmep[5]);
      fill_pair_mass(rcl[gg_excl], hmgg[5]);
      fill_mctruth(rcl[iep_excl], rcl[gg_excl], 5);
      hmtot_mconst_cut->Fill(_mtot);
      hcmoa_mconst_cut->Fill(_dth,_dph);
    }
  }
}

double AnaTdav2::err_mom_sq(RhoCandidate *rr) {
  TLorentzVector prec = rr->P4();
  RhoError err = rr->P4Err();
  double err_mom_fc_sq =
    prec.X()*prec.X()*err[0][0] +
    prec.Y()*prec.Y()*err[1][1] +
    prec.Z()*prec.Z()*err[2][2] +
    2*prec.X()*prec.Y()*err[0][1] +
    2*prec.X()*prec.Z()*err[0][2] +
    2*prec.Y()*prec.Z()*err[1][2];
  return err_mom_fc_sq/prec.Vect().Mag2();
}

double AnaTdav2::mom_pull_r(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->P4().Vect().Mag()-mc->P4().Vect().Mag())/sqrt(err_mom_sq(rf));
}

double AnaTdav2::mom_pull_f(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->GetFit()->P4().Vect().Mag()-mc->P4().Vect().Mag())/sqrt(err_mom_sq(rf));
}

double AnaTdav2::px_pull_r(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->P4().X()-mc->P4().X())/sqrt(rf->P4Err()[0][0]);
}

double AnaTdav2::px_pull_f(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->GetFit()->P4().X()-mc->P4().X())/sqrt(rf->P4Err()[0][0]);
}

double AnaTdav2::py_pull_r(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->P4().Y()-mc->P4().Y())/sqrt(rf->P4Err()[1][1]);
}

double AnaTdav2::py_pull_f(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->GetFit()->P4().Y()-mc->P4().Y())/sqrt(rf->P4Err()[1][1]);
}

double AnaTdav2::pz_pull_r(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->P4().Z()-mc->P4().Z())/sqrt(rf->P4Err()[2][2]);
}

double AnaTdav2::pz_pull_f(RhoCandidate* rf, RhoCandidate *mc) {
  return (rf->GetFit()->P4().Z()-mc->P4().Z())/sqrt(rf->P4Err()[2][2]);
}

void AnaTdav2::kin_fit_4c() {
  assert(rcl[iep_excl].GetLength()<=1 and
	 rcl[gg_excl].GetLength()<=1);
  if (rcl[iep_excl].GetLength()==1 and
      rcl[gg_excl].GetLength()==1) {

    RhoCandList pi0jpsi;
    pi0jpsi.Combine(rcl[iep_excl],rcl[gg_excl]);

    PndKinFitter fitter(pi0jpsi[0]);
    fitter.Add4MomConstraint(p4sys);
    fitter.Fit();

    double sig_chi2_4c = fitter.GetChi2();	// get chi2 of fit
    double sig_prob_4c = fitter.GetProb();	// access probability of fit
    double sig_pull_4c = fitter.GetPull();	// pull ?
    double sig_chi2diff4c = fitter.Chi2Diff();	// chi2diff ?
    //int ndf_4c = fitter.GetNdf();
    //cout << "sig_chi2_4c/NDF = " << sig_chi2_4c << "/" << ndf_4c << " sig_prob_4c = " << sig_prob_4c << endl;

    hpi0jpsi_chi24c->Fill(sig_chi2_4c);
    hpi0jpsi_chi24c_c->Fill(sig_chi2_4c);
    hpi0jpsi_prob4c->Fill(sig_prob_4c);
    hpi0jpsi_pull4c->Fill(sig_pull_4c);
    hpi0jpsi_chi2diff4c->Fill(sig_chi2diff4c);

    // Fill correlation histograms between chi2 of 4c and Mtot, cm opening angle
    double raw_mtot = (rcl[iep_excl][0]->P4()+rcl[gg_excl][0]->P4()).M();
    double fit_mtot = (rcl[iep_excl][0]->GetFit()->P4()+rcl[gg_excl][0]->GetFit()->P4()).M();
    double fit_mtot2 = pi0jpsi[0]->GetFit()->P4().M();

    double raw_cm_dph, raw_cm_dth;
    double fit_cm_dph, fit_cm_dth;
    dth_dph_cm(rcl[gg_excl][0],rcl[iep_excl][0],raw_cm_dth,raw_cm_dph);
    dth_dph_cm(rcl[gg_excl][0]->GetFit(),rcl[iep_excl][0]->GetFit(),fit_cm_dth,fit_cm_dph);
    hpi0jpsi_chi24c_vs_mtot_r->Fill(sig_chi2_4c,raw_mtot);
    hpi0jpsi_chi24c_vs_cm_dth_r->Fill(sig_chi2_4c,raw_cm_dth);
    hpi0jpsi_chi24c_vs_cm_dph_r->Fill(sig_chi2_4c,raw_cm_dph);
    hpi0jpsi_chi24c_vs_mtot_f->Fill(sig_chi2_4c,fit_mtot);
    hpi0jpsi_chi24c_vs_cm_dth_f->Fill(sig_chi2_4c,fit_cm_dth);
    hpi0jpsi_chi24c_vs_cm_dph_f->Fill(sig_chi2_4c,fit_cm_dph);

    int idp_r = rcl[iep_excl][0]->Daughter(0)->Charge()>0?0:1;
    int idm_r = rcl[iep_excl][0]->Daughter(0)->Charge()<0?0:1;
    assert(idp_r!=idm_r);
    int id_jpsi = (pi0jpsi[0]->Daughter(0)->Daughter(0)->PdgCode()!=22) ? 0:1;
    int idp_f = pi0jpsi[0]->Daughter(id_jpsi)->Daughter(0)->Charge()>0?0:1;
    int idm_f = pi0jpsi[0]->Daughter(id_jpsi)->Daughter(0)->Charge()<0?0:1;
    assert(idp_f!=idm_f);

    hmom_pull_ep_r->Fill(mom_pull_r(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hmom_pull_ep_f->Fill(mom_pull_f(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hmom_pull_em_r->Fill(mom_pull_r(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));
    hmom_pull_em_f->Fill(mom_pull_f(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));
    hpx_pull_ep_r->Fill(px_pull_r(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hpx_pull_ep_f->Fill(px_pull_f(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hpx_pull_em_r->Fill(px_pull_r(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));
    hpx_pull_em_f->Fill(px_pull_f(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));
    hpy_pull_ep_r->Fill(py_pull_r(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hpy_pull_ep_f->Fill(py_pull_f(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hpy_pull_em_r->Fill(py_pull_r(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));
    hpy_pull_em_f->Fill(py_pull_f(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));
    hpz_pull_ep_r->Fill(pz_pull_r(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hpz_pull_ep_f->Fill(pz_pull_f(rcl[iep_excl][0]->Daughter(idp_r),mcList[mc_posit]));
    hpz_pull_em_r->Fill(pz_pull_r(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));
    hpz_pull_em_f->Fill(pz_pull_f(rcl[iep_excl][0]->Daughter(idm_r),mcList[mc_elec]));

    double _t_gg = t_gg(rcl[gg_excl][0]);
    double _u_gg = u_gg(rcl[gg_excl][0]);
    bool t_ok = (tmin[iplab] < _t_gg && _t_gg < tmax[iplab]);
    bool u_ok = (tmin[iplab] < _u_gg && _u_gg < tmax[iplab]);
    bool _valid = t_ok || u_ok;

    fill_pair_mass(rcl[iep_excl], hmep[5]);
    fill_pair_mass(rcl[gg_excl], hmgg[5]);
    if (_valid) fill_pair_mass(rcl[iep_excl], hmep_valid[5]);
    if (_valid) fill_pair_mass(rcl[gg_excl], hmgg_valid[5]);
    fill_mctruth(rcl[iep_excl], rcl[gg_excl], 5);

    // loop over all gg pair candiates, and check if there is
    // a pair that has gives better pi0pi0jpsi than the pi0jpsi hypthesis
    RhoCandList pi0pi0jpsi;
    RhoCandList gg_mcut;
    // First narrow down the pi0 candiate list using just mass cut (no Energy vs OA cut)
    // for better acceptance of multi pi0 background events
    for (int igg=0; igg<rcl[gg].GetLength(); ++igg) {
      if (rcl[gg][igg] == rcl[gg_excl][0]) continue;
      if (0.1 < rcl[gg][igg]->M() && rcl[gg][igg]->M() < 0.17)
	gg_mcut.Append(rcl[gg][igg]);
    }

    pi0pi0jpsi.Combine(pi0jpsi,gg_mcut);
    double bg_chi2_4c = 1e9;
    bool xtra_pi0_found = false;
    for (int ibg=0; ibg<pi0pi0jpsi.GetLength(); ++ibg) {
      PndKinFitter bg_fitter(pi0pi0jpsi[ibg]);
      bg_fitter.Add4MomConstraint(p4sys);
      bg_fitter.Fit();
      if ( bg_fitter.GetChi2() < bg_chi2_4c) {
	bg_chi2_4c = bg_fitter.GetChi2();
	xtra_pi0_found = true;
      }
    }

    if (xtra_pi0_found) {
      hpi0pi0jpsi_chi24c->Fill(bg_chi2_4c);
      hpi0pi0jpsi_chi24c_c->Fill(bg_chi2_4c);
      hpi0vs2pi0_chi24c->Fill(sig_chi2_4c, bg_chi2_4c);
      hpi0vs2pi0_chi24c_c->Fill(sig_chi2_4c, bg_chi2_4c);
    }

    // if bg_chi2_4c = 1e9 at this point, it means there was no other pi0
    // candidate to test bg hypthesis. This automatically qualifies the event
    // as signal if chi2 passes the chi2 cut. If bg_chi2_4c < 1e9, then it means
    // another pi0 in the event was xtra_pi0_found, in this case, reject the event if the
    // chi2 of multi pi0 background hypothesis is better than signal hypothesis
    if (sig_chi2_4c<chi2_cut[iplab]) {
      rcl[iep_kinc].Append(rcl[iep_excl][0]);
      rcl[gg_kinc].Append(rcl[gg_excl][0]);
      fill_pair_mass(rcl[iep_excl], hmep[6]);
      fill_pair_mass(rcl[gg_excl], hmgg[6]);
      if (_valid) fill_pair_mass(rcl[iep_excl], hmep_valid[6]);
      if (_valid) fill_pair_mass(rcl[gg_excl], hmgg_valid[6]);
      fill_mctruth(rcl[iep_excl], rcl[gg_excl], 6);
      if (!xtra_pi0_found || (xtra_pi0_found && bg_chi2_4c > sig_chi2_4c)) {
    	rcl[iep_kinc_bg].Append(rcl[iep_excl][0]);
    	rcl[gg_kinc_bg].Append(rcl[gg_excl][0]);
    	fill_pair_mass(rcl[iep_excl], hmep[7]);
    	fill_pair_mass(rcl[gg_excl], hmgg[7]);
	if (_valid) fill_pair_mass(rcl[iep_excl], hmep_valid[7]);
	if (_valid) fill_pair_mass(rcl[gg_excl], hmgg_valid[7]);
	fill_mctruth(rcl[iep_excl], rcl[gg_excl], 7);
    	if (ng20mev<4) {
    	  rcl[iep_ngcut].Append(rcl[iep_excl][0]);
    	  rcl[gg_ngcut].Append(rcl[gg_excl][0]);
    	  fill_pair_mass(rcl[iep_excl], hmep[8]);
    	  fill_pair_mass(rcl[gg_excl], hmgg[8]);
	  if (_valid) fill_pair_mass(rcl[iep_excl], hmep_valid[8]);
	  if (_valid) fill_pair_mass(rcl[gg_excl], hmgg_valid[8]);
	  fill_mctruth(rcl[iep_excl], rcl[gg_excl], 8);
    	}
      }
    }

    //// Put everything in
    //if (ng20mev<3) {
    //  rcl[iep_ngcut].Append(rcl[iep_excl][0]);
    //  rcl[gg_ngcut].Append(rcl[gg_excl][0]);
    //  fill_pair_mass(rcl[iep_excl], hmep[6]);
    //  if (sig_chi2_4c<chi2_cut[iplab]) {
    //	rcl[iep_kinc].Append(rcl[iep_excl][0]);
    //	rcl[gg_kinc].Append(rcl[gg_excl][0]);
    //	fill_pair_mass(rcl[iep_excl], hmep[7]);
    //  }
    //}

  }

}

double AnaTdav2::pi0cost_cm(RhoCandidate* pi0) {
  TLorentzVector p4gg(pi0->P4());
  p4gg.Boost(boost_to_cm);
  return p4gg.CosTheta();
}

double AnaTdav2::pi0theta_cm(RhoCandidate* pi0) {
  TLorentzVector p4gg(pi0->P4());
  p4gg.Boost(boost_to_cm);
  return p4gg.Theta();
}

int AnaTdav2::find_bin(double val, const vector<double> &binning) {
  for (int ii = 0; ii < binning.size()-1; ++ii) {
    if (binning[ii] < val and
	val < binning[ii+1]) {
      return ii;
    }
  }
  return -1; // Crash
}

void AnaTdav2::fill_bins_kinc() {
  fill_bins(rcl[iep_kinc], rcl[gg_kinc]);
}

void AnaTdav2::fill_bins_kinc_bg() {
  fill_bins(rcl[iep_kinc_bg], rcl[gg_kinc_bg]);
}

void AnaTdav2::fill_bins_ngcut() {
  fill_bins(rcl[iep_ngcut], rcl[gg_ngcut]);
}

void AnaTdav2::fill_bins_excl() {
  fill_bins(rcl[iep_excl], rcl[gg_excl]);
}

TLorentzVector AnaTdav2::boost_transf(const TLorentzVector& vect_in, const TVector3& boost) {
  TLorentzVector vect_out(vect_in);
  vect_out.Boost(boost);
  return vect_out;
}
double AnaTdav2::cost_b(const TLorentzVector& v, const TVector3& boost) { return boost_transf(v,boost).CosTheta(); }
double AnaTdav2::the_b(const TLorentzVector& v, const TVector3& boost) { return TMath::RadToDeg()*(boost_transf(v,boost).Vect().Theta()); }
TLorentzVector AnaTdav2::get_p4ep(RhoCandidate* _epem) { return _epem->Daughter(0)->Charge()>0? _epem->Daughter(0)->P4(): _epem->Daughter(1)->P4(); }
TLorentzVector AnaTdav2::get_p4em(RhoCandidate* _epem) { return _epem->Daughter(0)->Charge()<0? _epem->Daughter(0)->P4(): _epem->Daughter(1)->P4(); }
int AnaTdav2::comb_bins(int nx, int ix, int iy) { return iy*nx + ix; }

void AnaTdav2::fill_bins(RhoCandList& rclep, RhoCandList& rclgg) {
  assert(rclep.GetLength()<=1 and
	 rclgg.GetLength()<=1);
  if (rclep.GetLength()==1 and
      rclgg.GetLength()==1) {

    double trecgg = t_gg(rclgg[0]);
    double urecgg = u_gg(rclgg[0]);
    double trecep = t_ep(rclep[0]);
    double urecep = u_ep(rclep[0]);
    double pi0cost_cm_rec = pi0cost_cm(rclgg[0]);
    double pi0theta_rec = rclgg[0]->P4().Theta();

    hpi0th->Fill(pi0theta_rec);
    hpi0cost_cm->Fill(pi0cost_cm_rec);
    htrecgg->Fill(trecgg);
    hurecgg->Fill(urecgg);
    htrecep->Fill(trecep);
    hurecep->Fill(urecep);

    if (rclep[0]->M()>jpsi_m_3sig_min&&rclep[0]->M()<jpsi_m_3sig_max) {
      htrecgg_mcut->Fill(trecgg);
      hurecgg_mcut->Fill(urecgg);
      htrecep_mcut->Fill(trecep);
      hurecep_mcut->Fill(urecep);
      hpi0th_mcut->Fill(pi0theta_rec);
      hpi0cost_cm_mcut->Fill(pi0cost_cm_rec);
    }

    if (event_pi0costh_cm>0) {
      htresgg->Fill(trecgg-event_t);
      htresep->Fill(trecep-event_t);
    } else {
      huresgg->Fill(urecgg-event_u);
      huresep->Fill(urecep-event_u);
    }


    int itbin = find_bin(trecgg,tu_binning);
    int iubin = find_bin(urecgg,tu_binning);
    if (itbin>=0)fill_pair_mass(rclep, hmept[itbin]);
    if (iubin>=0)fill_pair_mass(rclep, hmepu[iubin]);

    int ibin_pi0th = find_bin(pi0theta_rec,pi0th_binning);
    int ibin_pi0cost_cm = find_bin(pi0cost_cm_rec,pi0cost_cm_binning);
    fill_pair_mass(rclep, hmep_pi0th[ibin_pi0th]);
    fill_pair_mass(rclep, hmep_pi0cost_cm[ibin_pi0cost_cm]);


    /// below the positron angle dependnt stuff
    double epcosth_jpsi_rec = cost_b(get_p4ep(rclep[0]), -rclep[0]->P4().BoostVector());
    hepcosth_res->Fill(epcosth_jpsi_rec-event_epcosth_jpsi);
    int itbin_2d = find_bin(trecgg, tu_binning_2d);
    int iubin_2d = find_bin(urecgg, tu_binning_2d);
    int ibin_costh_2d = find_bin(epcosth_jpsi_rec, costh_binning_2d);
    int itu2d = itbin_2d>=0?itbin_2d:(iubin_2d>=0?2+iubin_2d:-1);
    if (itbin_2d>=0&&ibin_costh_2d>=0) fill_pair_mass(rclep, hmeptcth[comb_bins(tu_binning_2d.size()-1,itbin_2d,ibin_costh_2d)]);
    if (iubin_2d>=0&&ibin_costh_2d>=0) fill_pair_mass(rclep, hmepucth[comb_bins(tu_binning_2d.size()-1,iubin_2d,ibin_costh_2d)]);
    if (itbin_2d>=0&&ibin_costh_2d>=0) fill_pair_mass(rclep, hmeptcth0[comb_bins(tu_binning_2d.size()-1,itbin_2d,ibin_costh_2d)], m_epcth_wt0);
    if (iubin_2d>=0&&ibin_costh_2d>=0) fill_pair_mass(rclep, hmepucth0[comb_bins(tu_binning_2d.size()-1,iubin_2d,ibin_costh_2d)], m_epcth_wt0);
    if (itbin_2d>=0&&ibin_costh_2d>=0) fill_pair_mass(rclep, hmeptcth1[comb_bins(tu_binning_2d.size()-1,itbin_2d,ibin_costh_2d)], m_epcth_wt1);
    if (iubin_2d>=0&&ibin_costh_2d>=0) fill_pair_mass(rclep, hmepucth1[comb_bins(tu_binning_2d.size()-1,iubin_2d,ibin_costh_2d)], m_epcth_wt1);
    double ep_the_lab = get_p4ep(rclep[0]).Vect().Theta();
    double em_the_lab = get_p4em(rclep[0]).Vect().Theta();
    hepcosth_jpsi_rec_all->Fill(epcosth_jpsi_rec);
    hepcosth_jpsi_vs_epthlab_rec_all->Fill(epcosth_jpsi_rec,TMath::RadToDeg()*ep_the_lab);
    hepcosth_jpsi_vs_emthlab_rec_all->Fill(epcosth_jpsi_rec,TMath::RadToDeg()*em_the_lab);
    if (itu2d>=0){
      hepcosth_jpsi_rec[itu2d]->Fill(epcosth_jpsi_rec);
      hepcosth_jpsi_vs_epthlab_rec[itu2d]->Fill(epcosth_jpsi_rec,TMath::RadToDeg()*ep_the_lab);
      hepcosth_jpsi_vs_emthlab_rec[itu2d]->Fill(epcosth_jpsi_rec,TMath::RadToDeg()*em_the_lab);
    }

    // This is a bit crazy but replicate the analysis with fitted pi0 mom instead of reco jpsi mom to boost e+ angle
    //TLorentzVector _p4piz = -rclgg[0]->GetFit()->P4();
    TLorentzVector _p4jpsi = rcl[iep_excl][0]->GetFit()->P4();
    //double f_epcosth_jpsi_rec = cost_b(get_p4ep(rclep[0]), -_p4jpsi.BoostVector());
    double f_epcosth_jpsi_rec = cost_b(get_p4ep(rcl[iep_excl][0]->GetFit()), -_p4jpsi.BoostVector());
    f_hepcosth_res->Fill(f_epcosth_jpsi_rec-event_epcosth_jpsi);
    int f_ibin_costh_2d = find_bin(f_epcosth_jpsi_rec, costh_binning_2d);
    if (itbin_2d>=0&&f_ibin_costh_2d>=0) fill_pair_mass(rclep, f_hmeptcth[comb_bins(tu_binning_2d.size()-1,itbin_2d,f_ibin_costh_2d)]);
    if (iubin_2d>=0&&f_ibin_costh_2d>=0) fill_pair_mass(rclep, f_hmepucth[comb_bins(tu_binning_2d.size()-1,iubin_2d,f_ibin_costh_2d)]);
    if (itbin_2d>=0&&f_ibin_costh_2d>=0) fill_pair_mass(rclep, f_hmeptcth0[comb_bins(tu_binning_2d.size()-1,itbin_2d,f_ibin_costh_2d)], m_epcth_wt0);
    if (iubin_2d>=0&&f_ibin_costh_2d>=0) fill_pair_mass(rclep, f_hmepucth0[comb_bins(tu_binning_2d.size()-1,iubin_2d,f_ibin_costh_2d)], m_epcth_wt0);
    if (itbin_2d>=0&&f_ibin_costh_2d>=0) fill_pair_mass(rclep, f_hmeptcth1[comb_bins(tu_binning_2d.size()-1,itbin_2d,f_ibin_costh_2d)], m_epcth_wt1);
    if (iubin_2d>=0&&f_ibin_costh_2d>=0) fill_pair_mass(rclep, f_hmepucth1[comb_bins(tu_binning_2d.size()-1,iubin_2d,f_ibin_costh_2d)], m_epcth_wt1);
    //double ep_the_lab = get_p4ep(rclep[0]).Vect().Theta();
    //double em_the_lab = get_p4em(rclep[0]).Vect().Theta();
    f_hepcosth_jpsi_rec_all->Fill(f_epcosth_jpsi_rec);
    f_hepcosth_jpsi_vs_epthlab_rec_all->Fill(f_epcosth_jpsi_rec,TMath::RadToDeg()*ep_the_lab);
    f_hepcosth_jpsi_vs_emthlab_rec_all->Fill(f_epcosth_jpsi_rec,TMath::RadToDeg()*em_the_lab);
    if (itu2d>=0){
      f_hepcosth_jpsi_rec[itu2d]->Fill(f_epcosth_jpsi_rec);
      f_hepcosth_jpsi_vs_epthlab_rec[itu2d]->Fill(f_epcosth_jpsi_rec,TMath::RadToDeg()*ep_the_lab);
      f_hepcosth_jpsi_vs_emthlab_rec[itu2d]->Fill(f_epcosth_jpsi_rec,TMath::RadToDeg()*em_the_lab);
    }

  }
}

void AnaTdav2::Exec(Option_t* opt) {
  if (verb>1 or nevt%100==0)
    cout << "======== AnaTdav2::Exec evt " << nevt << " ======== " << endl;
  fAna->GetEvent();
  nevt++;
  fill_lists();
  hnevt->Fill(0);
  if (!calc_true_tu()) return;
  hnevt->Fill(1);
  //  calc_true_tu();
  calc_evt_wt();

  nocut_ref();

  bool old = false;
  if (old) {
    ep_uniq();
    pi0_sel();
    ep_pi0_asso();
    kin_excl();
    kin_fit();
    fill_bins_excl();
  } else {
    ep_all();
    pi0_sel();
    ep_pi0_asso_all();
    kin_excl_all();
    kin_fit_4c();
    //fill_bins_kinc_bg();
    fill_bins_ngcut();
    //fill_bins_kinc();
  }

}

void AnaTdav2::FinishTask() {
  cout << "AnaTdav2::FinishTask" << endl;
  cout << "Total #Evt= " << nevt << endl;
  fAna->Reset();
  write_hists();
}

void AnaTdav2::write_hists() {
  const char *root_dir = gDirectory->GetPath();

  TVectorD v(tu_binning.size());
  for (int ib=0; ib < tu_binning.size(); ++ib) {
    v[ib] = tu_binning[ib];
  }
  v.Write("tu_binning");

  TVectorD vv(tu_binning.size());
  for (int ib=0; ib < tu_binning.size(); ++ib) {
    vv[ib] = pi0th_binning[ib];
  }
  vv.Write("pi0th_binning");

  TVectorD vvv(tu_binning.size());
  for (int ib=0; ib < tu_binning.size(); ++ib) {
    vvv[ib] = pi0cost_cm_binning[ib];
  }
  vvv.Write("pi0cost_cm_binning");

  hnevt->Write();
  hwt->Write();

  pi_eff->Write();
  heff_epm->Write();

  for (int is = 0; is< nstep; ++is) {
    hmep[is]->Write();
    hmgg[is]->Write();
    hmep_valid[is]->Write();
    hmgg_valid[is]->Write();
    hnep[is]->Write();
    hngg[is]->Write();
    hnpi0jpsi[is]->Write();
  }

  hmtot->Write();
  hmmiss->Write();
  hmmiss2->Write();
  hcmoa->Write();
  hmep_mconst->Write();
  hmtot_mconst->Write();
  hcmoa_mconst->Write();
  hmtot_mconst_cut->Write();
  hcmoa_mconst_cut->Write();

  hpi0jpsi_chi24c->Write();
  hpi0jpsi_chi24c_c->Write();
  hpi0jpsi_prob4c->Write();
  hpi0jpsi_pull4c->Write();
  hpi0jpsi_chi2diff4c->Write();

  hpi0jpsi_chi24c_vs_mtot_r->Write();
  hpi0jpsi_chi24c_vs_cm_dth_r->Write();
  hpi0jpsi_chi24c_vs_cm_dph_r->Write();
  hpi0jpsi_chi24c_vs_mtot_f->Write();
  hpi0jpsi_chi24c_vs_cm_dth_f->Write();
  hpi0jpsi_chi24c_vs_cm_dph_f->Write();

  hpi0pi0jpsi_chi24c->Write();
  hpi0pi0jpsi_chi24c_c->Write();
  hpi0vs2pi0_chi24c->Write();
  hpi0vs2pi0_chi24c_c->Write();

  hng->Write();
  hng20mev->Write();
  hnch->Write();

  gDirectory->mkdir("epeff");
  gDirectory->cd("epeff");

  hepcosth_res->Write();
  f_hepcosth_res->Write();

  hepcosth_jpsi_mc_all->Write();
  hepcosth_jpsi_mc_all_wt0->Write();
  hepcosth_jpsi_mc_all_wt1->Write();

  hepcosth_jpsi_vs_epthlab_mc_all->Write();
  hepcosth_jpsi_vs_emthlab_mc_all->Write();

  hepcosth_jpsi_rec_all->Write();
  hepcosth_jpsi_vs_epthlab_rec_all->Write();
  hepcosth_jpsi_vs_emthlab_rec_all->Write();

  f_hepcosth_jpsi_rec_all->Write();
  f_hepcosth_jpsi_vs_epthlab_rec_all->Write();
  f_hepcosth_jpsi_vs_emthlab_rec_all->Write();

  for (int ii=0; ii < 4; ++ii) {
    hepcosth_jpsi_rec[ii]->Write();
    hepcosth_jpsi_vs_epthlab_rec[ii]->Write();
    hepcosth_jpsi_vs_emthlab_rec[ii]->Write();

    f_hepcosth_jpsi_rec[ii]->Write();
    f_hepcosth_jpsi_vs_epthlab_rec[ii]->Write();
    f_hepcosth_jpsi_vs_emthlab_rec[ii]->Write();

    hepcosth_jpsi_mc[ii]->Write();
    hepcosth_jpsi_vs_epthlab_mc[ii]->Write();
    hepcosth_jpsi_vs_emthlab_mc[ii]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("tcosthbins");
  gDirectory->cd("tcosthbins");
  for (int ib=0; ib < hmeptcth.size(); ++ib) {
    hmeptcth[ib]->Write();
    hmeptcth0[ib]->Write();
    hmeptcth1[ib]->Write();
    f_hmeptcth[ib]->Write();
    f_hmeptcth0[ib]->Write();
    f_hmeptcth1[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("ucosthbins");
  gDirectory->cd("ucosthbins");
  for (int ib=0; ib < hmepucth.size(); ++ib) {
    hmepucth[ib]->Write();
    hmepucth0[ib]->Write();
    hmepucth1[ib]->Write();
    f_hmepucth[ib]->Write();
    f_hmepucth0[ib]->Write();
    f_hmepucth1[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("pi0cost_cm_bins");
  gDirectory->cd("pi0cost_cm_bins");
  hpi0cost_cm->Write();
  hpi0cost_cm_mcut->Write();
  for (int ib=0; ib<tu_binning.size()-1; ++ib) {
    hmep_pi0cost_cm[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("pi0th_bins");
  gDirectory->cd("pi0th_bins");
  hpi0th->Write();
  hpi0th_mcut->Write();
  for (int ib=0; ib<tu_binning.size()-1; ++ib) {
    hmep_pi0th[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("tu_bins");
  gDirectory->cd("tu_bins");
  hpi0th->Write();
  for (int ib=0; ib<tu_binning.size()-1; ++ib) {
    hmept[ib]->Write();
    hmepu[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("mct");
  gDirectory->cd("mct");
  for (int is=0; is < nstep; ++is) {
    hmep_mct[is]->Write();
    hmep_non_mct[is]->Write();
    hthe_ep_mct[is]->Write();
    hthe_ep_mct_fwd[is]->Write();
    hthe_ep_mct_bwd[is]->Write();
    hmgg_mct[is]->Write();
    hmgg_non_mct[is]->Write();
    hthe_gg_mct[is]->Write();
    hthe_gg_mct_fwd[is]->Write();
    hthe_gg_mct_bwd[is]->Write();
    hoa_gg_mct[is]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("tu");
  gDirectory->cd("tu");
  htrecgg->Write();
  hurecgg->Write();
  htrecep->Write();
  hurecep->Write();
  httrumc->Write();
  hutrumc->Write();
  htrecgg_mcut->Write();
  hurecgg_mcut->Write();
  htrecep_mcut->Write();
  hurecep_mcut->Write();
  httrumc_vb->Write();
  hutrumc_vb->Write();
  htresgg->Write();
  huresgg->Write();
  htresep->Write();
  huresep->Write();
  htrupi0thcm->Write();
  htrupi0costhcm->Write();
  htrupi0thlab->Write();
  htrupi0thcm_mcut->Write();
  htrupi0costhcm_mcut->Write();
  htrupi0thlab_mcut->Write();
  htrupi0thcm_tcut->Write();
  htrupi0costhcm_tcut->Write();
  htrupi0thlab_tcut->Write();
  htrupi0thcm_tcut_mcut->Write();
  htrupi0costhcm_tcut_mcut->Write();
  htrupi0thlab_tcut_mcut->Write();
  htrupi0thcm_vs_m->Write();
  htrupi0costhcm_vs_m->Write();
  htrupi0thlab_vs_m->Write();
  htrupi0thcm_mcut_vs_m->Write();
  htrupi0costhcm_mcut_vs_m->Write();
  htrupi0thlab_mcut_vs_m->Write();
  gDirectory->cd(root_dir);

  gDirectory->mkdir("pull");
  gDirectory->cd("pull");
  hmom_pull_ep_r->Write();
  hmom_pull_ep_f->Write();
  hmom_pull_em_r->Write();
  hmom_pull_em_f->Write();
  hpx_pull_ep_r->Write();
  hpx_pull_ep_f->Write();
  hpx_pull_em_r->Write();
  hpx_pull_em_f->Write();
  hpy_pull_ep_r->Write();
  hpy_pull_ep_f->Write();
  hpy_pull_em_r->Write();
  hpy_pull_em_f->Write();
  hpz_pull_ep_r->Write();
  hpz_pull_ep_f->Write();
  hpz_pull_em_r->Write();
  hpz_pull_em_f->Write();
  gDirectory->cd(root_dir);

}

double AnaTdav2::dph_cm(TLorentzVector p4gg, TLorentzVector p4pm) {
  double retval = 0, dummy = 0;
  dth_dph_cm(p4gg,p4pm,dummy,retval);
  return retval;
}

double AnaTdav2::dth_cm(TLorentzVector p4gg, TLorentzVector p4pm) {
  double retval = 0, dummy = 0;
  dth_dph_cm(p4gg,p4pm,retval,dummy);
  return retval;
}

// We want to copy TLVs here before applying boost to CM. Thta way the original vector is preserved in
// case it is needed in calling function for further processing
void AnaTdav2::dth_dph_cm(TLorentzVector p4gg, TLorentzVector p4pm, double &_dth, double &_dph) {
  p4gg.Boost(boost_to_cm);
  p4pm.Boost(boost_to_cm);
  _dth = fabs(p4gg.Vect().Theta() + p4pm.Vect().Theta());
  _dph = p4gg.Vect().Phi() - p4pm.Vect().Phi();
  if (_dph<0) _dph *= -1;
}

void AnaTdav2::dth_dph_cm(RhoCandidate* _gg, RhoCandidate *_epem, double &_dth, double &_dph ) {
  TLorentzVector p4gg = _gg->P4();
  TLorentzVector p4epair = _epem->P4();
  //TLorentzVector p4ep = _epem->Daughter(0)->Charge()>0? _epem->Daughter(0)->P4(): _epem->Daughter(1)->P4();
  //p4gg.Boost(boost_to_cm);
  //p4epair.Boost(boost_to_cm);
  //_dth = fabs(p4gg.Vect().Theta() + p4epair.Vect().Theta());
  //_dph = p4gg.Vect().Phi() - p4epair.Vect().Phi();
  //if (_dph<0) _dph *= -1;
  dth_dph_cm(p4gg,p4epair,_dth,_dph);
}

bool AnaTdav2::bayes_pid(RhoCandidate* cand) {
  const int itrk = cand->GetTrackNumber();
  m_prob_drc = (PndPidProbability*) m_drc_array->At(itrk);
  m_prob_disc = (PndPidProbability*) m_disc_array->At(itrk);
  m_prob_mvd = (PndPidProbability*) m_mvd_array->At(itrk);
  m_prob_stt = (PndPidProbability*) m_stt_array->At(itrk);
  m_prob_emcb = (PndPidProbability*) m_emcb_array->At(itrk);
  double prob_comb = get_comb_prob(&PndPidProbability::GetElectronPidProb);
  return prob_comb>eid_prob_min;
}

double AnaTdav2::get_comb_prob(prob_func func) {
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

TEfficiency* AnaTdav2::rebin2d(TEfficiency *eff, int rebin) {
  assert(eff->GetDimension()==2);
  TEfficiency *tmp = (TEfficiency*) eff->Clone(Form("%s_rebinned",eff->GetName()));
  tmp->SetName(Form("%s (Rebinned, %d)", eff->GetName(), rebin));
  if (rebin<2) return eff;
  TH2F* tmpN = (TH2F*)tmp->GetPassedHistogram();
  TH2F* tmpD = (TH2F*)tmp->GetTotalHistogram();
  tmpN->RebinX(rebin);
  tmpN->RebinY(rebin);
  tmpD->RebinX(rebin);
  tmpD->RebinY(rebin);
  TEfficiency *retval = new TEfficiency(*tmpN,*tmpD);
  retval->SetTitle(tmp->GetTitle());
  return retval;
}

TH2F* AnaTdav2::smooth_hist2d(TH2F* h, int nbins= 400) {
  int nbinx = h->GetXaxis()->GetNbins();
  double xmax = h->GetXaxis()->GetBinUpEdge(nbinx);
  int nbiny = h->GetYaxis()->GetNbins();
  double ymax = h->GetYaxis()->GetBinUpEdge(nbiny);
  TH2F *smooth = new TH2F(Form("%s_smooth",h->GetName()), Form("%s (smoothed)",h->GetTitle()), nbins, 0, xmax, nbins, 0, ymax);
  cout << "xmax= " << xmax << " ymax= " << ymax << endl;
  smooth->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  smooth->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  TH1F *smoothx = (TH1F*) smooth->ProjectionX();
  TH1F *smoothy = (TH1F*) smooth->ProjectionY();
  for (int i = 1; i <= nbins; i++) {
    for (int j = 1; j <= nbins; j++) {
      double xx = smoothx->GetBinCenter(i);
      double yy = smoothy->GetBinCenter(j);
      double eff= h->Interpolate(xx,yy);
      smooth->SetBinContent(i,j,eff);
    }
  }
  return smooth;
}

TEfficiency* AnaTdav2::smooth_eff2d(TEfficiency *eff, int nbins=400)  {
  cout << "Smoothing 2d Eff" << endl;
  assert(eff->GetDimension()==2);
  TH2F* num = smooth_hist2d((TH2F*)eff->GetPassedHistogram(),nbins);
  TH2F* den = smooth_hist2d((TH2F*)eff->GetTotalHistogram(),nbins);
  TEfficiency *smooth = new TEfficiency(*num,*den);
  smooth->SetTitle(eff->GetTitle());
  cout << "Done smoothing 2d Eff" << endl;
  return smooth;
}

TH1F* AnaTdav2::smooth_hist1d(TH1F* h) {
  TH1F* hret = (TH1F*)h->Clone(Form("%s_smooth1d", h->GetName()));
  hret->SetTitle(h->GetTitle());
  hret->Smooth(10);
  return hret;
}

TEfficiency* AnaTdav2::smooth_eff1d(TEfficiency *eff)  {
  assert(eff->GetDimension()==1);
  TH1F* num = smooth_hist1d((TH1F*)eff->GetPassedHistogram());
  TH1F* den = smooth_hist1d((TH1F*)eff->GetTotalHistogram());
  TEfficiency *smooth = new TEfficiency(*num,*den);
  smooth->SetTitle(Form("%s(1d Smoothed)", eff->GetTitle()));
  return smooth;
}

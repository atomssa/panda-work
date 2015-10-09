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
  tmin{-0.092, -1.3, -2.85},
  tmax{0.59, 0.43, 0.3},
  //nevt_sim_bg{81874.0, 224120.0, 189015.0},
  nevt_sim_bg{816807.0,888292.0,898721.0},
  //nevt_xsect_bg{4.0e11, 1e11, 1e9},
  nevt_xsect_bg{4.0e11, 1e11, 2e10}, // xsect={0.2mb, 0.05mb, 0.02mb}
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
  jpsi_m_3sig_min(2.9),
  jpsi_m_3sig_max(3.3)
{
  assert(iplab>=0&&iplab<=2);
  rcl.resize(nrcl);
  gRandom->SetSeed();
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

  for (int is = 0; is< nstep; ++is) {
    hmep[is] = new TH1F(Form("hmep_%d",is),Form("hmep_%d",is),200,0,5);
    hnep[is] = new TH1F(Form("hnep_%d",is),Form("hnep_%d",is),6,-0.5,5.5);
    hngg[is] = new TH1F(Form("hngg_%d",is),Form("hngg_%d",is),26,-0.5,25.5);
    hnpi0jpsi[is] = new TH1F(Form("hnpi0jpsi_%d",is),Form("hnpi0jpsi_%d",is),6,-0.5,5.5);
  }
  hpi0th = new TH1F("hpi0th", "hpi0th", 1000, 0, TMath::Pi());
  hpi0cost_cm = new TH1F("hpi0cost_cm", "hpi0cost_cm", 1100, -1.1, 1.1);
  hpi0th_mc = new TH1F("hpi0th_mc", "hpi0th_mc", 1000, 0, TMath::Pi());
  hpi0cost_cm_mc = new TH1F("hpi0cost_cm_mc", "hpi0cost_cm_mc", 1100, -1.1, 1.1);
  for (int ib=0; ib<tu_binning.size()-1; ++ib) {
    hmep_pi0cost_cm.push_back( new TH1F(Form("hmep_pi0cost_cm_%d", ib), Form("hmep_pi0cost_cm_%d", ib), 200, 0, 5));
    hmep_pi0th.push_back( new TH1F(Form("hmep_pi0th_%d", ib), Form("hmep_pi0th_%d", ib), 200, 0, 5));
    hmept.push_back(new TH1F(Form("hmept%d", ib), Form("%3.1f < t < %3.1f;M_{inv}", tu_binning[ib], tu_binning[ib+1]), 200, 0, 5));
    hmepu.push_back(new TH1F(Form("hmepu%d", ib), Form("%3.1f < u < %3.1f;M_{inv}", tu_binning[ib], tu_binning[ib+1]), 200, 0, 5));
  }
  hmmiss = new TH1F("hmmiss", "hmmiss", 800, -4, 4);
  hmmiss2 = new TH1F("hmmiss2", "hmmiss2", 2000, -10, 10);
  hmtot = new TH1F("hmtot", "hmtot", 200, 0, 8);
  hcmoa = new TH2F("hcmoa", "hcmoa", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());

  hmep_mconst = new TH1F(Form("hmep_mconst"),Form("hmep_mconst"),200,0,5);
  hmtot_mconst = new TH1F("hmtot_mconst", "hmtot_mconst", 200, 0, 8);
  hcmoa_mconst = new TH2F("hcmoa_mconst", "hcmoa_mconst", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());
  hmtot_mconst_cut = new TH1F("hmtot_mconst_cut", "hmtot_mconst_cut", 200, 0, 8);
  hcmoa_mconst_cut = new TH2F("hcmoa_mconst_cut", "hcmoa_mconst_cut", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());

  double _tmin = iplab==0?-1.9:(iplab==1?-6.5:-14);
  if (mc_type!=1) {
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
  htrecgg_mc = new TH1F("htrecgg_mc", "htrecgg_mc", 1000, _tmin, 1);
  hurecgg_mc = new TH1F("hurecgg_mc", "hurecgg_mc", 1000, _tmin, 1);
  htrecep_mc = new TH1F("htrecep_mc", "htrecep_mc", 1000, _tmin, 1);
  hurecep_mc = new TH1F("hurecep_mc", "hurecep_mc", 1000, _tmin, 1);
  htresgg = new TH1F("htresgg", "htresgg", 1000, -3, 3);
  huresgg = new TH1F("huresgg", "huresgg", 1000, -3, 3);
  htresep = new TH1F("htresep", "htresep", 1000, -3, 3);
  huresep = new TH1F("huresep", "huresep", 1000, -3, 3);

  htrupi0thcm = new TH1F("htrupi0thcm", "htrupi0thch", 1000, 0., TMath::Pi());
  htrupi0costhcm = new TH1F("htrupi0costhcm", "htrupi0costhcm", 1100, -1.1, 1.1);
  htrupi0thlab = new TH1F("htrupi0thlab", "htrupi0thlab", 1000, 0., TMath::Pi());
  htrupi0thcm_mc = new TH1F("htrupi0thcm_mc", "htrupi0thch_mc", 1000, 0., TMath::Pi());
  htrupi0costhcm_mc = new TH1F("htrupi0costhcm_mc", "htrupi0costhcm_mc", 1100, -1.1, 1.1);
  htrupi0thlab_mc = new TH1F("htrupi0thlab_mc", "htrupi0thlab_mc", 1000, 0., TMath::Pi());
  htrupi0thcm_tc = new TH1F("htrupi0thcm_tc", "htrupi0thch_tc", 1000, 0., TMath::Pi());
  htrupi0costhcm_tc = new TH1F("htrupi0costhcm_tc", "htrupi0costhcm_tc", 1100, -1.1, 1.1);
  htrupi0thlab_tc = new TH1F("htrupi0thlab_tc", "htrupi0thlab_tc", 1000, 0., TMath::Pi());
  htrupi0thcm_tc_mc = new TH1F("htrupi0thcm_tc_mc", "htrupi0thch_tc_mc", 1000, 0., TMath::Pi());
  htrupi0costhcm_tc_mc = new TH1F("htrupi0costhcm_tc_mc", "htrupi0costhcm_tc_mc", 1100, -1.1, 1.1);
  htrupi0thlab_tc_mc = new TH1F("htrupi0thlab_tc_mc", "htrupi0thlab_tc_mc", 1000, 0., TMath::Pi());

  htrupi0thcm_vs_m = new TH2F("htrupi0thcm_vs_m", "htrupi0thch_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());
  htrupi0costhcm_vs_m = new TH2F("htrupi0costhcm_vs_m", "htrupi0costhcm_vs_m", 200, 0, 5, 1100, -1.1, 1.1);
  htrupi0thlab_vs_m = new TH2F("htrupi0thlab_vs_m", "htrupi0thlab_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());
  htrupi0thcm_mc_vs_m = new TH2F("htrupi0thcm_mc_vs_m", "htrupi0thch_mc_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());
  htrupi0costhcm_mc_vs_m = new TH2F("htrupi0costhcm_mc_vs_m", "htrupi0costhcm_mc_vs_m", 200, 0, 5, 1100, -1.1, 1.1);
  htrupi0thlab_mc_vs_m = new TH2F("htrupi0thlab_mc_vs_m", "htrupi0thlab_mc_vs_m", 200, 0, 5, 1000, 0., TMath::Pi());

  hnevt =  new TH1F("hnevt","hnevt", 10,0,10);
  hnevt->SetBinContent(3, nevt_sim_bg[iplab]);
  hnevt->SetBinContent(4, nevt_xsect_bg[iplab]);

  hpi0jpsi_chi24c = new TH1F("hpi0jpsi_chi24c","hpi0jpsi_chi24c",100,0,40);
  hpi0jpsi_prob4c  = new TH1F("hpi0jpsi_prob4c ","hpi0jpsi_prob4c ",100,0,1);

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

  //double ftmin[3] = {-0.45, -2.0, -5.0};
  //double delta[3] = {0.1, 0.2, 0.5};
  double ftmin[3] = {-0.1, -1.3, -2.8};
  double delta[3] = {0.08, 0.2, 0.4};

  cout << "double range[iplab=" << iplab << "][]= {";
  for (int i=0; i<15; ++i) {
    //double xx = ftmin[iplab] + delta[iplab]*i;
    double xx = i==0 ? tmin[iplab]: (ftmin[iplab] + delta[iplab]*i);
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

void AnaTdav2::fill_dth_dph_cm(RhoCandList& _ep, RhoCandList& _gg, TH2F* dest) {
  double _dth_cm=0, _dph_cm=0;
  for (int jep = 0; jep < _ep.GetLength(); ++jep)
    for (int jgg = 0; jgg < _gg.GetLength(); ++jgg) {
      dth_dph_cm(_gg[jgg], _ep[jep], _dth_cm, _dph_cm);
      dest->Fill(_dth_cm, _dph_cm);
    }
}

void AnaTdav2::fill_pair_mass(RhoCandList& org, TH1F* dest) {
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->M(),m_evt_wt);
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
double AnaTdav2::dist_chpi_match(RhoCandidate *rec, RhoCandidate *mc) {
  //const double dmom = (rec->P3().Mag()-mc->P3().Mag())/1.3e-2;
  const double _dth = (rec->P3().Theta()-mc->P3().Theta())/1.54e-3;
  const double _dph = (rec->P3().Phi()-mc->P3().Phi())/3.95e-3;
  //return TMath::Hypot(dmom, TMath::Hypot(dth, dph));
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
    if (mcList[j]->PdgCode()==111) {
      RhoCandidate *mcmother = mcList[j]->TheMother();
      int muid = mcmother? mcmother->GetTrackNumber(): -1;
      if ((mc_type!=1 and muid == -1) or
	  (mc_type==1 and muid == 0)) {
	event_t = t_gg(mcList[j]);
	event_u = u_gg(mcList[j]);
	// Apply "true" t and u cuts here. For the higher energies, it doesn't make
	// Sense to look at the whole range in t and u and gives wrong impression in S/B comparisons
	//if ( (tmin[iplab] < event_t and event_t < tmax[iplab]) or
	//     (tmin[iplab] < event_u and event_u < tmax[iplab]) ) {
	//  httrumc->Fill(event_t);
	//  hutrumc->Fill(event_u);
	//  htrupi0thcm->Fill(pi0theta_cm(mcList[j]));
	//  htrupi0costhcm->Fill(pi0cost_cm(mcList[j]));
	//  htrupi0thlab->Fill(mcList[j]->P4().Theta());
	//  return true;
	//} else {
	//  return false;
	//}

	httrumc->Fill(event_t);
	hutrumc->Fill(event_u);
	httrumc_vb->Fill(event_t);
	hutrumc_vb->Fill(event_u);
	htrupi0thcm->Fill(pi0theta_cm(mcList[j]));
	htrupi0costhcm->Fill(pi0cost_cm(mcList[j]));
	htrupi0thlab->Fill(mcList[j]->P4().Theta());

	bool t_ok = (tmin[iplab] < event_t && event_t < tmax[iplab]);
	bool u_ok = (tmin[iplab] < event_u && event_u < tmax[iplab]);

	if (t_ok || u_ok) {
	  htrupi0thcm_tc->Fill(pi0theta_cm(mcList[j]));
	  htrupi0costhcm_tc->Fill(pi0cost_cm(mcList[j]));
	  htrupi0thlab_tc->Fill(mcList[j]->P4().Theta());
	}

	// MC pi0 angluar distributions with jpsi mass cut on the pippim pair
	if (mc_type==0||mc_type==2) {
	  int kpip = -1, kpim=-1;
	  for (int k=0; k<mcList.GetLength(); ++k) {
	    if (mcList[k]->PdgCode()==211) { kpip = k; }
	    if (mcList[k]->PdgCode()==-211) { kpim = k; }
	    if (kpip>0&&kpim>0) {
	      double mpippim = (mcList[kpip]->P4()+mcList[kpim]->P4()).M();
	      htrupi0thcm_vs_m->Fill(mpippim, pi0theta_cm(mcList[j]));
	      htrupi0costhcm_vs_m->Fill(mpippim, pi0cost_cm(mcList[j]));
	      htrupi0thlab_vs_m->Fill(mpippim, mcList[j]->P4().Theta());
	      if ( mpippim>jpsi_m_3sig_min && mpippim<jpsi_m_3sig_max) {
		htrupi0thcm_mc->Fill(pi0theta_cm(mcList[j]));
		htrupi0costhcm_mc->Fill(pi0cost_cm(mcList[j]));
		htrupi0thlab_mc->Fill(mcList[j]->P4().Theta());

		if (t_ok || u_ok) {
		  htrupi0thcm_tc_mc->Fill(pi0theta_cm(mcList[j]));
		  htrupi0costhcm_tc_mc->Fill(pi0cost_cm(mcList[j]));
		  htrupi0thlab_tc_mc->Fill(mcList[j]->P4().Theta());
		}

		htrupi0thcm_mc_vs_m->Fill(mpippim, pi0theta_cm(mcList[j]));
		htrupi0costhcm_mc_vs_m->Fill(mpippim, pi0cost_cm(mcList[j]));
		htrupi0thlab_mc_vs_m->Fill(mpippim, mcList[j]->P4().Theta());
	      }
	      break;
	    }
	  }
	} else {
	  if (t_ok || u_ok) {
	    htrupi0thcm_tc_mc->Fill(pi0theta_cm(mcList[j]));
	    htrupi0costhcm_tc_mc->Fill(pi0cost_cm(mcList[j]));
	    htrupi0thlab_tc_mc->Fill(mcList[j]->P4().Theta());
	  }
	}

	return true;

      }
    }
  }
  cout << "MCtrue Pi0 not found. Not normal This shouldn't happen" << endl;
  return false;
}

void AnaTdav2::calc_evt_wt() {
  if (mc_type==0) {
    bool pip_found = false, pim_found = false;
    double pip_wt = 0.0, pim_wt = 0.0;
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
    m_evt_wt = nevt_xsect_bg[iplab]*m_pip_wt*m_pim_wt/nevt_sim_bg[iplab];
    //m_evt_wt = 1.0; // debug tmp
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

  if (mc_type!=1) {
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
}

void AnaTdav2::nocut_ref() {
  double tmp_eid_wt = m_evt_wt*nevt_sim_bg[iplab]/nevt_xsect_bg[iplab];
  if (mc_type==0) m_evt_wt = nevt_xsect_bg[iplab]/nevt_sim_bg[iplab];
  rcl[gg].Combine(rcl[g],rcl[g]);
  rcl[ep].Combine(rcl[e],rcl[p]);
  fill_pair_mass(rcl[ep], hmep[0]);
  fill_count_hists(gg,ep,0);
  if (mc_type==0) m_evt_wt *= tmp_eid_wt;

  rcl[iep].Combine(rcl[ie],rcl[ip]);

  fill_pair_mass(rcl[iep], hmep[1]);
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
  fill_pair_mass(rcl[iep_excl], hmep[4]);
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

  //static int noagree = 0;
  //if (jgg_most_btb!=jgg_dmtot_best||
  //    jep_most_btb!=jep_dmtot_best
  //    ) {
  //  noagree++;
  //} else {
  //  cout << "Evt" << nevt << "most_btb and best_mtot agree!" << endl;
  //}

  rcl[iep_excl].Append(rcl[iep_asso_all][jep_most_btb]);
  rcl[gg_excl].Append(rcl[gg_sel][jgg_most_btb]);

  fill_dth_dph_cm(rcl[iep_excl],rcl[gg_excl], hcmoa);
  fill_mtot(rcl[iep_excl],rcl[gg_excl], hmtot);
  fill_mmiss(rcl[iep_excl],rcl[gg_excl], hmmiss, hmmiss2);

  fill_pair_mass(rcl[iep_excl], hmep[4]);
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
      hmtot_mconst_cut->Fill(_mtot);
      hcmoa_mconst_cut->Fill(_dth,_dph);
    }
  }
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

    double chi2_4c = fitter.GetChi2();	// get chi2 of fit
    double prob_4c = fitter.GetProb();	// access probability of fit
    hpi0jpsi_chi24c->Fill(chi2_4c);
    hpi0jpsi_prob4c->Fill(prob_4c);
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

void AnaTdav2::fill_bins() {
  assert(rcl[iep_excl].GetLength()<=1 and
	 rcl[gg_excl].GetLength()<=1);
  if (rcl[iep_excl].GetLength()==1 and
      rcl[gg_excl].GetLength()==1) {

    double trecgg = t_gg(rcl[gg_excl][0]);
    double urecgg = u_gg(rcl[gg_excl][0]);
    double trecep = t_ep(rcl[iep_excl][0]);
    double urecep = u_ep(rcl[iep_excl][0]);
    double pi0cost_cm_rec = pi0cost_cm(rcl[gg_excl][0]);
    double pi0theta_rec = rcl[gg_excl][0]->P4().Theta();

    int ibin_pi0th = find_bin(pi0theta_rec,pi0th_binning);
    int ibin_pi0cost_cm = find_bin(pi0cost_cm_rec,pi0cost_cm_binning);
    int itbin = find_bin(trecgg,tu_binning);
    int iubin = find_bin(urecgg,tu_binning);

    if (itbin>=0)fill_pair_mass(rcl[iep_excl], hmept[itbin]);
    if (iubin>=0)fill_pair_mass(rcl[iep_excl], hmepu[iubin]);
    fill_pair_mass(rcl[iep_excl], hmep_pi0th[ibin_pi0th]);
    fill_pair_mass(rcl[iep_excl], hmep_pi0cost_cm[ibin_pi0cost_cm]);

    hpi0th->Fill(pi0theta_rec);
    hpi0cost_cm->Fill(pi0cost_cm_rec);
    htrecgg->Fill(trecgg);
    hurecgg->Fill(urecgg);
    htrecep->Fill(trecep);
    hurecep->Fill(urecep);

    if (rcl[iep_excl][0]->M()>jpsi_m_3sig_min&&rcl[iep_excl][0]->M()<jpsi_m_3sig_max) {
      htrecgg_mc->Fill(trecgg);
      hurecgg_mc->Fill(urecgg);
      htrecep_mc->Fill(trecep);
      hurecep_mc->Fill(urecep);
      hpi0th_mc->Fill(pi0theta_rec);
      hpi0cost_cm_mc->Fill(pi0cost_cm_rec);
    }

    htresgg->Fill(trecgg-event_t);
    htresep->Fill(trecep-event_t);
    huresgg->Fill(urecgg-event_u);
    huresep->Fill(urecep-event_u);
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
    fill_bins();
  } else {
    ep_all();
    pi0_sel();
    ep_pi0_asso_all();
    kin_excl_all();
    kin_fit_4c();
    fill_bins();
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
  hpi0jpsi_prob4c->Write();

  gDirectory->mkdir("pi0cost_cm_bins");
  gDirectory->cd("pi0cost_cm_bins");
  hpi0cost_cm->Write();
  hpi0cost_cm_mc->Write();
  for (int ib=0; ib<tu_binning.size()-1; ++ib) {
    hmep_pi0cost_cm[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("pi0th_bins");
  gDirectory->cd("pi0th_bins");
  hpi0th->Write();
  hpi0th_mc->Write();
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

  gDirectory->mkdir("tu");
  gDirectory->cd("tu");
  htrecgg->Write();
  hurecgg->Write();
  htrecep->Write();
  hurecep->Write();
  httrumc->Write();
  hutrumc->Write();
  htrecgg_mc->Write();
  hurecgg_mc->Write();
  htrecep_mc->Write();
  hurecep_mc->Write();
  httrumc_vb->Write();
  hutrumc_vb->Write();
  htresgg->Write();
  huresgg->Write();
  htresep->Write();
  huresep->Write();
  htrupi0thcm->Write();
  htrupi0costhcm->Write();
  htrupi0thlab->Write();
  htrupi0thcm_mc->Write();
  htrupi0costhcm_mc->Write();
  htrupi0thlab_mc->Write();
  htrupi0thcm_tc->Write();
  htrupi0costhcm_tc->Write();
  htrupi0thlab_tc->Write();
  htrupi0thcm_tc_mc->Write();
  htrupi0costhcm_tc_mc->Write();
  htrupi0thlab_tc_mc->Write();
  htrupi0thcm_vs_m->Write();
  htrupi0costhcm_vs_m->Write();
  htrupi0thlab_vs_m->Write();
  htrupi0thcm_mc_vs_m->Write();
  htrupi0costhcm_mc_vs_m->Write();
  htrupi0thlab_mc_vs_m->Write();



}

void AnaTdav2::dth_dph_cm(RhoCandidate* _gg, RhoCandidate *_epem, double &_dth, double &_dph ) {
  TLorentzVector p4gg = _gg->P4();
  TLorentzVector p4epair = _epem->P4();
  TLorentzVector p4ep = _epem->Daughter(0)->Charge()>0? _epem->Daughter(0)->P4(): _epem->Daughter(1)->P4();
  p4gg.Boost(boost_to_cm);
  p4epair.Boost(boost_to_cm);
  _dth = fabs(p4gg.Vect().Theta() + p4epair.Vect().Theta());
  _dph = p4gg.Vect().Phi() - p4epair.Vect().Phi();
  if (_dph<0) _dph *= -1;
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

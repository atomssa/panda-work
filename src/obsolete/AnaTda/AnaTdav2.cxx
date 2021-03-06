#include "AnaTdav2.h"

#include "FairTask.h"
#include "RhoCandList.h"
#include "PndKinFitter.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

AnaTdav2::AnaTdav2(const int& _iplab, const int& itype, const int& brem, const int& eid_meth):
  verb(false),
  nevt(0),
  brem_corr(brem!=0),
  bg_mc(itype==0),
  eid_param_method(eid_meth!=0),
  iplab((_iplab>=0&&_iplab<3)?_iplab:0),
  plab{5.513,8.,12.},
  mom_antip(plab[iplab]),
  boost_to_cm(),
  boost_to_lab(),
  p4pbar(),
  p4targ(),
  event_t(-9999.),
  event_u(-9999.),
  tmin{-0.443789, -0.5, -0.5},
  tmax{0.616486, 0.457248, 0.31538},
  eff_file_name("eff/effic_smooth.root"),
  eff_hist_name("eff_ep_em_rad"),
  eff_hist_rad(true),
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
  require_exclusivity(false)
{
  assert(iplab>=0&&iplab<=2);
  rcl.resize(nrcl);
  gRandom->SetSeed();
}

AnaTdav2::~AnaTdav2() {

}

void AnaTdav2::init_hists() {
  eff_file = TFile::Open(eff_file_name.c_str());
  heff_epm = (TH2F*) eff_file->Get(eff_hist_name.c_str())->Clone("heff_epm");
  for (int is = 0; is< nstep; ++is) {
    hmep[is] = new TH1F(Form("hmep_%d",is),Form("hmep_%d",is),200,0,5);
    hnep[is] = new TH1F(Form("hnep_%d",is),Form("hnep_%d",is),6,-0.5,5.5);
    hngg[is] = new TH1F(Form("hngg_%d",is),Form("hngg_%d",is),26,-0.5,25.5);
  }
  hpi0th = new TH1F("hpi0th", "hpi0th", 200, 0, TMath::Pi());
  hpi0cost_cm = new TH1F("hpi0cost_cm", "hpi0cost_cm", 200, -1., 1.);
  for (int ib=0; ib<nbinth; ++ib) {
    hmep_pi0cost_cm[ib] = new TH1F(Form("hmep_pi0cost_cm_%d", ib), Form("hmep_pi0cost_cm_%d", ib), 200, 0, 5);
    hmep_pi0th[ib] = new TH1F(Form("hmep_pi0th_%d", ib), Form("hmep_pi0th_%d", ib), 200, 0, 5);
    hmept[ib] = new TH1F(Form("hmept%d", ib), Form("%3.1f < t < %3.1f;M_{inv}", tu_binning[ib], tu_binning[ib+1]), 200, 0, 5);
    hmepu[ib] = new TH1F(Form("hmepu%d", ib), Form("%3.1f < u < %3.1f;M_{inv}", tu_binning[ib], tu_binning[ib+1]), 200, 0, 5);
  }
  hmtot = new TH1F("hmtot", "hmtot", 200, 0, 8);
  hcmoa = new TH2F("hcmoa", "hcmoa", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());

  hmep_mconst = new TH1F(Form("hmep_mconst"),Form("hmep_mconst"),200,0,5);
  hmtot_mconst = new TH1F("hmtot_mconst", "hmtot_mconst", 200, 0, 8);
  hcmoa_mconst = new TH2F("hcmoa_mconst", "hcmoa_mconst", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());
  hmtot_mconst_cut = new TH1F("hmtot_mconst_cut", "hmtot_mconst_cut", 200, 0, 8);
  hcmoa_mconst_cut = new TH2F("hcmoa_mconst_cut", "hcmoa_mconst_cut", 200, 0, 2*TMath::Pi(), 200, 0, 2*TMath::Pi());

  double _tmin = iplab==0?-1.9:(iplab==1?-6.5:-14);
  if (bg_mc) {
    _tmin = iplab==0?-11:(iplab==1?-16:-24);
  }
  httrumc = new TH1F("httrumc", "httrumc", 200, _tmin, 1);
  hutrumc = new TH1F("hutrumc", "hutrumc", 200, _tmin, 1);
  htrecgg = new TH1F("htrecgg", "htrecgg", 200, _tmin, 1);
  hurecgg = new TH1F("hurecgg", "hurecgg", 200, _tmin, 1);
  htrecep = new TH1F("htrecep", "htrecep", 200, _tmin, 1);
  hurecep = new TH1F("hurecep", "hurecep", 200, _tmin, 1);
  htresgg = new TH1F("htresgg", "htresgg", 200, -3, 3);
  huresgg = new TH1F("huresgg", "huresgg", 200, -3, 3);
  htresep = new TH1F("htresep", "htresep", 200, -3, 3);
  huresep = new TH1F("huresep", "huresep", 200, -3, 3);

  htrupi0thcm = new TH1F("htrupi0thcm", "htrupi0thch", 200, 0., TMath::Pi());
  htrupi0costhcm = new TH1F("htrupi0costhcm", "htrupi0costhcm", 220, -1.1, 1.1);
  htrupi0thlab = new TH1F("htrupi0thlab", "htrupi0thlab", 200, 0., TMath::Pi());
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

  // Equal subdivisions in costh_cm, boost to lab for th bins
  TLorentzVector pi0;
  TVector3 u;
  for (int ii = 0; ii < nbinth+1; ++ii) {
    pi0cost_cm_binning[ii] = -1.0 + ( 2.0*ii/nbinth );
    u.SetMagThetaPhi(1,TMath::ACos(pi0cost_cm_binning[ii]),0);
    pi0.SetPxPyPzE(u.Px(),u.Py(),u.Pz(),TMath::Hypot(0.1349766,1.0));
    pi0.Boost(boost_to_lab);
    pi0th_binning[nbinth-ii] = pi0.Theta();
  }
  pi0th_binning[nbinth]+=0.0001; // safety for numerical error
  print_binning(pi0cost_cm_binning,"pi0cost_cm_binning");
  print_binning(pi0th_binning,"pi0th_binning");

  // t and u binnings. This is the tricky part.
  // Use bins of 0.05 on either side of 0, and use partial bin for the outermost ones
  // for iplab==0 (5.513 GeV/c), go from -0.443789(@90deg) in steps of del=0.05
  // for iplab==1 (8 GeV/c), go from -0.5(@??deg)by bins of del=0.1
  // for iplab==2 (8 GeV/c), go from -0.5(@??deg)by bins of del=0.1
  const double step = iplab==0?0.1:0.1;
  const double start = iplab==0?tmin[iplab]:-0.5;
  for (int ibin=0; ibin<nbinth+1; ++ibin) {
    tu_binning[ibin] = start + (step*ibin);
  }
  print_binning(tu_binning, "tu_binning");
}

void AnaTdav2::print_binning(double *b, const char* n) {
  cout << n << ": {";
  for (int ii = 0; ii < nbinth+1; ++ii) cout << b[ii] << ", ";
  cout << "}" << endl;
}

InitStatus AnaTdav2::Init() {
  cout << "AnaTdav2::Init" << endl;
  fAna = new PndAnalysis();
  beam_cond();
  init_hists();
  return kSUCCESS;
}

void AnaTdav2::fill_mtot(RhoCandList& _ep, RhoCandList& _gg, TH1F* dest) {
  for (int jep = 0; jep < _ep.GetLength(); ++jep)
    for (int jgg = 0; jgg < _gg.GetLength(); ++jgg)
      dest->Fill(m(_gg[jgg], _ep[jep]));
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
  for (int j = 0; j < org.GetLength(); ++j) dest->Fill(org[j]->M());
}

void AnaTdav2::fill_count_hists(int _gg, int _ep, int ihist) {
  hnep[ihist]->Fill(rcl[_ep].GetLength());
  hngg[ihist]->Fill(rcl[_gg].GetLength());
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
  const double dth = (rec->P3().Theta()-mc->P3().Theta())/1.54e-3;
  const double dph = (rec->P3().Phi()-mc->P3().Phi())/3.95e-3;
  //return TMath::Hypot(dmom, TMath::Hypot(dth, dph));
  return TMath::Hypot(dth, dph);
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
  for (int ip = 0; ip < inp_pi.GetLength(); ++ip) {
    for (int im = 0; im < inm_pi.GetLength(); ++im) {
      const double dp = dist_chpi_match(inp_pi[ip],mcList[mc_p]);
      const double dm = dist_chpi_match(inm_pi[im],mcList[mc_m]);
      const double dtot = TMath::Hypot(dp,dm);
      if (dtot < dmin) {
	match_p = ip;
	match_m = im;
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
      if ((bg_mc and muid == -1) or
	  (!bg_mc and muid == 0)) {
	event_t = t_gg(mcList[j]);
	event_u = u_gg(mcList[j]);
	// Apply "true" t and u cuts here. For the higher energies, it doesn't make
	// Sense to look at the whole range in t and u and gives wrong impression in S/B comparisons
	if ( (tmin[iplab] < event_t and event_t < tmax[iplab]) or
	     (tmin[iplab] < event_u and event_u < tmax[iplab]) ) {
	  httrumc->Fill(event_t);
	  hutrumc->Fill(event_u);
	  htrupi0thcm->Fill(pi0theta_cm(mcList[j]));
	  htrupi0costhcm->Fill(pi0cost_cm(mcList[j]));
	  htrupi0thlab->Fill(mcList[j]->P4().Theta());
	  return true;
	} else {
	  return false;
	}
	//if (bg_mc) {
	//  httrumc->Fill(event_t);
	//  hutrumc->Fill(event_u);
	//  htrupi0thcm->Fill(pi0theta_cm(mcList[j]));
	//  htrupi0costhcm->Fill(pi0cost_cm(mcList[j]));
	//  htrupi0thlab->Fill(mcList[j]->P4().Theta());
	//  return true;
	//} else {
	//}
      }
    }
  }
  cout << "MCtrue Pi0 not found. Not normal This shouldn't happen" << endl;
  return false;
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

  if (bg_mc) {
    //rcl[p].Cleanup();
    //rcl[e].Cleanup();
    //fAna->FillList(rcl[p], (brem_corr?"BremPionAllPlus":"PionAllPlus"));
    //fAna->FillList(rcl[e], (brem_corr?"BremPionAllMinus":"PionAllMinus"));

    fAna->FillList(rcl[pip], "PionAllPlus");
    fAna->FillList(rcl[pim], "PionAllMinus");
    //fAna->FillList(rcl[ip], (brem_corr?"BremElectronVeryTightPlus":"ElectronAllPlus"), "PidAlgoEmcBayes");
    //fAna->FillList(rcl[ie], (brem_corr?"BremElectronVeryTightMinus":"ElectronAllMinus"), "PidAlgoEmcBayes");
    //fAna->FillList(rcl[ip], (brem_corr?"BremElectronAllPlus":"ElectronAllPlus"));
    //fAna->FillList(rcl[ie], (brem_corr?"BremElectronAllMinus":"ElectronAllMinus"));
    charged_pion_filter(rcl[ip], rcl[ie], rcl[pip], rcl[pim], rcl[p], rcl[e]);
  } else {
    eid_filter(rcl[ip],rcl[p]);
    eid_filter(rcl[ie],rcl[e]);
  }
}

void AnaTdav2::nocut_ref() {
  rcl[gg].Combine(rcl[g],rcl[g]);
  rcl[ep].Combine(rcl[e],rcl[p]);
  fill_pair_mass(rcl[ep], hmep[0]);
  fill_count_hists(gg,ep,0);
  if (bg_mc) {
    rcl[iep].Combine(rcl[e],rcl[p]);
  } else {
    rcl[iep].Combine(rcl[ie],rcl[ip]);
  }
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
  fill_pair_mass(rcl[iep_excl], hmep[4]);
  fill_count_hists(gg_excl,iep_excl,3);
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

int AnaTdav2::find_bin(double val, double *binning) {
  for (int ii = 0; ii < nbinth; ++ii) {
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

  if (!calc_true_tu()) return;

  nocut_ref();
  ep_uniq();
  pi0_sel();
  ep_pi0_asso();
  kin_excl();
  kin_fit();
  fill_bins();
}

void AnaTdav2::Finish() {
  cout << "AnaTdav2::Exec" << endl;
  fAna->Reset();
  write_hists();
}

void AnaTdav2::write_hists() {
  const char *root_dir = gDirectory->GetPath();

  for (int is = 0; is< nstep; ++is) {
    hmep[is]->Write();
    hnep[is]->Write();
    hngg[is]->Write();
  }

  hmtot->Write();
  hcmoa->Write();
  hmep_mconst->Write();
  hmtot_mconst->Write();
  hcmoa_mconst->Write();
  hmtot_mconst_cut->Write();
  hcmoa_mconst_cut->Write();

  gDirectory->mkdir("pi0cost_cm_bins");
  gDirectory->cd("pi0cost_cm_bins");
  hpi0cost_cm->Write();
  for (int ib=0; ib<nbinth; ++ib) {
    hmep_pi0cost_cm[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("pi0th_bins");
  gDirectory->cd("pi0th_bins");
  hpi0th->Write();
  for (int ib=0; ib<nbinth; ++ib) {
    hmep_pi0th[ib]->Write();
  }
  gDirectory->cd(root_dir);

  gDirectory->mkdir("tu_bins");
  gDirectory->cd("tu_bins");
  hpi0th->Write();
  for (int ib=0; ib<nbinth; ++ib) {
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
  htresgg->Write();
  huresgg->Write();
  htresep->Write();
  huresep->Write();
  htrupi0thcm->Write();
  htrupi0costhcm->Write();
  htrupi0thlab->Write();


}

void AnaTdav2::dth_dph_cm(RhoCandidate* _gg, RhoCandidate *_epem, double &dth, double &dph ) {
  TLorentzVector p4gg = _gg->P4();
  TLorentzVector p4epair = _epem->P4();
  TLorentzVector p4ep = _epem->Daughter(0)->Charge()>0? _epem->Daughter(0)->P4(): _epem->Daughter(1)->P4();
  p4gg.Boost(boost_to_cm);
  p4epair.Boost(boost_to_cm);
  dth = fabs(p4gg.Vect().Theta() + p4epair.Vect().Theta());
  dph = p4gg.Vect().Phi() - p4epair.Vect().Phi();
  if (dph<0) dph *= -1;
}

// This is going to be a long function...
bool AnaTdav2::bayes_pid(RhoCandidate* cand) {
  return true;
}

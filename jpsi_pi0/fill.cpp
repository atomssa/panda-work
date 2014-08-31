#include "filler.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TClonesArray.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>
#include <fstream>
#include <stdlib.h>

using namespace std;

static const bool verb = false;

static const double rtd= TMath::RadToDeg();
//static const double m_2pi; // = 2.0*TMath::Pi();
//static const double m_pi; // = TMath::Pi();

//static const double ct_max = 1.1;
//static const double ct_min = -1.1;
//static const double th_min = -0.1;
//static const double th_max = m_pi+0.1;
//static const double mom_min = 0;
//static const double mom_max = 6;
//static const double oa_min = -0.2;
//static const double oa_max = m_pi+0.2;
//static const double t_0 = 2;
//static const double t_1 = -2;

static int nbins = 200;
static axis mom_bins(nbins, 0.0, 6.0);
static axis ene_bins(nbins, 0.0, 6.0);
static axis t_bins(nbins, -2.0, 2.0);
static axis the_bins(nbins, -0.1, m_pi+0.1);
static axis phi_bins(nbins, -m_pi, m_pi);
static axis cost_bins(nbins, -1.1, 1.1);
static axis oa_bins(nbins, -0.1, m_pi+0.1);

//static const double mass_jpsi = 3.096;
//static const double mass_pi0 = 0.135;
static const double mass_prot= 0.938;

bool compare_unsorted(const vector<int>&, const vector<int>&);
bool compare_sorted(const vector<int>&, const vector<int>&);

void print_event(const vector<int>&);
void print_log(int, int, vector<int>&);

bool accept_pc(const TLorentzVector &v);
bool accept_pi0(const TLorentzVector &v);
bool accept_elec(const TLorentzVector &v);

void fill_bg_hists(TChain *, TClonesArray*);
void fill_sig_hists(TChain *, TClonesArray*);

void print_p4(const char*, const TLorentzVector&);
void set_init_cond(TVector3&, TLorentzVector&, TLorentzVector&, TLorentzVector&,TLorentzVector&, TLorentzVector&, TLorentzVector&);

int main(const int argc, const char **argv) {

  if (argc<3) {
    std::cout << "Need two arguments [type(0->Bg,1->Sig), file-name]" << std::endl;
    return -1;
  }

  int type = atoi(argv[1]);

  TChain *data_in = new TChain("data","data_in");
  ifstream inf;
  inf.open(argv[2]);
  string root_file;
  while(true) {
    inf >> root_file;
    if (!inf.good()) break;
    data_in->AddFile(root_file.c_str());
  }

  TClonesArray *part_array = new TClonesArray("TParticle");
  data_in->SetBranchAddress("Particles",&part_array);

  if (type == 0 ) {
    fill_bg_hists(data_in, part_array);
  } else if (type == 1) {
    fill_sig_hists(data_in, part_array);
  } else {
    cout << "type index " << type << " out of range" << endl;
    return -1;
  }

  return 0;

}

void fill_bg_hists(TChain *data_in, TClonesArray *part_array) {

  std::vector<int> ref;
  ref.push_back(-211);
  ref.push_back(111);
  ref.push_back(211);

  TVector3 boost_to_cm;
  TLorentzVector p4pbar,p4p,p4pbarp,p4pbar_cm,p4p_cm,p4pbarp_cm;
  set_init_cond(boost_to_cm, p4pbar,p4p,p4pbarp_cm,p4pbar_cm,p4p,p4pbarp_cm);

  const double p_antip = 5.513;
  const double m_max = 1.1 * TMath::Sqrt(2*mass_prot + 2*p_antip*mass_prot);
  const axis mass_bins(nbins,0,m_max);

  enum {pbar=0, pi0=1, pim=2, pip=3, pipm=4};
  const char *pp[] = {"pbar", "pi0", "pim", "pip", "pipm"};
  const char *pp_title[] = {"#bar{p}", "#pi^{0}", "#pi^{+}", "#pi^{-}", "(#pi^{+}#pi^{-})"};
  std::vector<filler*> fillers_lab;

  const char *lab_frame[] = {"lab"," (Lab frame)"};
  // (pi0) mom and angles
  // (pi+) - (pi-) system mass, mom and angles
  fillers_lab.push_back(new var1d(pi0, "mom", mom_bins, lab_frame, pp, pp_title));
  fillers_lab.push_back(new var1d(pi0, "the", the_bins, lab_frame, pp, pp_title));
  fillers_lab.push_back(new var1d(pi0, "cost", cost_bins, lab_frame, pp, pp_title));
  fillers_lab.push_back(new var1d(pi0, "phi", phi_bins, lab_frame, pp, pp_title));

  // (pi+) - (pi-) system mass, mom and angles
  fillers_lab.push_back(new var1d(pim, pip, "mom", mom_bins, lab_frame, pp, pp_title));
  fillers_lab.push_back(new var1d(pim, pip, "the", the_bins, lab_frame, pp, pp_title));
  fillers_lab.push_back(new var1d(pim, pip, "phi", phi_bins, lab_frame, pp, pp_title));
  fillers_lab.push_back(new var1d(pim, pip, "oa", oa_bins, lab_frame, pp, pp_title));
  fillers_lab.push_back(new var1d(pi0, pipm, "oa", oa_bins, lab_frame, pp, pp_title));

  std::vector<filler*> fillers_inv;
  const char *invariant[] = {"inv", ""};
  fillers_inv.push_back(new var1d(pim, pip, "mass", mass_bins, invariant, pp, pp_title));
  fillers_inv.push_back(new var1d(pbar, pipm, "u", axis(200,-11,1), invariant, pp, pp_title));
  fillers_inv.push_back(new var1d(pbar, pi0, "t", axis(200,-11,1), invariant, pp, pp_title));
  // "Dalitz" plot
  fillers_inv.push_back(new var2d(pi0, pip, "mass", mass_bins, invariant, pip, pim, "mass", mass_bins, invariant, pp, pp_title) );

  std::vector<filler*> fillers_cm;
  const char *cm_frame[] = {"cm"," (CM frame)"};
  fillers_cm.push_back(new var1d(pi0, "mom", mom_bins, cm_frame, pp, pp_title));
  fillers_cm.push_back(new var1d(pi0, "the", the_bins, cm_frame, pp, pp_title));
  fillers_cm.push_back(new var1d(pim, pip, "mom", mom_bins, cm_frame, pp, pp_title));
  fillers_cm.push_back(new var1d(pim, pip, "the", the_bins, cm_frame, pp, pp_title));
  fillers_cm.push_back(new var1d(pim, pip, "oa", oa_bins, cm_frame, pp, pp_title));
  fillers_cm.push_back(new var1d(pi0, pipm, "oa", oa_bins, cm_frame, pp, pp_title));

  const int Nevt = data_in->GetEntries();

  for (int ievt=0; ievt<Nevt; ievt++) {

    if (ievt>100000) break;
    if (ievt%10000==0)
      std::cout << "event " << ievt << "/" << Nevt << std::endl;
    data_in->GetEntry(ievt);

    // veto events with more than three particles in final state
    int npart = part_array->GetEntriesFast();
    if (npart!=3) {
      std::cout << "npart= " << npart << std::endl;
    }

    TLorentzVector p4pi0, p4pip, p4pim, p4pipm;
    int cpi0=0,cpip=0,cpim=0;
    for (int ipart=0; ipart < npart; ++ipart) {
      const TParticle* part = (TParticle*) part_array->At(ipart);
      if (verb) {
	cout << "part= " << part << endl;
	std::cout << " part " << ipart << " pdg= " << part->GetPDG()->PdgCode() << std::endl;
	part->Print();
      }
      int pdg = part->GetPDG()->PdgCode();
      // This would insure ordering
      if ( pdg == -211 ) { ((TParticle*)part_array->At(ipart))->Momentum(p4pim); cpim++; }
      if ( pdg == 111  ) { ((TParticle*)part_array->At(ipart))->Momentum(p4pi0); cpi0++; }
      if ( pdg == 211  ) { ((TParticle*)part_array->At(ipart))->Momentum(p4pip); cpip++; }
    }
    p4pipm = p4pip + p4pim;

    // check event type
    if ( cpim!=1 || cpip!=1 || cpi0!=1 ) {
      std::cout << "Event wrong type " << cpim << " pi- " << cpip << " pi+ " << cpi0 << " pi0 found where 1 expected" << std::endl;
      continue;
    }

    // check acceptance
    //if (!accept_pc(p4pim)||!accept_pc(p4pip)||accept_pi0(p4pi0)) continue;
    vector<TLorentzVector> p4s;
    p4s.push_back(p4pbar);
    p4s.push_back(p4pi0);
    p4s.push_back(p4pim);
    p4s.push_back(p4pip);
    p4s.push_back(p4pipm);

    for (auto fill: fillers_lab) (*fill)(p4s);
    for (auto fill: fillers_inv) (*fill)(p4s);
    for (auto fill: fillers_cm) (*fill)(p4s,boost_to_cm);

  }

  TFile *fout = TFile::Open("out_bg.root","RECREATE");
  fout->cd();
  for (auto fill: fillers_lab) fill->Write();
  for (auto fill: fillers_cm) fill->Write();
  for (auto fill: fillers_inv) fill->Write();
  fout->Write();
  fout->Close();

}

void fill_sig_hists(TChain *data_in, TClonesArray *part_array) {

  TVector3 boost_to_cm;
  TLorentzVector p4pbar,p4p,p4pbarp,p4pbar_cm,p4p_cm,p4pbarp_cm;
  set_init_cond(boost_to_cm,p4pbar,p4p,p4pbarp,p4pbar_cm,p4p_cm,p4pbarp_cm);
  const double p_antip = 5.513;

  const double m_max = 1.1 * TMath::Sqrt(2*mass_prot + 2*p_antip*mass_prot);
  const axis mass_bins(200,0,m_max);

  std::vector<int> ref;
  ref.push_back(-11);
  ref.push_back(11);
  ref.push_back(22);
  ref.push_back(22);

  enum {pbar=0, g1=1, g2=2, ep=3, em=4, pi0=5, jpsi=6};
  //const int ipbar=0, ig1=1, ig2=2, iep=3, iem=4, ipi0=5, ijpsi=6;
  const char *pp[] = {"pbar", "g1", "g2", "ep", "em", "pi0", "jpsi"};
  const char *pp_title[] = {"#bar{p}", "#gamma_{1}", "#gamma_{2}", "e^{+}", "e^{-}", "#pi^{0}", "J/#psi"};

  /*
  // This can simplify the filler registration by doing loops over similar fillers lists
  std::vector<std::vector<filler*>> filler_vects;
  std::vector<filler*> invariant_fillers;
  filler_vects.push_back(invariant_fillers);
  std::vector<filler*> lab_fillers;
  filler_vects.push_back(lab_fillers);
  std::vector<filler*> cm_fillers;
  filler_vects.push_back(cm_fillers);
  */

  std::vector<filler*> fillers_inv;
  const char *invariant[] = {"inv", ""};
  fillers_inv.push_back( new var1d(pi0, jpsi, "mass", mass_bins, invariant, pp, pp_title ) );
  fillers_inv.push_back( new var1d(ep, em, "mass", mass_bins, invariant, pp, pp_title ) );
  fillers_inv.push_back( new var1d(g1, g2, "mass", mass_bins, invariant, pp, pp_title ) );
  fillers_inv.push_back( new var1d(pbar, jpsi, "u", t_bins, invariant, pp, pp_title ) );
  fillers_inv.push_back( new var1d(pbar, pi0, "t", t_bins, invariant, pp, pp_title ) );

  std::vector<filler*> lab_fillers;
  const char *lab_frame[] = {"lab"," (Lab frame)"};
  for (int ipart=1; ipart<= 1 /*6*/; ++ipart) {   // for all single particle
    lab_fillers.push_back(new var1d(ipart, "mom", mom_bins, lab_frame, pp, pp_title));
    lab_fillers.push_back(new var1d(ipart, "e", ene_bins, lab_frame, pp, pp_title));
    lab_fillers.push_back(new var1d(ipart, "the", the_bins, lab_frame, pp, pp_title));
    lab_fillers.push_back(new var1d(ipart, "cost", cost_bins, lab_frame, pp, pp_title));
  }
  lab_fillers.push_back(new var1d(ep, em, "oa", oa_bins, lab_frame, pp, pp_title )); // (e+-e-) OA
  lab_fillers.push_back(new var1d(g1, g2, "oa", oa_bins, lab_frame, pp, pp_title )); // (g1-g2) OA
  lab_fillers.push_back(new var1d(pi0, jpsi, "oa", oa_bins, lab_frame, pp, pp_title )); // (pi0-jpsi) OA

  std::vector<filler*> cm_fillers;
  const char *cm_frame[] = {"cm"," (CM frame)"};
  for (int ipart=1; ipart<=1 /*6*/; ++ipart) {
    cm_fillers.push_back(new var1d(ipart, "mom", mom_bins, cm_frame, pp, pp_title));
    cm_fillers.push_back(new var1d(ipart, "e", ene_bins, cm_frame, pp, pp_title));
    cm_fillers.push_back(new var1d(ipart, "the", the_bins, cm_frame, pp, pp_title));
    cm_fillers.push_back(new var1d(ipart, "cost", cost_bins, cm_frame, pp, pp_title));
  }
  cm_fillers.push_back(new var1d(ep, em, "oa", oa_bins, cm_frame, pp, pp_title )); // (e+-e-) OA
  cm_fillers.push_back(new var1d(g1, g2, "oa", oa_bins, cm_frame, pp, pp_title )); // (g1-g2) OA
  cm_fillers.push_back(new var1d(pi0, jpsi, "oa", oa_bins, cm_frame, pp, pp_title )); // (pi0-jpsi) OA

  std::vector<filler*> jpsif_fillers;
  const char *jpsi_frame[] {"jpsiframe"," (J/#psi frame)"};
  jpsif_fillers.push_back(new var1d(ep, "mom", mom_bins, jpsi_frame, pp, pp_title)); // (e+) mom
  jpsif_fillers.push_back(new var1d(ep, "cost", cost_bins, jpsi_frame, pp, pp_title)); // (e+) cos(theta)
  jpsif_fillers.push_back(new var1d(em, "mom", mom_bins, jpsi_frame, pp, pp_title)); // (e+) mom
  jpsif_fillers.push_back(new var1d(em, "cost", cost_bins, jpsi_frame, pp, pp_title)); // (e+) cos(theta)
  jpsif_fillers.push_back(new var1d(ep, em, "oa", cost_bins, jpsi_frame, pp, pp_title)); // (e+-e-) OA

  std::vector<filler*> pi0f_fillers;
  const char *pi0_frame[] {"pi0frame"," (#pi^{0} frame)"};
  pi0f_fillers.push_back(new var1d(g1, "mom", mom_bins, pi0_frame, pp, pp_title)); // (g1) mom
  pi0f_fillers.push_back(new var1d(g1, "cost", cost_bins, pi0_frame, pp, pp_title)); // (g1) cos(theta)
  pi0f_fillers.push_back(new var1d(g2, "mom", mom_bins, pi0_frame, pp, pp_title)); // (g2) mom
  pi0f_fillers.push_back(new var1d(g2, "cost", cost_bins, pi0_frame, pp, pp_title)); // (g2) cos(theta)
  pi0f_fillers.push_back(new var1d(g1, g2, "oa", oa_bins, pi0_frame, pp, pp_title)); // (g1-g2) OA

  std::vector<filler*> var2d_fillers;
  var2d_fillers.push_back(new var2d(g1, "mom", mom_bins, lab_frame, g2, "mom", mom_bins, lab_frame, pp, pp_title) );

  const int Nevt = data_in->GetEntries();

  for (int ievt=0; ievt<Nevt; ievt++) {

    //if (ievt>10) break;
    if (ievt%10000==0)
      std::cout << "event " << ievt << "/" << Nevt << std::endl;
    data_in->GetEntry(ievt);

    // veto events with more than three particles in final state
    int npart = part_array->GetEntriesFast();
    if (npart!=4) {
      std::cout << "npart= " << npart << std::endl;
    }

    bool first = true;
    int cep=0,cem=0,cg=0;
    TLorentzVector p4g1, p4g2, p4ep, p4em, p4pi0, p4jpsi;
    for (int ipart=0; ipart < npart; ++ipart) {
      const TParticle* part = (TParticle*) part_array->At(ipart);
      if (verb) {
	cout << "part= " << part << endl;
	std::cout << " part " << ipart << " pdg= " << part->GetPDG()->PdgCode() << std::endl;
	part->Print();
      }
      int pdg = part->GetPDG()->PdgCode();
      // This would insure ordering
      if ( pdg == -11 ) { ((TParticle*)part_array->At(ipart))->Momentum(p4ep); cep++; }
      if ( pdg == 11  ) { ((TParticle*)part_array->At(ipart))->Momentum(p4em); cem++; }
      if ( pdg == 22 && !first ) { ((TParticle*)part_array->At(ipart))->Momentum(p4g2); cg++; }
      if ( pdg == 22 && first ) { ((TParticle*)part_array->At(ipart))->Momentum(p4g1); cg++; first = false; }
    }
    p4pi0 = p4g1 + p4g2;
    p4jpsi = p4ep + p4em;

    // check event type
    if ( cep!=1 || cem!=1 || cg!=2 ) {
      std::cout << "Event wrong type " << cep << " e+ " << cem << " e- " << cg << " gamma found where (1,1,2) expected" << std::endl;
      continue;
    }

    // check acceptance
    //if (!accept_pc(p4pim)||!accept_pc(p4pip)||accept_pi0(p4pi0)) continue;
    vector<TLorentzVector> p4s;
    p4s.push_back(p4pbar);
    p4s.push_back(p4g1);
    p4s.push_back(p4g2);
    p4s.push_back(p4ep);
    p4s.push_back(p4em);
    p4s.push_back(p4pi0);
    p4s.push_back(p4jpsi);

    for (auto fill: fillers_inv) (*fill)(p4s);
    for (auto fill: lab_fillers) (*fill)(p4s);
    for (auto fill: cm_fillers) (*fill)(p4s,boost_to_cm);
    for (auto fill: jpsif_fillers) (*fill)(p4s,-p4jpsi.BoostVector());
    for (auto fill: pi0f_fillers) (*fill)(p4s,-p4pi0.BoostVector());
    for (auto fill: var2d_fillers) (*fill)(p4s);
  }

  TFile *fout = TFile::Open("out_sig.root","RECREATE");
  fout->cd();

  for (auto fill: fillers_inv) fill->Write();
  for (auto fill: lab_fillers) fill->Write();
  for (auto fill: cm_fillers) fill->Write();
  for (auto fill: jpsif_fillers) fill->Write();
  for (auto fill: pi0f_fillers) fill->Write();
  for (auto fill: var2d_fillers) fill->Write();

  fout->Write();
  fout->Close();

}

bool accept_pc(const TLorentzVector &v) {
  return true;
}

bool accept_pi0(const TLorentzVector &v) {
  return true;
}

bool accept_elec(const TLorentzVector &v) {
  return true;
}

void print_p4(const char *n, const TLorentzVector &v) {
  std::cout << n << "= (" << v.Px() << ", " << v.Py() << ", " << v.Pz() << ", " << v.E() << ")";
  std::cout << " B= ( " << v.BoostVector().x() << ", " << v.BoostVector().y() << ", " << v.BoostVector().z() << " )" ;
  std::cout <<  std::endl;
}

void set_init_cond(TVector3 &boost_to_cm, TLorentzVector &pbar, TLorentzVector &p, TLorentzVector &pbarp,
		   TLorentzVector &pbar_cm, TLorentzVector &p_cm, TLorentzVector &pbarp_cm) {

  const double p_antip = 5.513;
  const double E_antip = TMath::Hypot(mass_prot, p_antip);
  const double beta_cm = p_antip/(E_antip + mass_prot);
  cout << "betac_cm = " << beta_cm << endl;
  boost_to_cm.SetZ(-beta_cm);
  boost_to_cm.Print();
  pbar.SetPxPyPzE(0,0,p_antip,E_antip);
  p.SetPxPyPzE(0,0,0,mass_prot);
  pbarp = pbar + p;
  cout << "======= LAB ======= " << endl;
  print_p4(" p4pbar_lab= ", pbar);
  print_p4("    p4p_lab= ", p);
  print_p4("p4pbarp_lab= ", pbarp);
  cout << "======================" << endl;
  pbar_cm = pbar;
  p_cm = p;
  pbarp_cm = pbarp;
  pbar_cm.Boost(boost_to_cm);
  p_cm.Boost(boost_to_cm);
  pbarp_cm.Boost(boost_to_cm);
  cout << "======= CM ======= " << endl;
  print_p4(" p4pbar_cm= ", pbar_cm);
  print_p4("    p4p_cm= ", p_cm);
  print_p4("p4pbarp_cm= ", pbarp_cm);
  cout << "======================" << endl;
  // minima and maxima for histogram filling... [TODO] place them somewhere better
  //const double mass_diff = pow( (pow(mass_jpsi,2)-pow(mass_pi0,2))/p4pbarp.M()/2.0, 2);
  //const double t_0 = mass_diff - (p4pbar_cm - p4p_cm).M2();
  //const double t_1 = mass_diff - (p4pbar_cm + p4p_cm).M2();
  //std::cout << " pdg p307: mass_diff= " <<  mass_diff << " t_0 = " << t_0 << " t_1 = " << t_1 << std::endl;
}

void print_event(const vector<int> &evt) {
    cout << "---> Event pdg codes : (";
    for (unsigned int ipdg=0; ipdg<evt.size(); ++ipdg) {
      cout << evt[ipdg];
      if (ipdg!=(evt.size()-1)) cout << ",";
    }
    cout << ")" << endl;
}

void print_log(int ievt, int npart, vector<int> &pdg_codes){
  std::cout << "Event(Type3Pi) " << ievt << " Tca.entries= " << npart << std::endl;
  print_event(pdg_codes);
  cout << endl;
}

bool compare_sorted(const vector<int> &evt, const vector<int> &ref) {
  if (evt.size()!=ref.size()) return false;
  for (unsigned int ievt=0; ievt<evt.size(); ++ievt) {
    if (evt[ievt]!=ref[ievt]) return false;
  }
  return true;
}

bool compare_unsorted(const vector<int> &evt, const vector<int> &ref) {

  if (evt.size()!=ref.size()) return false;
  vector<int> visited;
  for (unsigned int ievt=0; ievt<evt.size(); ++ievt) {

    // search in reference, and mark index visited if foun
    bool found = false;
    for (unsigned int iref=0; iref<ref.size(); ++iref) {
      if ( evt[ievt] == ref[iref]) {
	found = true;

	// check if this index of ref was already
	bool was_visited = false;
	for (unsigned int ivis=0; ivis<visited.size(); ++ivis) {
	  if (int(iref) == visited[ivis]) {
	    was_visited = true;
	    break;
	  }
	}
	if (was_visited) continue;

	visited.push_back(iref);
	break;
      }
    }

    if (!found) return false;

  }
  return true;
}

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
#include "TRandom.h"
#include "TF1.h"
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

static int nbins = 200;
static axis mom_bins(nbins, 0.0, 6.1);
static axis ene_bins(nbins, 0.0, 6.0);
static axis t_bins(nbins, -2.0, 2.0);
static axis the_bins(nbins, -1, 181);
static axis phi_bins(nbins, -181, 181);
static axis cost_bins(nbins, -1.1, 1.1);
static axis oa_bins(nbins, -1, 181);

//static const double mass_jpsi = 3.096;
static const double mass_pi0 = 0.1349766;
static const double mass_prot= 0.938;

static const string eff_tag[2] = {"_noeff", "_eff"};
static const string filt_tag[2] = {"_nofilt", ""};
static int type = 0;
static int ieff = 0;
static int ifilt = 1;

// Efficiency curve for e+ and e-
TH2F* eff_ep_em;
// Efficiency curve for pi0 vs. OA of g-g
TF1* eff_pi0 = new TF1("eff_pi0","0.812 * (1-exp (-0.398 * x) )", 0., 180.);

bool compare_unsorted(const vector<int>&, const vector<int>&);
bool compare_sorted(const vector<int>&, const vector<int>&);

void print_event(const vector<int>&);
void print_log(int, int, vector<int>&);

void decay_pi0(const TLorentzVector&, TLorentzVector&, TLorentzVector&);
bool in_jpsi_window(const TLorentzVector&, const TLorentzVector&);
bool in_t_window(const TLorentzVector&, const TLorentzVector&);
bool accept_pc(const TLorentzVector &v);
bool accept_pi0(const TLorentzVector&,const TLorentzVector&, const TLorentzVector&);
bool accept_elec(const TLorentzVector &v);
double weight_elec(const TLorentzVector &v);
double weight_pc(const TLorentzVector &v);
double weight_pi0(const TLorentzVector&,const TLorentzVector&, const TLorentzVector&);

void fill_bg_hists(TChain *, TClonesArray*);
void fill_sig_hists(TChain *, TClonesArray*);

void print_p4(const char*, const TLorentzVector&);
void set_init_cond(TVector3&, TLorentzVector&, TLorentzVector&, TLorentzVector&,TLorentzVector&, TLorentzVector&, TLorentzVector&);

int main(const int argc, const char **argv) {

  gRandom->SetSeed();

  if (argc<4) {
    std::cout << "Need three arguments [type(0->Bg,1->Sig), eff(0->yes,1->no), filt(0->yes, 1->no)]" << std::endl;
    return -1;
  }

  type = atoi(argv[1]);
  if (type>1) { cout << "type " << type << " out of range, filling BG (default)" << endl; type = 1;}
  ieff = atoi(argv[2]);
  if (ieff>1) { cout << "type " << ieff << " out of range, applying efficiency (default)" << endl; ieff = 1;}

  string list_file;
  ifilt = atoi(argv[3]);
  if (ifilt>1) { cout << "ifilt "<< ifilt << "out of range, using filtered input by default" << endl; ifilt = 1; }

  if (type==1) {
    list_file = (ifilt==0)?"lists/sig_unfilt.list":"lists/sig.list";
    cout << "Options: SIG: ieff= " << ieff << " ifilt = " << ifilt << endl;
  } else {
    list_file = "lists/bg.list";
  }

  // Grab the efficiency histograms and functions
  TFile *f_eff_ep_em = TFile::Open("eff/epem_smooth.root");
  eff_ep_em = (TH2F*) f_eff_ep_em->Get("eff_ep_em_rad")->Clone("eff_ep_em");
  //f_eff_ep_em->Close();

  TChain *data_in = new TChain("data","data_in");
  ifstream inf;
  inf.open(list_file.c_str());
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
  const axis mass_sq_bins(nbins,0,m_max*m_max);

  enum {pbar=0, pi0=1, pim=2, pip=3, pipm=4, g1=5, g2=6};
  const char *names[][2] = {
    {"pbar","#bar{p}"}, {"pi0","#pi^{0}"}, {"pim","#pi^{-}"}, {"pip","#pi^{+}"},
    {"pipm","(#pi^{+}#pi^{-})"}, {"g1","#gamma_{1}"}, {"g2","#gamma_{2}"}
  };
  std::vector<filler*> lab_fillers;

  const char *lab_frame[] = {"lab",""};
  // (pi0), (pi+) and (pi-) mass, mom and angles
  for (int ipart=pi0; ipart<=pipm; ++ipart) {
    lab_fillers.push_back(new var1d(ipart, "mom", mom_bins, lab_frame, names));
    lab_fillers.push_back(new var1d(ipart, "the", the_bins, lab_frame, names));
    lab_fillers.push_back(new var1d(ipart, "the", ipart==pi0?axis(200,0,60):the_bins, lab_frame, names));
    lab_fillers.push_back(new var1d(ipart, "cost", cost_bins, lab_frame, names));
    lab_fillers.push_back(new var1d(ipart, "phi", phi_bins, lab_frame, names));
  }
  // (pi+) - (pi-) system mass, mom and angles
  lab_fillers.push_back(new var1d(pip, pim, "mom", mom_bins, lab_frame, names));
  lab_fillers.push_back(new var1d(pip, pim, "the", the_bins, lab_frame, names));
  lab_fillers.push_back(new var1d(pip, pim, "phi", phi_bins, lab_frame, names));
  lab_fillers.push_back(new var1d(pip, pim, "oa", oa_bins, lab_frame, names));
  lab_fillers.push_back(new var1d(pi0, pipm, "oa", oa_bins, lab_frame, names));

  lab_fillers.push_back(new var1d(pip, pim, "mom_fwd", mom_bins, lab_frame, names)); // (e+-e-) OA
  lab_fillers.push_back(new var1d(pip, pim, "mom_bwd", mom_bins, lab_frame, names)); // (e+-e-) OA
  lab_fillers.push_back(new var1d(pip, pim, "the_fwd", the_bins, lab_frame, names)); // (e+-e-) OA
  lab_fillers.push_back(new var1d(pip, pim, "the_bwd", the_bins, lab_frame, names)); // (e+-e-) OA

  lab_fillers.push_back(new var2d(pip, pim, "mom_fwd", mom_bins, lab_frame, pip, pim, "the_fwd", the_bins, lab_frame, names));
  lab_fillers.push_back(new var2d(pip, pim, "mom_bwd", mom_bins, lab_frame, pip, pim, "the_bwd", the_bins, lab_frame, names));

  lab_fillers.push_back(new var1d(g1, g2, "oa", axis(200,0,180), lab_frame, names)); // (g1-g2) OA
  lab_fillers.push_back(new var2d(pi0, "e", axis(200,0,6), lab_frame, g1, g2, "oa", axis(200,0,60), lab_frame, names));

  lab_fillers.push_back(new var2d(pip, "mom", mom_bins, lab_frame, pip, "the", the_bins, lab_frame, names));
  lab_fillers.push_back(new var2d(pim, "mom", mom_bins, lab_frame, pim, "the", the_bins, lab_frame, names));

  std::vector<filler*> inv_fillers;
  const char *invariant[] = {"inv", ""};
  inv_fillers.push_back(new var2d(pi0, pip, "mass", mass_bins, invariant, pip, pim, "mass", mass_bins, invariant, names) );
  inv_fillers.push_back(new var2d(pi0, pip, "mass_sq", mass_sq_bins, invariant, pip, pim, "mass_sq", mass_sq_bins, invariant, names) );

  inv_fillers.push_back(new var1d(pip, pim, "mass", mass_bins, invariant, names));
  inv_fillers.push_back(new var1d(pi0, pip, "mass", mass_bins, invariant, names));
  inv_fillers.push_back(new var1d(pi0, pim, "mass", mass_bins, invariant, names));
  // Mandelstam variables
  inv_fillers.push_back(new var1d(pbar, pipm, "u", axis(200,-11,1), invariant, names));
  inv_fillers.push_back(new var1d(pbar, pi0, "t", axis(200,-11,1), invariant, names));
  // "Dalitz" plot
  inv_fillers.push_back(new var2d(pbar, pi0, "t", axis(200,-11,1), invariant, pi0, "the", the_bins, lab_frame, names));

  std::vector<filler*> cm_fillers;
  const char *cm_frame[] = {"cm"," (CM frame)"};
  cm_fillers.push_back(new var1d(pi0, "mom", mom_bins, cm_frame, names));
  cm_fillers.push_back(new var1d(pi0, "the", the_bins, cm_frame, names));
  cm_fillers.push_back(new var1d(pi0, "cost", cost_bins, cm_frame, names));
  cm_fillers.push_back(new var1d(pip, "the", the_bins, cm_frame, names));
  cm_fillers.push_back(new var1d(pim, "the", the_bins, cm_frame, names));
  cm_fillers.push_back(new var1d(pip, pim, "mom", mom_bins, cm_frame, names));
  cm_fillers.push_back(new var1d(pip, pim, "the", the_bins, cm_frame, names));
  cm_fillers.push_back(new var1d(pip, pim, "oa", axis(200,0,60), cm_frame, names));
  cm_fillers.push_back(new var1d(pi0, pipm, "oa", axis(200,0,60), cm_frame, names));
  cm_fillers.push_back(new var1d(pi0, pipm, "oa", axis(200,0,60), cm_frame, names));
  cm_fillers.push_back(new var1d(g1, g2, "oa", oa_bins, cm_frame, names)); // (g1-g2) OA

  std::vector<filler*> pi0f_fillers;
  const char *pi0_frame[] {"pi0f"," (#pi^{0} frame)"};
  pi0f_fillers.push_back(new var1d(g1, "mom", mom_bins, pi0_frame, names)); // (g1) mom
  pi0f_fillers.push_back(new var1d(g1, "cost", cost_bins, pi0_frame, names)); // (g1) cos(theta)
  pi0f_fillers.push_back(new var1d(g2, "mom", mom_bins, pi0_frame, names)); // (g2) mom
  pi0f_fillers.push_back(new var1d(g2, "cost", cost_bins, pi0_frame, names)); // (g2) cos(theta)
  pi0f_fillers.push_back(new var1d(g1, g2, "oa", oa_bins, pi0_frame, names)); // (g1-g2) OA

  std::vector<filler*> cmVlab_fillers;
  // think how to pass "no boost" lab frame in the same filler with another frame
  cmVlab_fillers.push_back(new var2d(g1, g2, "oa", oa_bins, pi0_frame, g1, g2, "oa", oa_bins, lab_frame, names)); // (g1-g2) OA

  const int Nevt = data_in->GetEntries();

  for (int ievt=0; ievt<Nevt; ievt++) {

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

    // phasespace decay pi0
    TLorentzVector p4g1, p4g2;
    decay_pi0(p4pi0,p4g1,p4g2);

    // check event type
    if ( cpim!=1 || cpip!=1 || cpi0!=1 ) {
      std::cout << "Event wrong type " << cpim << " pi- " << cpip << " pi+ " << cpi0 << " pi0 found where 1 expected" << std::endl;
      continue;
    }

    if (ifilt==1) {
      if (!in_jpsi_window(p4pip,p4pim) ) continue;
      if (!in_t_window(p4pbar, p4pi0)) continue;
    }

    // apply efficiency
    //if (ieff==1) {
    //  if (!accept_pc(p4pim)||!accept_pc(p4pip)||!accept_pi0(p4pi0,p4g1,p4g2)) continue;
    //}
    double weight = 1.0;
    if (ieff==1) weight = weight_pc(p4pim)*weight_pc(p4pip)*weight_pi0(p4pi0,p4g1,p4g2);

    vector<TLorentzVector> p4s;
    p4s.push_back(p4pbar);
    p4s.push_back(p4pi0);
    p4s.push_back(p4pim);
    p4s.push_back(p4pip);
    p4s.push_back(p4pipm);
    p4s.push_back(p4g1);
    p4s.push_back(p4g2);

    for (auto fill: lab_fillers) (*fill)(p4s,weight);
    for (auto fill: inv_fillers) (*fill)(p4s,weight);
    for (auto fill: cm_fillers) (*fill)(p4s,boost_to_cm,weight);
    for (auto fill: pi0f_fillers) (*fill)(p4s,-p4pi0.BoostVector(),weight);

  }


  TFile *fout = TFile::Open(Form("hists/bg%s%s.root",filt_tag[ifilt].c_str(),eff_tag[ieff].c_str()),"RECREATE");
  //TFile *fout = TFile::Open(Form("hists/bg%s.root",eff_tag[ieff].c_str()),"RECREATE");
  fout->cd();
  for (auto fill: lab_fillers) fill->Write();
  for (auto fill: cm_fillers) fill->Write();
  for (auto fill: inv_fillers) fill->Write();
  for (auto fill: pi0f_fillers) fill->Write();
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

  enum {pbar=0, g1=1, g2=2, ep=3, em=4, pi0=5, jpsi=6, g3=7, g4=8};
  const char *names[][2] = {
    {"pbar","#bar{p}"}, {"g1","#gamma_{1}"}, {"g2","#gamma_{2}"},
    {"ep","e^{+}"}, {"em","e^{-}"}, {"pi0","#pi^{0}"}, {"jpsi", "J/#psi"},
    {"g3","#gamma_{3}"}, {"g4","#gamma_{4}"},
  };

  std::vector<filler*> inv_fillers;
  const char *invariant[] = {"inv", ""};
  inv_fillers.push_back( new var1d(pi0, jpsi, "mass", mass_bins, invariant, names));
  inv_fillers.push_back( new var1d(ep, em, "mass", mass_bins, invariant, names));
  inv_fillers.push_back( new var1d(g1, g2, "mass", mass_bins, invariant, names));
  inv_fillers.push_back( new var1d(pbar, jpsi, "u", t_bins, invariant, names));
  inv_fillers.push_back( new var1d(pbar, pi0, "t", t_bins, invariant, names));

  std::vector<filler*> lab_fillers;
  const char *lab_frame[] = {"lab",""};
  for (int ipart=g1; ipart<= jpsi; ++ipart) {   // for all single particle
    lab_fillers.push_back(new var1d(ipart, "mom", mom_bins, lab_frame, names));
    lab_fillers.push_back(new var1d(ipart, "e", ene_bins, lab_frame, names));
    lab_fillers.push_back(new var1d(ipart, "the", the_bins, lab_frame, names));
    lab_fillers.push_back(new var1d(ipart, "cost", cost_bins, lab_frame, names));
  }

  lab_fillers.push_back(new var1d(ep, em, "oa", oa_bins, lab_frame, names)); // (e+-e-) OA
  lab_fillers.push_back(new var1d(g1, g2, "oa", axis(200,0,180), lab_frame, names)); // (g1-g2) OA
  lab_fillers.push_back(new var1d(g3, g4, "oa", axis(200,0,180), lab_frame, names)); // (g1-g2) OA
  lab_fillers.push_back(new var1d(pi0, jpsi, "oa", oa_bins, lab_frame, names)); // (pi0-jpsi) OA
  lab_fillers.push_back(new var1d(ep, em, "mom_fwd", mom_bins, lab_frame, names)); // (e+-e-) OA
  lab_fillers.push_back(new var1d(ep, em, "mom_bwd", mom_bins, lab_frame, names)); // (e+-e-) OA
  lab_fillers.push_back(new var1d(ep, em, "the_fwd", the_bins, lab_frame, names)); // (e+-e-) OA
  lab_fillers.push_back(new var1d(ep, em, "the_bwd", the_bins, lab_frame, names)); // (e+-e-) OA

  lab_fillers.push_back(new var2d(pi0, "e", axis(200,0,6), lab_frame, g1, g2, "oa", axis(200,0,60), lab_frame, names));
  lab_fillers.push_back(new var2d(pi0, "the", the_bins, lab_frame, g1, g2, "oa", axis(200,0,60), lab_frame, names));
  lab_fillers.push_back(new var2d(ep, "mom", mom_bins, lab_frame, ep, "the", the_bins, lab_frame, names));
  lab_fillers.push_back(new var2d(em, "mom", mom_bins, lab_frame, em, "the", the_bins, lab_frame, names));

  lab_fillers.push_back(new var2d(pi0, "the", the_bins, lab_frame, pi0, "mom", mom_bins, lab_frame, names));

  //lab_fillers.push_back(new var1d(ep, em, "the", the_bins, lab_frame, names)); // -- stupid sanity check
  lab_fillers.push_back(new var2d(pbar, pi0, "t", t_bins, invariant, pi0, "the", the_bins, lab_frame, names));
  lab_fillers.push_back(new var2d(ep, "the", the_bins, lab_frame, ep, "mom", mom_bins, lab_frame, names));
  lab_fillers.push_back(new var2d(em, "the", the_bins, lab_frame, em, "mom", mom_bins, lab_frame, names));

  lab_fillers.push_back(new var2d(ep, em, "mom_fwd", mom_bins, lab_frame, ep, em, "the_fwd", the_bins, lab_frame, names));
  lab_fillers.push_back(new var2d(ep, em, "mom_bwd", mom_bins, lab_frame, ep, em, "the_bwd", the_bins, lab_frame, names));

  std::vector<filler*> cm_fillers;
  const char *cm_frame[] = {"cm"," (CM frame)"};
  for (int ipart=g1; ipart<=jpsi; ++ipart) {
    cm_fillers.push_back(new var1d(ipart, "mom", mom_bins, cm_frame, names));
    cm_fillers.push_back(new var1d(ipart, "e", ene_bins, cm_frame, names));
    cm_fillers.push_back(new var1d(ipart, "the", the_bins, cm_frame, names));
    cm_fillers.push_back(new var1d(ipart, "cost", cost_bins, cm_frame, names));
  }
  cm_fillers.push_back(new var1d(ep, em, "oa", oa_bins, cm_frame, names)); // (e+-e-) OA
  cm_fillers.push_back(new var1d(g1, g2, "oa", axis(200,0,180), cm_frame, names)); // (g1-g2) OA
  cm_fillers.push_back(new var1d(pi0, jpsi, "oa", oa_bins, cm_frame, names)); // (pi0-jpsi) OA

  cm_fillers.push_back(new var2d(pbar, pi0, "t", t_bins, invariant, pi0, "the", the_bins, cm_frame, names));
  cm_fillers.push_back(new var2d(pbar, jpsi, "u", t_bins, invariant, pi0, "the", the_bins, cm_frame, names));
  cm_fillers.push_back(new var2d(pbar, pi0, "t", t_bins, invariant, pi0, "cost", cost_bins, cm_frame, names));
  cm_fillers.push_back(new var2d(pbar, jpsi, "u", t_bins, invariant, pi0, "cost", cost_bins, cm_frame, names));

  std::vector<filler*> jpsif_fillers;
  const char *jpsi_frame[] {"jpsif"," (J/#psi frame)"};
  jpsif_fillers.push_back(new var1d(ep, "mom", mom_bins, jpsi_frame, names)); // (e+) mom
  jpsif_fillers.push_back(new var1d(ep, "cost", cost_bins, jpsi_frame, names)); // (e+) cos(theta)
  jpsif_fillers.push_back(new var1d(em, "mom", mom_bins, jpsi_frame, names)); // (e+) mom
  jpsif_fillers.push_back(new var1d(em, "cost", cost_bins, jpsi_frame, names)); // (e+) cos(theta)
  jpsif_fillers.push_back(new var1d(ep, em, "oa", cost_bins, jpsi_frame, names)); // (e+-e-) OA

  std::vector<filler*> pi0f_fillers;
  const char *pi0_frame[] {"pi0f"," (#pi^{0} frame)"};
  pi0f_fillers.push_back(new var1d(g1, "mom", mom_bins, pi0_frame, names)); // (g1) mom
  pi0f_fillers.push_back(new var1d(g1, "cost", cost_bins, pi0_frame, names)); // (g1) cos(theta)
  pi0f_fillers.push_back(new var1d(g2, "mom", mom_bins, pi0_frame, names)); // (g2) mom
  pi0f_fillers.push_back(new var1d(g2, "cost", cost_bins, pi0_frame, names)); // (g2) cos(theta)
  pi0f_fillers.push_back(new var1d(g1, g2, "oa", oa_bins, pi0_frame, names)); // (g1-g2) OA

  //pi0f_fillers.push_back(new var1d(g3, "mom", mom_bins, pi0_frame, names)); // (g3) mom
  //pi0f_fillers.push_back(new var1d(g3, "cost", cost_bins, pi0_frame, names)); // (g3) cos(theta)
  //pi0f_fillers.push_back(new var1d(g4, "mom", mom_bins, pi0_frame, names)); // (g4) mom
  //pi0f_fillers.push_back(new var1d(g4, "cost", cost_bins, pi0_frame, names)); // (g4) cos(theta)
  //pi0f_fillers.push_back(new var1d(g3, g4, "oa", oa_bins, pi0_frame, names)); // (g3-g4) OA

  std::vector<filler*> var2d_fillers;
  var2d_fillers.push_back(new var2d(g1, "mom", mom_bins, lab_frame, g2, "mom", mom_bins, lab_frame, names) );

  const int Nevt = data_in->GetEntries();

  for (int ievt=0; ievt<Nevt; ievt++) {

    //if (ievt>100000) break;
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

    //cout << "pi0mass = " << p4pi0.M() << endl;;

    // phasespace decay pi0
    TLorentzVector p4g3, p4g4;
    decay_pi0(p4pi0,p4g3,p4g4);

    // check event type
    if ( cep!=1 || cem!=1 || cg!=2 ) {
      std::cout << "Event wrong type " << cep << " e+ " << cem << " e- " << cg << " gamma found where (1,1,2) expected" << std::endl;
      continue;
    }

    // check acceptance
    //if (ieff==1) {
    //  if (!accept_elec(p4em)||!accept_elec(p4ep)||!accept_pi0(p4pi0,p4g1,p4g2)) continue;
    //}
    double weight = 1.0;
    if (ieff==1) weight = weight_elec(p4em)*weight_elec(p4ep)*weight_pi0(p4pi0,p4g1,p4g2);

    vector<TLorentzVector> p4s;
    p4s.push_back(p4pbar);
    p4s.push_back(p4g1);
    p4s.push_back(p4g2);
    p4s.push_back(p4ep);
    p4s.push_back(p4em);
    p4s.push_back(p4pi0);
    p4s.push_back(p4jpsi);
    p4s.push_back(p4g3);
    p4s.push_back(p4g4);

    for (auto fill: inv_fillers) (*fill)(p4s,weight);
    for (auto fill: lab_fillers) (*fill)(p4s,weight);
    for (auto fill: cm_fillers) (*fill)(p4s,boost_to_cm,weight);
    for (auto fill: jpsif_fillers) (*fill)(p4s,-p4jpsi.BoostVector(),weight);
    for (auto fill: pi0f_fillers) (*fill)(p4s,-p4pi0.BoostVector(),weight);
    for (auto fill: var2d_fillers) (*fill)(p4s,weight);

  }

  TFile *fout = TFile::Open(Form("hists/sig%s%s.root",filt_tag[ifilt].c_str(),eff_tag[ieff].c_str()),"RECREATE");
  fout->cd();

  for (auto fill: inv_fillers) fill->Write();
  for (auto fill: lab_fillers) fill->Write();
  for (auto fill: cm_fillers) fill->Write();
  for (auto fill: jpsif_fillers) fill->Write();
  for (auto fill: pi0f_fillers) fill->Write();
  for (auto fill: var2d_fillers) fill->Write();

  fout->Write();
  fout->Close();

}

void decay_pi0(const TLorentzVector &pi0, TLorentzVector &g1, TLorentzVector &g2) {

  const TVector3 boost_to_lab = pi0.BoostVector();
  const double E = mass_pi0/2., p=E;
  const double cost = (2.0*gRandom->Uniform())-1.0;
  const double pz = p * cost;
  const double pt = sqrt(p*p-pz*pz);
  const double phi = gRandom->Uniform();
  const double py = pt * sin(phi);
  const double px = pt * cos(phi);

  // GO FIGURE -
  // TODO still not quite there but close enough for now
  TVector3 g(px,py,pz);
  TVector3 dir = -pi0.Vect().Unit();
  g.RotateUz(dir);
  g1.SetPxPyPzE(g.X(),g.Y(),g.Z(),E);
  g1.Boost(boost_to_lab);
  g2.SetPxPyPzE(-g.X(),-g.Y(),-g.Z(),E);
  g2.Boost(boost_to_lab);

}

bool in_jpsi_window(const TLorentzVector &vp, const TLorentzVector &vm) {
  const double m = (vp+vm).M();
  return ( 2.96 < m && m < 3.22 );
}

bool in_t_window(const TLorentzVector &vpbar, const TLorentzVector &vpi0) {
  const double t = (vpbar-vpi0).M2();
  return ( -0.5 < t && t < 0.6 );
}

bool accept_pc(const TLorentzVector &v) {
  return true;
}

double weight_pc(const TLorentzVector &v) {
  return 3e-4;
}

bool accept_pi0(const TLorentzVector &v, const TLorentzVector &vg1, const TLorentzVector &vg2) {
  const double prob = gRandom->Uniform();
  const double eff = eff_pi0->Eval(rtd*vg1.Vect().Angle(vg2.Vect()));
  return (prob<eff);
  return true;
}

double weight_pi0(const TLorentzVector &v, const TLorentzVector &vg1, const TLorentzVector &vg2) {
  return eff_pi0->Eval(rtd*vg1.Vect().Angle(vg2.Vect()));
}


bool accept_elec(const TLorentzVector &v) {
  const int ix = eff_ep_em->GetXaxis()->FindBin(v.Vect().Mag());
  const int iy = eff_ep_em->GetYaxis()->FindBin(v.Vect().Theta());
  const double eff = eff_ep_em->GetBinContent(ix,iy);
  const double prob = gRandom->Uniform();
  //cout << "prob= " << prob << " eff= " << eff << " acc?= " << (prob<eff) << endl;
  return (prob<eff);
}

double weight_elec(const TLorentzVector &v) {
  const int ix = eff_ep_em->GetXaxis()->FindBin(v.Vect().Mag());
  const int iy = eff_ep_em->GetYaxis()->FindBin(v.Vect().Theta());
  return eff_ep_em->GetBinContent(ix,iy);
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

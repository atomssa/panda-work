//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 10 20:37:37 2015 by ROOT version 5.34/05
// from TTree t/t
// found on file: ../../grid.out/esim_oct14_binsong_config8/runall.1.out/bremcorr.root
//////////////////////////////////////////////////////////

#ifndef t_h
#define t_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

static Double_t config[22][3] = {
  {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
  {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
  {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
  {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
  {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
  {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

static const int nch_max = 10;
static const int nmcb_max = 200;
static const int nsb_max = 100;
static const int npb_max = 100;
static const int nab_max = 100;

class t {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Int_t           nch;
  Int_t           charge[nch_max];   //[nch]
  Float_t         mom_mc[nch_max];   //[nch]
  Float_t         mom_rec[nch_max];   //[nch]
  Float_t         mom_cor[nch_max];   //[nch]
  Float_t         mom_wcor[nch_max];   //[nch]
  Float_t         mom_mrg[nch_max];   //[nch]
  Float_t         mom_sep[nch_max];   //[nch]
  Float_t         mom_stored[nch_max];   //[nch]
  Float_t         mom_sep_w[nch_max];   //[nch]
  Float_t         mom_sep_w_bf[nch_max];   //[nch]
  Float_t         mom_mrg_w[nch_max];   //[nch]
  Float_t         mom_mrg_w_bf[nch_max];   //[nch]
  Float_t         mom_mrg_pc[nch_max];   //[nch]
  Float_t         mom_mrg_w_pc[nch_max];   //[nch]
  Float_t         mom_mrg_w_bf_pc[nch_max];   //[nch]
  Float_t         mom_wbfcor[nch_max];   //[nch]
  Float_t         mom_wbfcor_mw_bf[nch_max];   //[nch]
  Float_t         mom_wbfcor_mw_bf_pc[nch_max];   //[nch]
  Float_t         phi[nch_max];   //[nch]
  Float_t         the[nch_max];   //[nch]
  Float_t         phi_mc[nch_max];   //[nch]
  Float_t         the_mc[nch_max];   //[nch]
  Int_t           pdg_mc[nch_max];   //[nch]
  Int_t           nphot_sep[nch_max];   //[nch]
  Int_t           nphot_mrg[nch_max];   //[nch]
  Int_t           nphot_mrg_pc[nch_max];   //[nch]
  Int_t           is_prim[nch_max];   //[nch]
  Int_t           _nmcb[nch_max];   //[nch]
  Int_t           imcb_s[nch_max];   //[nch]
  Int_t           imcb_e[nch_max];   //[nch]
  Int_t           nmcb;
  Float_t         mcb_phi[57];   //[nmcb]
  Float_t         mcb_the[57];   //[nmcb]
  Float_t         mcb_ene[57];   //[nmcb]
  Float_t         mcb_rad[57];   //[nmcb]
  Float_t         mcb_zed[57];   //[nmcb]
  Int_t           mcb_match[57];   //[nmcb]
  Int_t           mcb_score[57];   //[nmcb]
  Int_t           mcb_match_ab[57];   //[nmcb]
  Int_t           mcb_score_ab[57];   //[nmcb]
  Int_t           _nsb[nch_max];   //[nch]
  Int_t           isb_s[nch_max];   //[nch]
  Int_t           isb_e[nch_max];   //[nch]
  Int_t           nsb;
  Float_t         sb_idx[nsb_max];   //[nsb]
  Float_t         sb_phi[nsb_max];   //[nsb]
  Float_t         sb_the[nsb_max];   //[nsb]
  Float_t         sb_ene[nsb_max];   //[nsb]
  Float_t         sb_rcalc[nsb_max];   //[nsb]
  Float_t         sb_zcalc[nsb_max];   //[nsb]
  Int_t           sb_match[nsb_max];   //[nsb]
  Int_t           sb_score[nsb_max];   //[nsb]
  Int_t           _npb[nch_max];   //[nch]
  Int_t           ipb_s[nch_max];   //[nch]
  Int_t           ipb_e[nch_max];   //[nch]
  Int_t           npb;
  Int_t           pb_acc[npb_max];   //[npb]
  Float_t         pb_phi[npb_max];   //[npb]
  Float_t         pb_the[npb_max];   //[npb]
  Float_t         pb_ene[npb_max];   //[npb]
  Float_t         pb_rcalc[npb_max];   //[npb]
  Float_t         pb_zcalc[npb_max];   //[npb]
  Int_t           nab;
  Float_t         ab_phi[nab_max];   //[nab]
  Float_t         ab_the[nab_max];   //[nab]
  Float_t         ab_ene[nab_max];   //[nab]
  Int_t           ab_isb[nab_max];   //[nab]
  Int_t           ab_ich[nab_max];   //[nab]
  Int_t           ab_match[nsb_max];   //[nsb]
  Int_t           ab_score[nsb_max];   //[nsb]

  // List of branches
  TBranch        *b_nch;   //!
  TBranch        *b_charge;   //!
  TBranch        *b_mom_mc;   //!
  TBranch        *b_mom_rec;   //!
  TBranch        *b_mom_cor;   //!
  TBranch        *b_mom_wcor;   //!
  TBranch        *b_mom_mrg;   //!
  TBranch        *b_mom_sep;   //!
  TBranch        *b_mom_stored;   //!
  TBranch        *b_mom_sep_w;   //!
  TBranch        *b_mom_sep_w_bf;   //!
  TBranch        *b_mom_mrg_w;   //!
  TBranch        *b_mom_mrg_w_bf;   //!
  TBranch        *b_mom_mrg_pc;   //!
  TBranch        *b_mom_mrg_w_pc;   //!
  TBranch        *b_mom_mrg_w_bf_pc;   //!
  TBranch        *b_mom_wbfcor;   //!
  TBranch        *b_mom_wbfcor_mw_bf;   //!
  TBranch        *b_mom_wbfcor_mw_bf_pc;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_the;   //!
  TBranch        *b_phi_mc;   //!
  TBranch        *b_the_mc;   //!
  TBranch        *b_pdg_mc;   //!
  TBranch        *b_nphot_sep;   //!
  TBranch        *b_nphot_mrg;   //!
  TBranch        *b_nphot_mrg_pc;   //!
  TBranch        *b_is_prim;   //!
  TBranch        *b__nmcb;   //!
  TBranch        *b_imcb_s;   //!
  TBranch        *b_imcb_e;   //!
  TBranch        *b_nmcb;   //!
  TBranch        *b_mcb_phi;   //!
  TBranch        *b_mcb_the;   //!
  TBranch        *b_mcb_ene;   //!
  TBranch        *b_mcb_rad;   //!
  TBranch        *b_mcb_zed;   //!
  TBranch        *b_mcb_match;   //!
  TBranch        *b_mcb_score;   //!
  TBranch        *b_mcb_match_ab;   //!
  TBranch        *b_mcb_score_ab;   //!
  TBranch        *b__nsb;   //!
  TBranch        *b_isb_s;   //!
  TBranch        *b_isb_e;   //!
  TBranch        *b_nsb;   //!
  TBranch        *b_sb_idx;   //!
  TBranch        *b_sb_phi;   //!
  TBranch        *b_sb_the;   //!
  TBranch        *b_sb_ene;   //!
  TBranch        *b_sb_rcalc;   //!
  TBranch        *b_sb_zcalc;   //!
  TBranch        *b_sb_match;   //!
  TBranch        *b_sb_score;   //!
  TBranch        *b__npb;   //!
  TBranch        *b_ipb_s;   //!
  TBranch        *b_ipb_e;   //!
  TBranch        *b_npb;   //!
  TBranch        *b_pb_acc;   //!
  TBranch        *b_pb_phi;   //!
  TBranch        *b_pb_the;   //!
  TBranch        *b_pb_ene;   //!
  TBranch        *b_pb_rcalc;   //!
  TBranch        *b_pb_zcalc;   //!
  TBranch        *b_nab;   //!
  TBranch        *b_ab_phi;   //!
  TBranch        *b_ab_the;   //!
  TBranch        *b_ab_ene;   //!
  TBranch        *b_ab_isb;   //!
  TBranch        *b_ab_ich;   //!
  TBranch        *b_ab_match;   //!
  TBranch        *b_ab_score;   //!

  TString fout_name;

  t(TTree *tree=0);
  t(TString fname);
  virtual ~t();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual TH1F*    th1f(TString,TString,int,float,float,int);
  virtual TLegend* tleg(double,double,double,double);
  virtual void     set_style(TH2* h);
  virtual void     set_style(TH1* h, int col);
  virtual Double_t res_func(const double &m, const double &mref) { return (mref-m)/mref; }
  virtual Bool_t   check_invariants(Int_t, Bool_t);
  virtual double   dphi(double, int);
};

#endif

#ifdef t_cxx
t::t(TTree *tree) : fChain(0)
{
  if (tree == 0) {
    gDirectory->GetObject("t",tree);
    if (!tree) {
      // if parameter tree is not specified (or zero), connect the file
      // used to generate this class and read the Tree.
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../grid.out/psim_oct14_binsong_config8/runall.1/bremcorr.root");
      if (!f || !f->IsOpen()) {
	f = new TFile("../../grid.out/psim_oct14_binsong_config8/runall.1/bremcorr.root");
      }
      f->GetObject("t",tree);
    }
  }
  fout_name = "bremcorr_hists.root";
  Init(tree);
  Loop();
}

t::t(TString fname) : fChain(0)
{
  cout << endl;
  cout << endl;
  cout << "t::t(TString) Opening " << fname << endl;
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TTree *tree;
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
  if (!f || !f->IsOpen()) {
    f = new TFile(fname);
  }
  f->GetObject("t",tree);
  fout_name = fname;
  fout_name.ReplaceAll(".root","_hists.root");
  Init(tree);
  Loop();
}

t::~t()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

double t::dphi(double _phi, int ichtrk) {
  return charge[ichtrk]>0?(phi[ichtrk] - _phi):(_phi - phi[ichtrk]);
}

Int_t t::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t t::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void t::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("nch", &nch, &b_nch);
  fChain->SetBranchAddress("charge", charge, &b_charge);
  fChain->SetBranchAddress("mom_mc", mom_mc, &b_mom_mc);
  fChain->SetBranchAddress("mom_rec", mom_rec, &b_mom_rec);
  fChain->SetBranchAddress("mom_cor", mom_cor, &b_mom_cor);
  fChain->SetBranchAddress("mom_wcor", mom_wcor, &b_mom_wcor);
  fChain->SetBranchAddress("mom_mrg", mom_mrg, &b_mom_mrg);
  fChain->SetBranchAddress("mom_sep", mom_sep, &b_mom_sep);
  fChain->SetBranchAddress("mom_stored", mom_stored, &b_mom_stored);
  fChain->SetBranchAddress("mom_sep_w", mom_sep_w, &b_mom_sep_w);
  fChain->SetBranchAddress("mom_sep_w_bf", mom_sep_w_bf, &b_mom_sep_w_bf);
  fChain->SetBranchAddress("mom_mrg_w", mom_mrg_w, &b_mom_mrg_w);
  fChain->SetBranchAddress("mom_mrg_w_bf", mom_mrg_w_bf, &b_mom_mrg_w_bf);
  fChain->SetBranchAddress("mom_mrg_pc", mom_mrg_pc, &b_mom_mrg_pc);
  fChain->SetBranchAddress("mom_mrg_w_pc", mom_mrg_w_pc, &b_mom_mrg_w_pc);
  fChain->SetBranchAddress("mom_mrg_w_bf_pc", mom_mrg_w_bf_pc, &b_mom_mrg_w_bf_pc);
  fChain->SetBranchAddress("mom_wbfcor", mom_wbfcor, &b_mom_wbfcor);
  fChain->SetBranchAddress("mom_wbfcor_mw_bf", mom_wbfcor_mw_bf, &b_mom_wbfcor_mw_bf);
  fChain->SetBranchAddress("mom_wbfcor_mw_bf_pc", mom_wbfcor_mw_bf_pc, &b_mom_wbfcor_mw_bf_pc);
  fChain->SetBranchAddress("phi", phi, &b_phi);
  fChain->SetBranchAddress("the", the, &b_the);
  fChain->SetBranchAddress("phi_mc", phi_mc, &b_phi_mc);
  fChain->SetBranchAddress("the_mc", the_mc, &b_the_mc);
  fChain->SetBranchAddress("pdg_mc", pdg_mc, &b_pdg_mc);
  fChain->SetBranchAddress("nphot_sep", nphot_sep, &b_nphot_sep);
  fChain->SetBranchAddress("nphot_mrg", nphot_mrg, &b_nphot_mrg);
  fChain->SetBranchAddress("nphot_mrg_pc", nphot_mrg_pc, &b_nphot_mrg_pc);
  fChain->SetBranchAddress("is_prim", is_prim, &b_is_prim);
  fChain->SetBranchAddress("_nmcb", _nmcb, &b__nmcb);
  fChain->SetBranchAddress("imcb_s", imcb_s, &b_imcb_s);
  fChain->SetBranchAddress("imcb_e", imcb_e, &b_imcb_e);
  fChain->SetBranchAddress("nmcb", &nmcb, &b_nmcb);
  fChain->SetBranchAddress("mcb_phi", mcb_phi, &b_mcb_phi);
  fChain->SetBranchAddress("mcb_the", mcb_the, &b_mcb_the);
  fChain->SetBranchAddress("mcb_ene", mcb_ene, &b_mcb_ene);
  fChain->SetBranchAddress("mcb_rad", mcb_rad, &b_mcb_rad);
  fChain->SetBranchAddress("mcb_zed", mcb_zed, &b_mcb_zed);
  fChain->SetBranchAddress("mcb_match", mcb_match, &b_mcb_match);
  fChain->SetBranchAddress("mcb_score", mcb_score, &b_mcb_score);
  fChain->SetBranchAddress("mcb_match_ab", mcb_match_ab, &b_mcb_match_ab);
  fChain->SetBranchAddress("mcb_score_ab", mcb_score_ab, &b_mcb_score_ab);
  fChain->SetBranchAddress("_nsb", _nsb, &b__nsb);
  fChain->SetBranchAddress("isb_s", isb_s, &b_isb_s);
  fChain->SetBranchAddress("isb_e", isb_e, &b_isb_e);
  fChain->SetBranchAddress("nsb", &nsb, &b_nsb);
  fChain->SetBranchAddress("sb_idx", sb_idx, &b_sb_idx);
  fChain->SetBranchAddress("sb_phi", sb_phi, &b_sb_phi);
  fChain->SetBranchAddress("sb_the", sb_the, &b_sb_the);
  fChain->SetBranchAddress("sb_ene", sb_ene, &b_sb_ene);
  fChain->SetBranchAddress("sb_rcalc", sb_rcalc, &b_sb_rcalc);
  fChain->SetBranchAddress("sb_zcalc", sb_zcalc, &b_sb_zcalc);
  fChain->SetBranchAddress("sb_match", sb_match, &b_sb_match);
  fChain->SetBranchAddress("sb_score", sb_score, &b_sb_score);
  fChain->SetBranchAddress("_npb", _npb, &b__npb);
  fChain->SetBranchAddress("ipb_s", ipb_s, &b_ipb_s);
  fChain->SetBranchAddress("ipb_e", ipb_e, &b_ipb_e);
  fChain->SetBranchAddress("npb", &npb, &b_npb);
  fChain->SetBranchAddress("pb_acc", pb_acc, &b_pb_acc);
  fChain->SetBranchAddress("pb_phi", pb_phi, &b_pb_phi);
  fChain->SetBranchAddress("pb_the", pb_the, &b_pb_the);
  fChain->SetBranchAddress("pb_ene", pb_ene, &b_pb_ene);
  fChain->SetBranchAddress("pb_rcalc", pb_rcalc, &b_pb_rcalc);
  fChain->SetBranchAddress("pb_zcalc", pb_zcalc, &b_pb_zcalc);
  fChain->SetBranchAddress("nab", &nab, &b_nab);
  fChain->SetBranchAddress("ab_phi", ab_phi, &b_ab_phi);
  fChain->SetBranchAddress("ab_the", ab_the, &b_ab_the);
  fChain->SetBranchAddress("ab_ene", ab_ene, &b_ab_ene);
  fChain->SetBranchAddress("ab_isb", ab_isb, &b_ab_isb);
  fChain->SetBranchAddress("ab_ich", ab_ich, &b_ab_ich);
  fChain->SetBranchAddress("ab_match", ab_match, &b_ab_match);
  fChain->SetBranchAddress("ab_score", ab_score, &b_ab_score);
  Notify();

}

Bool_t t::check_invariants(Int_t ient, Bool_t debug) {
  Bool_t ok = true;
  int _nmcb_sum = 0, _nsb_sum = 0, _npb_sum = 0;
  for (int ich=0; ich<nch; ++ich) {
    if (imcb_e[ich] - imcb_s[ich] != _nmcb[ich] ) {
      if (debug) cout << "assert imcb_e[ich] - imcb_s[ich] == _nmcb[ich] failed for event " << ient;
      if (debug) cout << ": imcb_e["<<ich<<"]= " << imcb_e[ich] << " imcb_s["<<ich<<"]= " << imcb_s[ich] << " _nmcb["<<ich<<"]= " <<  _nmcb[ich] << endl;
      ok=false;
    }
    _nmcb_sum+=_nmcb[ich];
    if (isb_e[ich] - isb_s[ich] != _nsb[ich] ) {
      if (debug) cout << "assert isb_e[ich] - isb_s[ich] == _nsb[ich] failed for event " << ient;
      if (debug) cout << ": isb_e["<<ich<<"]= " << isb_e[ich] << " isb_s["<<ich<<"]= " << isb_s[ich] << " _nsb["<<ich<<"]= " <<  _nsb[ich] << endl;
      ok=false;
    }
    _nsb_sum+=_nsb[ich];
    if (ipb_e[ich] - ipb_s[ich] != _npb[ich] ) {
      if (debug) cout << "assert ipb_e[ich] - ipb_s[ich] == _npb[ich] failed for event " << ient;
      if (debug) cout << ": ipb_e["<<ich<<"]= " << ipb_e[ich] << " ipb_s["<<ich<<"]= " << ipb_s[ich] << " _npb["<<ich<<"]= " <<  _npb[ich] << endl;
      ok=false;
    }
    _npb_sum+=_npb[ich];
  }
  if (nmcb!=_nmcb_sum ) {
    if (debug) cout << "assert nmcb==_nmcb_sum failed for event " << ient;
    if (debug) cout << ": nmcb= " << nmcb << " _nmcb_sum= " << _nmcb_sum << endl;
    ok=false;
  }
  if (nsb!=_nsb_sum ) {
    if (debug) cout << "assert nsb==_nsb_sum failed for event " << ient;
    if (debug) cout << ": nsb= " << nsb << " _nsb_sum= " << _nsb_sum << endl;
    ok=false;
  }
  if (npb!=_npb_sum ) {
    if (debug) cout << "assert npb==_npb_sum failed for event " << ient;
    if (debug) cout << ": npb= " << npb << " _npb_sum= " << _npb_sum << endl;
    ok=false;
  }
  return ok;
}

Bool_t t::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void t::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t t::Cut(Long64_t entry)
{
  if (entry==0)
    std::cout << "Cut on entry= " << entry << std::endl;
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

TH1F* t::th1f(TString name, TString title, int nb, float min, float max, int col) {
  TH1F* tmp = new TH1F(name, title, nb, min, max);
  set_style(tmp,col);
  return tmp;
}

TLegend* t::tleg(double xmin, double ymin, double xmax, double ymax) {
  TLegend* tmp = new TLegend(xmin, ymin, xmax, ymax);
  tmp->SetFillStyle(0);
  tmp->SetBorderSize(0);
  tmp->SetTextSize(0.07);
  return tmp;
}

void t::set_style(TH2* h) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);

  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);

  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  //if (col>0) {
  //  h->SetLineWidth(3);
  //  h->SetLineColor(col);
  //}
}

void t::set_style(TH1* h, int col) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);

  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);

  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}



#endif // #ifdef t_cxx

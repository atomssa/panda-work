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


class t {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Int_t           nch;
  Int_t           charge[10];   //[nch]
  Float_t         mom_mc[10];   //[nch]
  Float_t         mom_rec[10];   //[nch]
  Float_t         mom_cor[10];   //[nch]
  Float_t         mom_wcor[10];   //[nch]
  Float_t         mom_mrg[10];   //[nch]
  Float_t         mom_sep[10];   //[nch]
  Float_t         mom_stored[10];   //[nch]
  Float_t         phi[10];   //[nch]
  Float_t         the[10];   //[nch]
  Int_t           nphot_sep[10];   //[nch]
  Int_t           nphot_mrg[10];   //[nch]
  Int_t           is_prim[10];   //[nch]
  Int_t           _nmcb[10];   //[nch]
  Int_t           imcb_s[10];   //[nch]
  Int_t           imcb_e[10];   //[nch]
  Int_t           nmcb;
  Float_t         mcb_phi[200];   //[nmcb]
  Float_t         mcb_the[200];   //[nmcb]
  Float_t         mcb_ene[200];   //[nmcb]
  Float_t         mcb_rad[200];   //[nmcb]
  Float_t         mcb_zed[200];   //[nmcb]
  Int_t           mcb_match[200];   //[nmcb]
  Int_t           mcb_score[200];   //[nmcb]
  Int_t           mcb_match_ab[200];   //[nmcb]
  Int_t           mcb_score_ab[200];   //[nmcb]
  Int_t           _nsb[10];   //[nch]
  Int_t           isb_s[10];   //[nch]
  Int_t           isb_e[10];   //[nch]
  Int_t           nsb;
  Int_t           sb_idx[100];   //[nsb]
  Float_t         sb_phi[100];   //[nsb]
  Float_t         sb_the[100];   //[nsb]
  Float_t         sb_ene[100];   //[nsb]
  Float_t         sb_rcalc[100];   //[nsb]
  Int_t           sb_match[100];   //[nsb]
  Int_t           sb_score[100];   //[nsb]
  Int_t           nab;
  Float_t         ab_phi[100];   //[nsb]
  Float_t         ab_the[100];   //[nsb]
  Float_t         ab_ene[100];   //[nsb]
  Int_t           ab_isb[100];   //[nsb]
  Int_t           ab_ich[100];   //[nsb]
  Int_t           ab_match[100];   //[nsb]
  Int_t           ab_score[100];   //[nsb]

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
  TBranch        *b_phi;   //!
  TBranch        *b_the;   //!
  TBranch        *b_nphot_sep;   //!
  TBranch        *b_nphot_mrg;   //!
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
  TBranch        *b_sb_match;   //!
  TBranch        *b_sb_score;   //!
  TBranch        *b_nab;
  TBranch        *b_ab_phi;   //[nsb]
  TBranch        *b_ab_the;   //[nsb]
  TBranch        *b_ab_ene;   //[nsb]
  TBranch        *b_ab_isb;   //[nsb]
  TBranch        *b_ab_ich;   //[nsb]
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
  fChain->SetBranchAddress("phi", phi, &b_phi);
  fChain->SetBranchAddress("the", the, &b_the);
  fChain->SetBranchAddress("nphot_sep", nphot_sep, &b_nphot_sep);
  fChain->SetBranchAddress("nphot_mrg", nphot_mrg, &b_nphot_mrg);
  fChain->SetBranchAddress("is_prim", is_prim, &b_is_prim);
  fChain->SetBranchAddress("_nmcb", _nmcb, &b__nmcb);
  fChain->SetBranchAddress("imcb_s", imcb_s, &b_imcb_s);
  fChain->SetBranchAddress("imcb_e", imcb_e, &b_imcb_e);
  fChain->SetBranchAddress("nmcb", &nmcb, &b_nmcb);
  fChain->SetBranchAddress("mcb_phi", mcb_phi, &b_mcb_phi);
  fChain->SetBranchAddress("mcb_the", mcb_the, &b_mcb_the);
  fChain->SetBranchAddress("mcb_ene", mcb_ene, &b_mcb_ene);
  fChain->SetBranchAddress("mcb_zed", mcb_zed, &b_mcb_zed);
  fChain->SetBranchAddress("mcb_rad", mcb_rad, &b_mcb_rad);
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
  fChain->SetBranchAddress("sb_match", sb_match, &b_sb_match);
  fChain->SetBranchAddress("sb_score", sb_score, &b_sb_score);
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

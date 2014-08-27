#include <TH1F.h>
#include "TLorentzVector.h"
#include <TH2F.h>
#include "TNamed.h"
#include <vector>
#include <string>
#include <cassert>
#include <iostream>

class TH1F;
class TLorentzVector;

using std::string;
using std::vector;
using std::cout;
using std::endl;

class filler {
 protected:
  TNamed *hist;
  vector<int> part_idx;
 public:

 filler(const char* h_name, const char* h_title, const int &npart): hist(new TNamed(h_name,h_title)), part_idx(npart){}

  ~filler() { delete hist; }

  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TLorentzVector&)=0;

};

class filler1d: filler {
 public:
 filler1d(const char* h_name, const char* h_title, const int &npart, const int &nbins, const float &min, const float &max ): filler(h_name,h_title,npart)
  {
    TH1F *htemp = dynamic_cast<TH1F*>(hist);
    htemp->SetBins(nbins, min, max);
  }
  /*
 filler1d(string h_name, string h_title, int npart, const int &nbins, const float &min, const float &max ):
  hist(),part_idx(npart)
  {
    hist.SetNameTitle(h_name.c_str(),h_title.c_str());
    hist.SetBins(nbins, min, max);
  }

 filler1d(string h_name, string h_title, const vector<int> &_part_idx, const int &nbins, const float &min, const float &max ):
  hist(),part_idx(_part_idx.size())
    {
      // constraints
      assert(part_idx.size()>0);
      // initialization of members
      part_idx = _part_idx; // copy ctor
      hist.SetNameTitle(h_name.c_str(),h_title.c_str());
      hist.SetBins(nbins, min, max);
    }
  */
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TLorentzVector&)=0;

  TH1F* getHist() {return dynamic_cast<TH1F*>(hist); }

};

class mom_filler1d: public filler1d {
 public:

 //mom_filler1d(vector<string> part_name, const vector<int> &_part_idx, const int &nbins, const float &min, const float &max):
 // filler1d(Form("mom_%s",part_name[_part_idx[0]].c_str()),
 //	   Form("mom_%s",part_name[_part_idx[0]].c_str()),
 //	   _part_idx, nbins, min, max) {
 //   assert( _part_idx.size()==1 );
 //   assert( part_name.size()==_part_idx.size() );
 // }

 mom_filler1d(const int &_part_idx, const char *part_name, const int &nbins, const float &min, const float &max):
  filler1d(Form("mom_%s",part_name), Form("mom_%s",part_name), 1, nbins, min, max) {
    part_idx[0] = _part_idx;
  }

  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=part_idx.size());
    TH1F* htemp = dynamic_cast<TH1F*> hist;
    htemp->Fill(p4s[part_idx[0]].Vect().Mag());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TLorentzVector &boost) {
    assert(p4s.size()>=part_idx.size());
    TH1F* htemp = dynamic_cast<TH1F*> hist;
    htemp.Fill(p4s[part_idx[0]].Vect().Mag());
  }
};


class mass_filler1d: public filler1d {
 public:

 //mass_filler1d(vector<string> part_name, const vector<int> &_part_idx, const int &nbins, const float &min, const float &max):
 // filler1d(Form("mass_%s_%s",part_name[_part_idx[0]].c_str(),part_name[_part_idx[1]].c_str()),
 //	   Form("mass_%s_%s",part_name[_part_idx[0]].c_str(),part_name[_part_idx[1]].c_str()),
 //	   _part_idx, nbins, min, max) {
 //   assert( _part_idx.size() == 2 );
 //   assert( part_name.size() == _part_idx.size() );
 // }

 mass_filler1d(const int &part0_idx, const int &part1_idx, const char* part1_name, const char *part2_name, const int &nbins, const float &min, const float &max):
  filler1d(Form("mass_%s_%s",part1_name,part2_name), Form("mass_%s_%s",part2_name,part2_name),2, nbins, min, max) {
    part_idx[0] = part0_idx;
    part_idx[1] = part1_idx;
  }

  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=part_idx.size());
    TLorentzVector pair = p4s[0] + p4s[1];
    TH1F* htemp = dynamic_cast<TH1F*> hist;
    htemp->Fill(pair.M());
  }
  // This doesn't make sense for mass, it should be uncallable if possible
  virtual void operator()(const vector<TLorentzVector> &p4s, const TLorentzVector &boost) {

    cout << "Mass with lorentz boost == mass without :P " << endl;
  }

};

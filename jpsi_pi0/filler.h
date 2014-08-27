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

 filler(const int &npart) : part_idx(npart) {
    cout << "begin filler ctor hist pointer= " << hist << " npart= " << npart << endl;
  }

  ~filler() { delete hist; }

  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TLorentzVector&)=0;

  void Write() { this->hist->Write(); }

};

class filler1d: public filler {

 public:
 filler1d(const char* h_name, const char* h_title, const int &npart,
	  const int &nbins, const float &min, const float &max ): filler(npart)
  {
    cout << "begin filler1d ctor hist pointer before new = " << hist << endl;
    hist = new TH1F(h_name, h_title, nbins, min, max);
    cout << "begin filler1d ctor hist pointer after new = " << hist << endl;
    ((TH1*)hist)->SetBins(nbins, min, max);
    cout << "oops" << endl;
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

  TH1* getHist() {return dynamic_cast<TH1*>(hist); }

};


class filler2d: public filler {
 public:
 filler2d(const char* h_name, const char* h_title, const int &npart,
	  const int &nbinsx, const float &xmin, const float &xmax,
	  const int &nbinsy, const float &ymin, const float &ymax): filler(npart)
  {
    cout << "begin filler1d ctor hist pointer before new = " << hist << endl;
    hist = new TH2F(h_name, h_title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    cout << "begin filler1d ctor hist pointer after new = " << hist << endl;
    ((TH2*)hist)->SetBins(nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    cout << "oops" << endl;
  }
  /*
 filler2d(string h_name, string h_title, int npart, const int &nbins, const float &min, const float &max ):
  hist(),part_idx(npart)
  {
    hist.SetNameTitle(h_name.c_str(),h_title.c_str());
    hist.SetBins(nbins, min, max);
  }

 filler2d(string h_name, string h_title, const vector<int> &_part_idx, const int &nbins, const float &min, const float &max ):
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

  TH2* getHist() {return dynamic_cast<TH2*>(hist); }

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
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[part_idx[0]].Vect().Mag());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TLorentzVector &boost) {
    assert(p4s.size()>=part_idx.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[part_idx[0]].Vect().Mag());
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
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(pair.M());
  }
  // This doesn't make sense for mass, it should be uncallable if possible
  virtual void operator()(const vector<TLorentzVector> &p4s, const TLorentzVector &boost) {

    cout << "Mass with lorentz boost == mass without :P " << endl;
  }

};

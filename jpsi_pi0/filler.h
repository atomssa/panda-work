#include <TH1F.h>
#include "TLorentzVector.h"
#include "TVector3.h"
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
using std::move;

//____________
class filler {
 protected:
  TNamed *hist;
  vector<int> pi;
 public:
 filler(const int &npart) : pi(npart) {
    cout << "begin filler ctor hist pointer= " << hist << " npart= " << npart << endl;
  }
  ~filler() { }
  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
  void Write() { this->hist->Write(); }
  void Write(const char* prefix) {
    this->hist->SetName( Form("%s_%s", prefix, hist->GetName()) );
    this->hist->Write();
  }
  TLorentzVector&& boost_transf(const TLorentzVector &vect_in, TVector3 boost) {
    TLorentzVector vect_out(vect_in);
    vect_out.Boost(boost);
    return move(vect_out);
  }
};

//_____________________________
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
  ~filler1d() {delete hist; }
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
  TH1* getHist() {return dynamic_cast<TH1*>(hist); }
};

//_____________________________
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
  ~filler2d(){delete hist; }
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
  TH2* getHist() {return dynamic_cast<TH2*>(hist); }
};

class mom_filler1d: public filler1d {
 public:
 mom_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f="lab", const char* t[]=nullptr):
  filler1d(Form("mom_%s",names[_i]),
	   t!=nullptr?Form("%s momentum (%s frame);p_{%s}[GeV/c^2]",t[_i],f,t[_i]):Form("mom_%s",names[_i]),
	   1, nbins, min, max) {
    pi[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[pi[0]].Vect().Mag());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[pi[0]], boost).Vect().Mag() );
  }
};

//___________________________________
class the_filler1d: public filler1d {
 public:
 the_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f="lab", const char* t[]=nullptr):
  filler1d(Form("the_%s",names[_i]),
	   t!=nullptr?Form("%s #theta (%s frame);#theta_{%s}[rad]",t[_i],f,t[_i]):Form("the_%s",names[_i]),
	   1, nbins, min, max) {
    pi[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[pi[0]].Vect().Theta());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[pi[0]], boost).Vect().Theta() );
  }
};

//___________________________________
class phi_filler1d: public filler1d {
 public:
 phi_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f="lab", const char* t[]=nullptr):
  filler1d(Form("phi_%s",names[_i]),
	   t!=nullptr?Form("%s #phi (%s frame);#phi_{%s}[rad]",t[_i],f,t[_i]):Form("phi_%s",names[_i]),
	   1, nbins, min, max) {
    pi[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[pi[0]].Vect().Phi());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[pi[0]], boost).Vect().Phi() );
  }
};

//_______________________________________
class pair_mass_filler1d: public filler1d {
 public:
 pair_mass_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f, const char* t[]=nullptr):
  filler1d(Form("mass_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair invariant mass;M^{inv}_{%s-%s}",t[_i0],t[_i1],t[_i0],t[_i1]):Form("mass_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]]+p4s[pi[1]]).M());
  }
  // This doesn't make sense for mass, it should be uncallable if possible
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    cout << "Mass with lorentz boost == mass without :P " << endl;
  }
};

//_________________________________________
class pair_the_filler1d: public filler1d {
 public:
 pair_the_filler1d(const int &_i0, const int &_i1, const char *names[], const int &nbins, const float &min, const float &max, const char *f, const char* t[]=nullptr):
  filler1d(Form("p_the_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair polar angle #theta (%s frame);#theta_{%s-%s}",t[_i0],t[_i1],f,t[_i0],t[_i1]):Form("the_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]] + p4s[pi[0]]).Vect().Theta() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf( ( p4s[pi[0]] + p4s[pi[1]] ), boost).Vect().Theta() );
  }
};

//_________________________________________
class pair_phi_filler1d: public filler1d {
 public:
 pair_phi_filler1d(const int &_i0, const int &_i1, const char *names[], const int &nbins, const float &min, const float &max, const char *f, const char* t[]=nullptr):
  filler1d(Form("p_phi_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair azimutal angle #phi (%s frame);#phi_{%s-%s}",t[_i0],t[_i1],f,t[_i0],t[_i1]):Form("the_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]] + p4s[pi[0]]).Vect().Phi() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf( ( p4s[pi[0]] + p4s[pi[1]] ), boost).Vect().Phi() );
  }
};

//_________________________________________
class pair_mom_filler1d: public filler1d {
 public:
 pair_mom_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f, const char* t[]):
  filler1d(Form("p_mom_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair momentum (%s frame);p_{%s-%s}",t[_i0],t[_i1],f,t[_i0],t[_i1]):Form("p_mom_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( ( p4s[pi[0]] + p4s[pi[0]] ).Vect().Mag() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf( ( p4s[pi[0]] + p4s[pi[1]] ), boost).Vect().Mag() );
  }
};

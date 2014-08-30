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
 filler(const int &npart) : pi(npart) {  }
  ~filler() { }
  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
  void Write() { this->hist->Write(); }
  void Write(const char* prefix) {
    this->hist->SetName( Form("%s_%s", prefix, hist->GetName()) );
    this->hist->Write();
  }
  TLorentzVector boost_transf(const TLorentzVector &vect_in, const TVector3 &boost) const {
    TLorentzVector vect_out(vect_in);
    vect_out.Boost(boost);
    return vect_out;
  }

  double mass(const TLorentzVector &v1, const TLorentzVector &v2) const { return (v1+v2).M(); }
  double the(const TLorentzVector &v) const { return v.Vect().Theta(); }
  double the(const TLorentzVector &v, TVector3 boost) const { return boost_transf(v,boost).Vect().Theta(); }

};

//_____________________________
class filler1d: public filler {
 public:
 filler1d(const char* h_name, const char* h_title, const int &npart,
	  const int &nbins, const float &min, const float &max ): filler(npart)
  {
    hist = new TH1F(h_name, h_title, nbins, min, max);
    ((TH1*)hist)->SetBins(nbins, min, max);
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
    hist = new TH2F(h_name, h_title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    ((TH2*)hist)->SetBins(nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  }
  ~filler2d(){delete hist; }
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
  TH2* getHist() {return dynamic_cast<TH2*>(hist); }
};

//____________________________________
class mom_filler1d: public filler1d {
 public:
 mom_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_mom_%s",f[0], names[_i]),
	   t!=nullptr?Form("%s momentum (%s frame);p_{%s}[GeV/c]",t[_i],f[1],t[_i]):Form("mom_%s",names[_i]),
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

//____________________________________
class pt_filler1d: public filler1d {
 public:
 pt_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_pt_%s",f[0], names[_i]),
	   t!=nullptr?Form("%s transverse momentum (%s frame);p_{%s}[GeV/c]",t[_i],f[1],t[_i]):Form("mom_%s",names[_i]),
	   1, nbins, min, max) {
    pi[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[pi[0]].Pt());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[pi[0]], boost).Pt(boost) );
  }
};

//___________________________________
class energy_filler1d: public filler1d {
 public:
 energy_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_e_%s",f[0], names[_i]),
	   t!=nullptr?Form("%s energy (%s frame);p_{%s}[GeV]",t[_i],f[1],t[_i]):Form("energy_%s",names[_i]),
	   1, nbins, min, max) {
    pi[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[pi[0]].E());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[pi[0]], boost).E() );
  }
};

//___________________________________
class the_filler1d: public filler1d {
 public:
 the_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_the_%s",f[0],names[_i]),
	   t!=nullptr?Form("%s #theta (%s frame);#theta_{%s}[rad]",t[_i],f[1],t[_i]):Form("the_%s",names[_i]),
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
    //htemp->Fill( boost_transf(p4s[pi[0]], boost).Vect().Theta() );
    htemp->Fill( boost_transf(p4s[pi[0]], boost).Vect().Angle(boost) );
  }
};

//___________________________________
class cost_filler1d: public filler1d {
 public:
 cost_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_cost_%s",f[0],names[_i]),
	   t!=nullptr?Form("%s cos(#theta_{%s}) (%s frame);cos(#theta_{%s})",t[_i],t[_i],f[1],t[_i]):Form("cost_%s",names[_i]),
	   1, nbins, min, max) {
    pi[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[pi[0]].Vect().CosTheta());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    //htemp->Fill( boost_transf(p4s[pi[0]], boost).Vect().Theta() );
    htemp->Fill( TMath::Cos( boost_transf(p4s[pi[0]], boost).Vect().Angle(boost) ));
  }
};

//___________________________________
class phi_filler1d: public filler1d {
 public:
 phi_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_phi_%s",f[0],names[_i]),
	   t!=nullptr?Form("%s #phi (%s frame);#phi_{%s}[rad]",t[_i],f[1],t[_i]):Form("phi_%s",names[_i]),
	   1, nbins, min, max) {
    pi[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[pi[0]].Vect().Phi());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    //assert(p4s.size()>=pi.size());
    //TH1* htemp = dynamic_cast<TH1*>(hist);
    //htemp->Fill( boost_transf(p4s[pi[0]], boost).Vect().Phi() );
    cout << "Phi with lorentz boost in the frame of boost doesn't make sense i think... " << endl;
    cout << "   because the boost has only a direction, not a frame wrt which phi can be calculated " << endl;
    cout << "      so all the opening angle goes to theta " << endl;
  }
};

//_______________________________________
class pair_mass_filler1d: public filler1d {
 public:
 pair_mass_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr):
  filler1d(Form("mass_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair invariant mass;M^{inv}_{%s-%s}[GeV/c^{2}]",t[_i0],t[_i1],t[_i0],t[_i1]):Form("mass_%s_%s",names[_i0],names[_i1]),
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

//_______________________________________
class mand_s_filler1d: public filler1d {
 public:
 mand_s_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr):
  filler1d(Form("s_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("Mandelstam s (p_{%s}+p_{%s})^{2};s[(GeV/c)^{2}]",t[_i1],t[_i0]):Form("s_%s_%s",names[_i1],names[_i0]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]]+p4s[pi[1]]).M2() );
  }
  // This doesn't make sense for mass, it should be uncallable if possible
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    cout << "Mandelstam variable with lorentz boost == without :P " << endl;
  }
};

//_______________________________________
class mand_t_filler1d: public filler1d {
 public:
 mand_t_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr):
  filler1d(Form("t_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("Mandelstam t (p_{%s}-p_{%s})^{2};t[(GeV/c)^{2}]",t[_i1],t[_i0]):Form("t_%s_%s",names[_i1],names[_i0]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]]-p4s[pi[1]]).M2() );
  }
  // This doesn't make sense for mass, it should be uncallable if possible
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    cout << "Mandelstam variable with lorentz boost == without :P " << endl;
  }
};

//_______________________________________
class mand_u_filler1d: public filler1d {
 public:
 mand_u_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr):
  filler1d(Form("u_%s_%s",names[_i0],names[_i1]),
	   t!=nullptr?Form("Mandelstam u (p_{%s}-p_{%s})^{2};u[(GeV/c)^{2}]",t[_i1],t[_i0]):Form("t_%s_%s",names[_i1],names[_i0]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]]-p4s[pi[1]]).M2() );
  }
  // This doesn't make sense for mass, it should be uncallable if possible
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    cout << "Mandelstam variable with lorentz boost == without :P " << endl;
  }
};

//_________________________________________
class pair_the_filler1d: public filler1d {
 public:
 pair_the_filler1d(const int &_i0, const int &_i1, const char *names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr):
  filler1d(Form("%s_p_the_%s_%s",f[0],names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair polar angle #theta (%s frame);#theta_{%s-%s}[rad]",t[_i0],t[_i1],f[1],t[_i0],t[_i1]):Form("the_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]] + p4s[pi[1]]).Vect().Theta() );
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
 pair_phi_filler1d(const int &_i0, const int &_i1, const char *names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr):
  filler1d(Form("%s_p_phi_%s_%s",f[0],names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair azimutal angle #phi (%s frame);#phi_{%s-%s}[rad]",t[_i0],t[_i1],f[1],t[_i0],t[_i1]):Form("the_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]] + p4s[pi[1]]).Vect().Phi() );
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
 pair_mom_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]):
  filler1d(Form("%s_p_mom_%s_%s",f[0],names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair momentum (%s frame);p_{%s-%s}[GeV/c]",t[_i0],t[_i1],f[1],t[_i0],t[_i1]):Form("p_mom_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( ( p4s[pi[0]] + p4s[pi[1]] ).Vect().Mag() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf( ( p4s[pi[0]] + p4s[pi[1]] ), boost).Vect().Mag() );
  }
};

//_______________________________________
class pair_oa_filler1d: public filler1d {
 public:
 pair_oa_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]):
  filler1d(Form("%s_p_oa_%s_%s",f[0],names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair opening angle (%s frame);OA_{%s-%s}[rad]",t[_i0],t[_i1],f[1],t[_i0],t[_i1]):Form("p_oa_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( p4s[pi[0]].Vect().Angle( p4s[pi[1]].Vect() )  );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    TLorentzVector bp4_0 = boost_transf(p4s[pi[0]], boost);
    TLorentzVector bp4_1 = boost_transf(p4s[pi[1]], boost);
    htemp->Fill( bp4_0.Vect().Angle( bp4_1.Vect() ) );
  }
};

//_______________________________________
class dalitz_filler2d: public filler2d {
 public:
 dalitz_filler2d(const int &_i0, const int &_i1, const int &_i2, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]):
  filler2d(Form("dalitz_%s_%s_%s",names[_i0],names[_i1],names[_i2]),
	   t!=nullptr?Form("%s-%s-%s dalitz; M^{inv}_{%s-%s}[GeV/c^{2}]; M^{inv}_{%s-%s}[GeV/c^{2}]",t[_i0],t[_i1],t[_i2],t[_i0],t[_i1],t[_i1],t[_i2]):
	   Form("dalitz_%s_%s_%s",names[_i0],names[_i1],names[_i2]),
	   3, nbins, min, max, nbins, min, max) {
    pi[0] = _i0;
    pi[1] = _i1;
    pi[2] = _i2;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=pi.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[pi[0]]+p4s[pi[1]]).M(), (p4s[pi[1]]+p4s[pi[2]]).M() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    cout << "Dalitz plot with lorentz boost == dalitz plot without :P " << endl;
  }
};

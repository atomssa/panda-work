#include <TH1F.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TH2F.h>
#include "TNamed.h"
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <iostream>
#include <boost/assert.hpp>

class TH1F;
class TLorentzVector;

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;
using std::map;
using std::make_pair;

struct binning {
binning(int _nbins, double _min, double _max ):nbins(_nbins), min(_min), max(_max) {}
  int nbins;
  double min;
  double max;
};

//____________
class filler {
 protected:

  typedef double(filler::*func)(const TLorentzVector&) const;
  typedef double(filler::*func_boost)(const TLorentzVector&, const TVector3&) const;
  typedef double(filler::*pair_func)(const TLorentzVector&, const TLorentzVector&) const;
  typedef double(filler::*pair_func_boost)(const TLorentzVector&, const TLorentzVector&, const TVector3&) const;

  map<const char*, const char*> vart; // variable tex-form for titles
  map<const char*, const char*> varst; // variable short tex-form for axes
  map<const char*, const char*> varu; // variable units

  map<const char*, func> func_dict;
  map<const char*, func_boost> func_boost_dict;
  map<const char*, pair_func> pair_func_dict;
  map<const char*, pair_func_boost> pair_func_boost_dict;

  TNamed *hist;
  vector<int> ix;
  vector<int> iy;

  vector<func> _func;
  vector<func_boost> _func_b;
  vector<pair_func> _func_p;
  vector<pair_func_boost> _func_p_b;

  //vector<TVector3> dummy;
  //vector<double> dummy_d;
  //vector<int> dummy_i;

  void set_name(const char*name) { hist->SetName(name); }
  void set_title(const char*title) { hist->SetTitle(title); }

 public:

  void init_var_titles( ) {
    vart.insert(make_pair("mass", "Mass")); varst.insert(make_pair("mass", "M^{inv}")); varu.insert(make_pair("mass", "[GeV/c^{2}]"));
    vart.insert(make_pair("mom", "Momentum")); varst.insert(make_pair("mom", "p")); varu.insert(make_pair("mom", "[GeV/c]"));
    vart.insert(make_pair("pt", "p_{T}")); varst.insert(make_pair("pt", "p_{T}")); varu.insert(make_pair("pt", "[GeV/c]"));
    vart.insert(make_pair("e", "Energy")); varst.insert(make_pair("e", "E")); varu.insert(make_pair("e", "[GeV]"));
    vart.insert(make_pair("the", "#theta")); varst.insert(make_pair("the", "#theta")); varu.insert(make_pair("the", "[rad]"));
    vart.insert(make_pair("cost", "cos(#theta)")); varst.insert(make_pair("cost", "cos(#theta)")); varu.insert(make_pair("cost", ""));
    vart.insert(make_pair("phi", "#phi")); varst.insert(make_pair("phi", "#phi")); varu.insert(make_pair("phi", "[rad]"));
    vart.insert(make_pair("oa", "Opening Angle")); varst.insert(make_pair("oa", "OA")); varu.insert(make_pair("oa", "[rad]"));
    vart.insert(make_pair("u", "Mandelstam u")); varst.insert(make_pair("u", "u")); varu.insert(make_pair("u", "[(GeV/c^{2})^{2}]"));
    vart.insert(make_pair("s", "Mandelstam s")); varst.insert(make_pair("s", "s")); varu.insert(make_pair("s", "[(GeV/c^{2})^{2}]"));
    vart.insert(make_pair("t", "Mandelstam t")); varst.insert(make_pair("t", "t")); varu.insert(make_pair("t", "[(GeV/c^{2})^{2}]"));
  }

  void init_pair_func_dict() {
    // invariants
    pair_func_dict.insert(make_pair("mass", &filler::mass));
    pair_func_dict.insert(make_pair("mass_sq", &filler::mass_sq));
    pair_func_dict.insert(make_pair("u", &filler::mand_u));
    pair_func_dict.insert(make_pair("s", &filler::mand_s));
    pair_func_dict.insert(make_pair("t", &filler::mand_t));

    // doesn't make sense in boosted frame where x and y are arbitrary
    pair_func_dict.insert(make_pair("phi", &filler::phi_p));

    // These have a boosted counterpart xxx_p -> xxx_p_b
    pair_func_dict.insert(make_pair("mom", &filler::mom_p));
    pair_func_dict.insert(make_pair("pt", &filler::pt_p));
    pair_func_dict.insert(make_pair("e", &filler::ene_p));
    pair_func_dict.insert(make_pair("the", &filler::the_p));
    pair_func_dict.insert(make_pair("cost", &filler::cost_p));

    // This makes sense only for pairs and has boosted counterpart
    pair_func_dict.insert(make_pair("oa", &filler::oa_p));
  }

  void init_pair_func_boost_dict() {
    pair_func_boost_dict.insert(make_pair("mom", &filler::mom_p_b));
    pair_func_boost_dict.insert(make_pair("pt", &filler::pt_p_b));
    pair_func_boost_dict.insert(make_pair("e", &filler::ene_p_b));
    pair_func_boost_dict.insert(make_pair("the", &filler::the_p_b));

    pair_func_boost_dict.insert(make_pair("cost", &filler::cost_p_b));

    pair_func_boost_dict.insert(make_pair("oa", &filler::oa_p_b));
  }

  void init_func_dict() {
    func_dict.insert(make_pair("mom", &filler::mom));
    func_dict.insert(make_pair("pt", &filler::pt));
    func_dict.insert(make_pair("e", &filler::ene));
    func_dict.insert(make_pair("the", &filler::the));
    func_dict.insert(make_pair("cost", &filler::cost));

    func_dict.insert(make_pair("phi", &filler::phi));
  }

  void init_func_boost_dict() {
    func_boost_dict.insert(make_pair("mom", &filler::mom_b));
    func_boost_dict.insert(make_pair("pt", &filler::pt_b));
    func_boost_dict.insert(make_pair("e", &filler::ene_b));
    func_boost_dict.insert(make_pair("the", &filler::the_b));
    func_boost_dict.insert(make_pair("cost", &filler::cost_b));
  }

  void init_all() {
    init_var_titles();
    init_func_dict();
    init_func_boost_dict();
    init_pair_func_dict();
    init_pair_func_boost_dict();
  }

  filler (TNamed *_hist, const int &npartx) : hist(_hist), ix(npartx) {
    init_all();
  }

  filler (TNamed *_hist, const int &npartx, const int &nparty)
    : hist(_hist), ix(npartx), iy(nparty), _func(2), _func_b(2), _func_p(2), _func_p_b(2) {
    init_all();
  }

 filler(const int &npart) : ix(npart) {
    init_all();
  }

  ~filler() { }
  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
  // if x axis and y axis need to be boosted by different vecotrs -
  // this function should do nothing and communicate annoyingly if called on a 1d filler
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&, const TVector3&)=0;

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
  double mass_sq(const TLorentzVector &v1, const TLorentzVector &v2) const { return (v1+v2).M2(); }
  double mand_s(const TLorentzVector &v1, const TLorentzVector &v2) const { return (v1+v2).M2(); }
  double mand_u(const TLorentzVector &v1, const TLorentzVector &v2) const { return (v2-v1).M2(); }
  double mand_t(const TLorentzVector &v1, const TLorentzVector &v2) const { return (v2-v1).M2(); }

  double mom(const TLorentzVector &v) const { return v.Vect().Mag(); }
  double mom_b(const TLorentzVector &v, const TVector3 &b) const { return mom( boost_transf(v,b) ); }
  double mom_p(const TLorentzVector &v1, const TLorentzVector &v2) const { return mom(v1+v2); }
  double mom_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const { return mom_b(v1+v2,b); }

  double pt(const TLorentzVector &v) const { return v.Vect().Pt(); }
  double pt_b(const TLorentzVector &v, const TVector3 &boost) const { return pt( boost_transf(v,boost) ); }
  double pt_p(const TLorentzVector &v1, const TLorentzVector &v2) const { return pt(v1+v2); }
  double pt_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const { return pt_b(v1+v2,b); }

  double ene(const TLorentzVector &v) const { return v.E(); }
  double ene_b(const TLorentzVector &v, const TVector3 &boost) const { return ene(boost_transf(v,boost)); }
  double ene_p(const TLorentzVector &v1, const TLorentzVector &v2) const { return ene(v1+v2); }
  double ene_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const { return ene_b(v1+v2,b); }

  double the(const TLorentzVector &v) const { return v.Vect().Theta(); }
  double the_b(const TLorentzVector &v, const TVector3 &boost) const { return boost_transf(v,boost).Angle(boost); }
  double the_p(const TLorentzVector &v1, const TLorentzVector &v2) const { return the(v1+v2); }
  double the_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const { return the_b(v1+v2,b); }

  double phi(const TLorentzVector &v) const { return v.Vect().Phi(); }
  double phi_p(const TLorentzVector &v1, const TLorentzVector &v2) const { return phi(v1+v2); }

  double cost(const TLorentzVector &v) const { return v.Vect().CosTheta(); }
  double cost_b(const TLorentzVector &v, const TVector3 &boost) const { return TMath::Cos(boost_transf(v,boost).Vect().Angle(boost)); }
  double cost_p(const TLorentzVector &v1, const TLorentzVector &v2) const { return cost(v1+v2); }
  double cost_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const { return cost_b(v1+v2,b); }

  double oa_p(const TLorentzVector &v1, const TLorentzVector &v2) const { return v1.Vect().Angle(v2.Vect());  }
  double oa_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const { return boost_transf(v1,b).Vect().Angle( boost_transf(v2,b).Vect()); }

};

//_____________________________
class filler1d: public filler {

 protected:
  func _func;
  func_boost _func_b;
  pair_func _func_p;
  pair_func_boost _func_p_b;
  const char* func_tag;
  const char* filler_tag;

 filler1d(const char* h_name, const char* h_title, const int &npart,
	  const int &nbins, const float &min, const float &max ): filler(npart)
  {
    hist = new TH1F(h_name, h_title, nbins, min, max);
  }

 filler1d(const int &npart): filler(new TH1F(),npart){}

  void set_bins(const int &nbins, const float &min, const float &max) {
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->SetBins(nbins, min, max);
  }

  ~filler1d() { delete hist; }

 public:
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
    TH1* htemp = dynamic_cast<TH1*>(hist);
    if (ix.size() == 1) {
      BOOST_ASSERT_MSG(_func != nullptr, filler_tag);
      htemp->Fill( (this->*_func)(p4s[ix[0]]));
    } else if (ix.size() == 2) {
      BOOST_ASSERT_MSG(_func_p != nullptr, filler_tag);
      htemp->Fill( (this->*_func_p)(p4s[ix[0]], p4s[ix[1]]));
    }
  }

  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3& boost) {
    cout<< filler_tag << " --- filer1d::operator() _func = " << _func << " _func_b= " << _func_b << " _func_p= " << _func_p << " _func_p_b= " << _func_p_b << endl;
    BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
    TH1* htemp = dynamic_cast<TH1*>(hist);
    if (ix.size() == 1 ) {
      BOOST_ASSERT_MSG(_func_b != nullptr, filler_tag);
      htemp->Fill( (this->*_func_b)( p4s[ix[0]], boost));
    } else if (ix.size() == 2 ) {
      BOOST_ASSERT_MSG(_func_p_b != nullptr, filler_tag);
      htemp->Fill( (this->*_func_p_b)( p4s[ix[0]] , p4s[ix[1]], boost ) );
    }
  }

  virtual void operator()(const vector<TLorentzVector>&, const TVector3&, const TVector3&) {
    cout << "Calling operator() on an 1d filler with two boost vectors, doesn't make any sense, check your code!" << endl;
  }

  TH1* getHist() {return dynamic_cast<TH1*>(hist); }

};

//____________________________________
class var1d: public filler1d {
 public:
 var1d(const int &_i, const char *var, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(1) {
    ix[0] = _i;

    func_tag = var;
    filler_tag = Form("Filler func tag= %s", func_tag);
    if (func_dict.find(var) == func_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    if (func_boost_dict.find(var) == func_boost_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    _func = func_dict[var];
    _func_b = func_boost_dict[var];

    cout << filler_tag << " ---- set funcs to _func= " << _func << " _func_b= " << _func_b << endl;
    set_name(Form("%s_%s_%s",f[0], var, names[_i]));
    if (t==nullptr)
      set_title(Form("%s_%s",var,names[_i]));
    else
      set_title(Form("%s %s (%s frame);%s_{%s}%s", t[_i], vart[var], f[1], varst[var], t[_i], varu[var] ));
    set_bins(nbins, min, max);
  }
};

//____________________________________
class inv_var1d: public filler1d {
 public:
 inv_var1d(const int &_i, const char *var, const char* names[], const int &nbins, const float &min, const float &max, const char* t[]=nullptr):
  filler1d(1) {
    ix[0] = _i;

    func_tag = var;
    filler_tag = Form("Filler func tag= %s", func_tag);
    if (func_dict.find(var) == func_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    _func = func_dict[var];

    set_name(Form("%s_%s",var, names[_i]));
    if (t==nullptr)
      set_title(Form("%s_%s",var, names[_i]));
    else
      set_title(Form("%s %s;%s_{%s}%s", t[_i], vart[var], varst[var], varu[var], t[_i]));
    set_bins(nbins, min, max);
  }
};

//_________________________________________
class pair_var1d: public filler1d {
 public:
 pair_var1d(const int &_i0, const int &_i1, const char *var, const char *names[],
		   const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr): filler1d(2) {
    ix[0] = _i0;
    ix[1] = _i1;

    func_tag = var;
    filler_tag = Form("Filler func tag= %s", func_tag);
    if (pair_func_dict.find(var) == pair_func_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    if (pair_func_boost_dict.find(var) == pair_func_boost_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    _func_p = pair_func_dict[var];
    _func_p_b = pair_func_boost_dict[var];

    set_name(Form("%s_p_%s_%s_%s",f[0],var, names[_i0],names[_i1]));
    if (t==nullptr)
      set_title(Form("%s_%s_%s",var, names[_i0],names[_i1]));
    else
      set_title(Form("%s-%s pair %s (%s frame);%s_{%s-%s}%s",t[_i0],t[_i1], vart[var], f[1], varst[var], t[_i0],t[_i1], varu[var]));
    set_bins(nbins, min, max);
  }
};

//_________________________________________
class pair_inv_var1d: public filler1d {
 public:
 pair_inv_var1d(const int &_i0, const int &_i1, const char *var, const char *names[],
		   const int &nbins, const float &min, const float &max, const char* t[]=nullptr): filler1d(2) {
    ix[0] = _i0;
    ix[1] = _i1;

    func_tag = var;
    filler_tag = Form("Filler func tag= %s", func_tag);
    if (pair_func_dict.find(var) == pair_func_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    _func_p = pair_func_dict[var];

    set_name(Form("_%s_%s_%s", var, names[_i0],names[_i1]));
    if (t==nullptr)
      set_title(Form("%s_%s_%s",var, names[_i0],names[_i1]));
    else
      set_title(Form("%s-%s pair %s;%s _{%s-%s}%s", t[_i0], t[_i1], vart[var], varst[var], t[_i0], t[_i1], varu[var] ));
    set_bins(nbins, min, max);
  }
};

////_____________________________
//class filler2d: public filler {
// public:
// filler2d(const char* h_name, const char* h_title, const int &npart,
//	  const int &nbinsx, const float &xmin, const float &xmax,
//	  const int &nbinsy, const float &ymin, const float &ymax): filler(npart)
//  {
//    hist = new TH2F(h_name, h_title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
//    ((TH2*)hist)->SetBins(nbinsx, xmin, xmax, nbinsy, ymin, ymax);
//  }
//  ~filler2d(){delete hist; }
//  // Virtual overlaod function call operators to pass the 4momenta to be filled
//  virtual void operator()(const vector<TLorentzVector>&)=0;
//  // if needs to be boosted by some
//  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
//  TH2* getHist() {return dynamic_cast<TH2*>(hist); }
//};

//_____________________________
class filler2d: public filler {

 protected:

  func wtf_func_x, wtf_func_y;
  func_boost wtf_func_b_x, wtf_func_b_y;
  pair_func wtf_func_p_x, wtf_func_p_y;
  pair_func_boost wtf_func_p_b_x,wtf_func_p_b_y;

  const char *func_tag_x, *func_tag_y;
  const char *filler_tag;

 filler2d(const char* h_name, const char* h_title, const int &npart,
	  const int &nbinsx, const float &xmin, const float &xmax,
	  const int &nbinsy, const float &ymin, const float &ymax): filler(npart) {
    hist = new TH2F(h_name, h_title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  }

 filler2d(const int &npartx, const int &nparty): filler(new TH2F(),npartx,nparty) {  }

  void set_bins(const int &nbinsx, const float &xmin, const float &xmax,
		const int &nbinsy, const float &ymin, const float &ymax) {
    TH2* htemp = dynamic_cast<TH2*>(hist);
    htemp->SetBins(nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  }

  void set_bins(const binning &x, const binning &y) {
    TH2* htemp = dynamic_cast<TH2*>(hist);
    htemp->SetBins(x.nbins, x.min, y.max, y.nbins, y.min, y.max);
  }

  ~filler2d() { delete hist; }

 private:

  double value(func __func, const TLorentzVector &v) {
    cout << "filler2d::value(func, TLV) trying to get a ref to function " << __func << endl;
    BOOST_ASSERT_MSG(__func != nullptr, filler_tag);
    return (this->*__func)(v);
  }

  double value(func_boost __func, const TLorentzVector &v, const TVector3 &b) {
    cout << "filler2d::value(func_boost, TLV, TV3) trying to get a ref to function " << __func << endl;
    BOOST_ASSERT_MSG(__func != nullptr, filler_tag);
    return (this->*__func)(v,b);
  }

  double value(pair_func __func, const TLorentzVector &v1, const TLorentzVector &v2) {
    cout << "filler2d::value(pair_func, TLV1, TLV2) trying to get a ref to function " << __func << endl;
    BOOST_ASSERT_MSG(__func != nullptr, filler_tag);
    return (this->*__func)(v1,v2);
  }

  double value(pair_func_boost __func, const TLorentzVector &v1, const TLorentzVector &v2, const TVector3 &b) {
    cout << "filler2d::value(pair_func_boost, TLV1, TLV2, b) trying to get a ref to function " << __func << endl;
    BOOST_ASSERT_MSG(__func != nullptr, filler_tag);
    return (this->*__func)(v1,v2,b);
  }

  double get_value(const int &iaxis, const vector<int> &ii, const vector<TLorentzVector> &p4s) {
    if (ii.size() == 1) {
      cout << "filler2d::get_value(axisidx=1, part_idx_list, part_idx) " << endl;
      return value(_func[iaxis], p4s[ii[0]]);
      //return iaxis==0?value(wtf_func_x, p4s[ii[0]]):value(wtf_func_y, p4s[ii[0]]);
    } else if (ii.size() == 2) {
      cout << "filler2d::get_value(axisidx=2, part_idx_list, part_idx) " << endl;
      return value(_func_p[iaxis], p4s[ii[0]], p4s[ii[1]]);
      //return iaxis==0?value(wtf_func_p_x, p4s[ii[0]], p4s[ii[1]] ):value(wtf_func_p_y, p4s[ii[0]], p4s[ii[1]]);
    }
    cout << "blah" << endl;
  }

  double get_value(const int &iaxis, const vector<int> &ii, const vector<TLorentzVector> &p4s, const TVector3 &b) {
    if (ii.size() == 1) {
      //cout << "filler2d::get_value(axis=" << iaxis << ", part_idx_list, part_idx, boost_vect) npart=1 calling get_value with _func_b[" << iaxis << "]" << endl;
      //cout << "_func_b[" << iaxis << "] = " << _func_b[iaxis] << endl;
      //cout << "ftag_x= " << func_tag_x << " ftag_y=  " << func_tag_y << " _func.size()= " << _func.size()  << " _func_b.size()= " << _func_b.size()
      // << " _func_p.size()= " << _func_p.size()  << " _func_p_b.size()= " << _func_p_b.size() << endl;
      cout << "filler2d::get_value: "
	//<<"dummy_d.size()= " << dummy_d.size() << " dummy.size()= " << dummy.size() << " dummy_i.size()= " << dummy_i.size()
	   << " ix.size()= " << ix.size() << " iy.size()= " << iy.size()
	//<< " ftag_x= " << func_tag_x << " ftag_y=  " << func_tag_y
      	   << " _func.size()= " << _func.size()  << " _func_b.size()= " << _func_b.size()
      	   << " _func_p.size()= " << _func_p.size()  << " _func_p_b.size()= " << _func_p_b.size() << endl;
      cout << "filler2d::get_value(private) with func wtf_func_b_x= " << wtf_func_b_x << " wtf_Func_b_y= " << wtf_func_b_y << endl;
      //return 0.0;
      return value(_func_b[iaxis], p4s[ii[0]], b);
      //return iaxis==0?value(wtf_func_b_x, p4s[ii[0]], b):value(wtf_func_b_y, p4s[ii[0]], b);
    } else if (ii.size() == 2) {
      cout << "filler2d::get_value(axis=" << iaxis << ", part_idx_list, part_idx, boost_vect) npart=2 calling get_value with _func_b[" << iaxis << "]" << endl;
      return value(_func_p_b[iaxis], p4s[ii[0]], p4s[ii[1]], b);
      //return iaxis==0?value(wtf_func_p_b_x, p4s[ii[0]], p4s[ii[1]], b):value(wtf_func_p_b_y, p4s[ii[0]], p4s[ii[1]], b);
    }
  }

 public:
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    BOOST_ASSERT_MSG(p4s.size()>=ix.size() , filler_tag);
    TH2* htemp = dynamic_cast<TH2*>(hist);
    htemp->Fill(get_value(0,ix,p4s),get_value(1,iy,p4s));
  }

  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3& b) {
    cout << filler_tag << "  operator() (public) with func "
    	 << "wtf_func_x= " << wtf_func_x << " wtf_func_y= " << wtf_func_y
	 << "wtf_func_p_x= " << wtf_func_p_x << " wtf_func_p_y= " << wtf_func_p_y
	 << "wtf_func_b_x= " << wtf_func_b_x << " wtf_Func_b_y= " << wtf_func_b_y
	 << "wtf_func_p_b_x= " << wtf_func_p_b_x << " wtf_func_p_b_y= " << wtf_func_p_b_y
    	 << endl;
    BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
    TH2* htemp = dynamic_cast<TH2*>(hist);
    htemp->Fill(get_value(0,ix,p4s,b), get_value(1,iy,p4s,b));
  }

  virtual void operator()(const vector<TLorentzVector>&p4s, const TVector3&b1, const TVector3&b2) {
    BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
    TH2* htemp = dynamic_cast<TH2*>(hist);
    htemp->Fill(get_value(0,ix,p4s,b1), get_value(1,iy,p4s,b2));
  }

  TH2* getHist() {return dynamic_cast<TH2*>(hist); }

};

//____________________________________
class var2d: public filler2d {
   protected:
  void set_funcs(const int &iaxis, const char* var) {
    if (func_dict.find(var) == func_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    if (func_boost_dict.find(var) == func_boost_dict.end()) { cout << "WARNING: Cant find function with tag " << var << endl; return; }
    //_func.push_back(func_dict[var]);
    _func[iaxis] = func_dict[var];
    _func_b[iaxis] = func_boost_dict[var];

    if (iaxis==0) {
      wtf_func_x = func_dict[var];
      wtf_func_b_x = func_boost_dict[var];
    } else {
      wtf_func_y = func_dict[var];
      wtf_func_b_y = func_boost_dict[var];
    }
    //dummy.push_back(TVector3());
    //dummy_d.push_back(0.0);

    cout << "var2d::set_funcs     "
      //<< "dummy_d.size()= " << dummy_d.size() << " dummy.size()= " << dummy.size() << " dummy_i.size()= " << dummy_i.size()
  	 << " ix.size()= " << ix.size() << " iy.size()= " << iy.size()
      //<< " ftag_x= " << func_tag_x << " ftag_y=  " << func_tag_y
  	 << " _func.size()= " << _func.size()  << " _func_b.size()= " << _func_b.size()
  	 << " _func_p.size()= " << _func_p.size()  << " _func_p_b.size()= " << _func_p_b.size() << endl;
  }
  void register_vars(const char* varx, const char* vary) {
    cout <<"registering vars" << endl;
    func_tag_x = varx;
    func_tag_y = vary;
    set_funcs(0,varx); // make sure x is called before y (otherwise axes will be inversed)
    set_funcs(1,vary); // make sure x is called before y (otherwise axes will be inversed)
  }

 public:

 var2d(const int &_ix, const char *varx, const int &nbinsx, const double &xmin, const double &xmax, const char* framex[],
       const int &_iy, const char *vary, const int &nbinsy, const double &ymin, const double &ymax, const char* framey[],
       const char* names[], const char* t[]=nullptr): filler2d(1,1) {
    filler_tag = Form("Filler func tagx = %s tagy = %s", varx, vary);
    ix[0] = _ix;
    iy[0] = _iy;
    register_vars(varx, vary);
    set_name(Form("%s_%s_%s_%s_%s_%s",framex[0], varx, names[_ix],framey[0], vary, names[_iy]));
    if (t==nullptr) {
      set_title(hist->GetName());
    } else {
      set_title(Form("%s %s%s vs %s %s%s", t[_iy], vart[vary], framey[2], t[_ix], vart[varx], framex[2] ));
    }
    set_bins(nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  }

 var2d(const int &_ix, const char *varx, const binning &x, const char* framex[],
       const int &_iy, const char *vary, const binning &y, const char* framey[],
       const char* names[], const char* t[]=nullptr): filler2d(1,1){

    //var2d(_ix,varx,x.nbins,x.min,x.max,framex,_iy,vary,y.nbins,y.min,y.max,framey,names,t); -- why didn't this work?
    filler_tag = Form("Filler func tagx = %s tagy = %s", varx, vary);
    ix[0] = _ix;
    iy[0] = _iy;
    register_vars(varx, vary);
    set_name(Form("%s_%s_%s_%s_%s_%s",framex[0], varx, names[_ix],framey[0], vary, names[_iy]));
    if (t==nullptr) {
      set_title(hist->GetName());
    } else {
      set_title(Form("%s %s%s vs %s %s%s", t[_iy], vart[vary], framey[2], t[_ix], vart[varx], framex[2] ));
    }
    set_bins(x, y);
  }

};














//____________________________________
class mom_filler1d: public filler1d {
 public:
 mom_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_mom_%s",f[0], names[_i]),
	   t!=nullptr?Form("%s momentum (%s frame);p_{%s}[GeV/c]",t[_i],f[1],t[_i]):Form("mom_%s",names[_i]),
	   1, nbins, min, max) {
    ix[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[ix[0]].Vect().Mag());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[ix[0]], boost).Vect().Mag() );
  }
};

//____________________________________
class pt_filler1d: public filler1d {
 public:
 pt_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_pt_%s",f[0], names[_i]),
	   t!=nullptr?Form("%s transverse momentum (%s frame);p_{%s}[GeV/c]",t[_i],f[1],t[_i]):Form("mom_%s",names[_i]),
	   1, nbins, min, max) {
    ix[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[ix[0]].Pt());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[ix[0]], boost).Pt(boost) );
  }
};

//___________________________________
class energy_filler1d: public filler1d {
 public:
 energy_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_e_%s",f[0], names[_i]),
	   t!=nullptr?Form("%s energy (%s frame);p_{%s}[GeV]",t[_i],f[1],t[_i]):Form("energy_%s",names[_i]),
	   1, nbins, min, max) {
    ix[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[ix[0]].E());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf(p4s[ix[0]], boost).E() );
  }
};

//___________________________________
class the_filler1d: public filler1d {
 public:
 the_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_the_%s",f[0],names[_i]),
	   t!=nullptr?Form("%s #theta (%s frame);#theta_{%s}[rad]",t[_i],f[1],t[_i]):Form("the_%s",names[_i]),
	   1, nbins, min, max) {
    ix[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[ix[0]].Vect().Theta());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    //htemp->Fill( boost_transf(p4s[ix[0]], boost).Vect().Theta() );
    htemp->Fill( boost_transf(p4s[ix[0]], boost).Vect().Angle(boost) );
  }
};

//___________________________________
class cost_filler1d: public filler1d {
 public:
 cost_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_cost_%s",f[0],names[_i]),
	   t!=nullptr?Form("%s cos(#theta_{%s}) (%s frame);cos(#theta_{%s})",t[_i],t[_i],f[1],t[_i]):Form("cost_%s",names[_i]),
	   1, nbins, min, max) {
    ix[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[ix[0]].Vect().CosTheta());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    //htemp->Fill( boost_transf(p4s[ix[0]], boost).Vect().Theta() );
    htemp->Fill( TMath::Cos( boost_transf(p4s[ix[0]], boost).Vect().Angle(boost) ));
  }
};

//___________________________________
class phi_filler1d: public filler1d {
 public:
 phi_filler1d(const int &_i, const char* names[], const int &nbins, const float &min, const float &max, const char* f[], const char* t[]=nullptr):
  filler1d(Form("%s_phi_%s",f[0],names[_i]),
	   t!=nullptr?Form("%s #phi (%s frame);#phi_{%s}[rad]",t[_i],f[1],t[_i]):Form("phi_%s",names[_i]),
	   1, nbins, min, max) {
    ix[0] = _i;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill(p4s[ix[0]].Vect().Phi());
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    //assert(p4s.size()>=ix.size());
    //TH1* htemp = dynamic_cast<TH1*>(hist);
    //htemp->Fill( boost_transf(p4s[ix[0]], boost).Vect().Phi() );
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
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[ix[0]]+p4s[ix[1]]).M());
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
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[ix[0]]+p4s[ix[1]]).M2() );
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
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[ix[0]]-p4s[ix[1]]).M2() );
    cout << "org t = " << (p4s[ix[0]]-p4s[ix[1]]).M2() << endl;
    cout << "======================================" << endl;
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
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[ix[0]]-p4s[ix[1]]).M2() );
    cout << "org u = " << (p4s[ix[0]]-p4s[ix[1]]).M2() << endl;
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
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[ix[0]] + p4s[ix[1]]).Vect().Theta() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf( ( p4s[ix[0]] + p4s[ix[1]] ), boost).Vect().Theta() );
  }
};

//_________________________________________
class pair_phi_filler1d: public filler1d {
 public:
 pair_phi_filler1d(const int &_i0, const int &_i1, const char *names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]=nullptr):
  filler1d(Form("%s_p_phi_%s_%s",f[0],names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair azimutal angle #phi (%s frame);#phi_{%s-%s}[rad]",t[_i0],t[_i1],f[1],t[_i0],t[_i1]):Form("the_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[ix[0]] + p4s[ix[1]]).Vect().Phi() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf( ( p4s[ix[0]] + p4s[ix[1]] ), boost).Vect().Phi() );
  }
};

//_________________________________________
class pair_mom_filler1d: public filler1d {
 public:
 pair_mom_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]):
  filler1d(Form("%s_p_mom_%s_%s",f[0],names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair momentum (%s frame);p_{%s-%s}[GeV/c]",t[_i0],t[_i1],f[1],t[_i0],t[_i1]):Form("p_mom_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( ( p4s[ix[0]] + p4s[ix[1]] ).Vect().Mag() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( boost_transf( ( p4s[ix[0]] + p4s[ix[1]] ), boost).Vect().Mag() );
  }
};

//_______________________________________
class pair_oa_filler1d: public filler1d {
 public:
 pair_oa_filler1d(const int &_i0, const int &_i1, const char* names[], const int &nbins, const float &min, const float &max, const char *f[], const char* t[]):
  filler1d(Form("%s_p_oa_%s_%s",f[0],names[_i0],names[_i1]),
	   t!=nullptr?Form("%s-%s pair opening angle (%s frame);OA_{%s-%s}[rad]",t[_i0],t[_i1],f[1],t[_i0],t[_i1]):Form("p_oa_%s_%s",names[_i0],names[_i1]),
	   2, nbins, min, max) {
    ix[0] = _i0;
    ix[1] = _i1;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( p4s[ix[0]].Vect().Angle( p4s[ix[1]].Vect() )  );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    TLorentzVector bp4_0 = boost_transf(p4s[ix[0]], boost);
    TLorentzVector bp4_1 = boost_transf(p4s[ix[1]], boost);
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
    ix[0] = _i0;
    ix[1] = _i1;
    ix[2] = _i2;
  }
  virtual void operator()(const vector<TLorentzVector> &p4s) {
    assert(p4s.size()>=ix.size());
    TH1* htemp = dynamic_cast<TH1*>(hist);
    htemp->Fill( (p4s[ix[0]]+p4s[ix[1]]).M(), (p4s[ix[1]]+p4s[ix[2]]).M() );
  }
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3 &boost) {
    cout << "Dalitz plot with lorentz boost == dalitz plot without :P " << endl;
  }
};

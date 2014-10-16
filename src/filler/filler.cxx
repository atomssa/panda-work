#include "filler.h"
#include <TH1F.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TH2F.h>
#include "TNamed.h"
#include <vector>
#include <map>
#include <cassert>
#include <iostream>
#include <boost/assert.hpp>

class TH1F;
class TLorentzVector;

typedef const TLorentzVector& cr_tlv;
typedef const TVector3& cr_tv3;
typedef const std::vector<TLorentzVector>& cr_vec_tlv;

using std::vector;
using std::cout;
using std::endl;
using std::move;
using std::map;
using std::make_pair;

static const double m_2pi = 2.0*TMath::Pi();
static const double m_pi = TMath::Pi();
static const double rtd= TMath::RadToDeg();

filler::filler(TNamed *_hist, const vector<int> &_ix):
  hist(_hist), ix(_ix), _func(1), _func_b(1), _func_p(1), _func_p_b(1) {
  init_all();
}

filler::filler(TNamed *_hist, const vector<int> &_ix, const vector<int> &_iy):
  hist(_hist), ix(_ix), iy(_iy), _func(2), _func_b(2), _func_p(2), _func_p_b(2) {
  init_all();
}

// LEGACY
filler::filler (TNamed *_hist, const int &npart) : hist(_hist), ix(npart) {
  init_all();
}

// LEGACY
filler::filler(const int &npart) : ix(npart) {
  init_all();
}

filler::~filler() { }

void filler::Write() { this->hist->Write(); }

void filler::Write(const char* prefix) {
  this->hist->SetName( Form("%s_%s", prefix, hist->GetName()) );
  this->hist->Write();
}

void filler::set_name(const char*name) { hist->SetName(name); }
void filler::set_title(const char*title) { hist->SetTitle(title); }

const char* filler::warning_msg(const char* var, const char* dict) {
  const char *line1 = Form( "WARNING: Can not find function with tag %s in %s", var, dict);
  if (strstr(dict,"boost")==nullptr)
    return Form("%s\n If you call this filler without a boost expect a crash ... ", line1);
  else
    return Form("%s\n If you call this filler with a boost, expect a crash ... ", line1);
}

void filler::init_var(const char* var, const char* title, const char *stitle,
	      const char*unit, const axis& bins) {
  vart.insert(make_pair(var, title));
  varst.insert(make_pair(var, stitle));
  varu.insert(make_pair(var, unit));
  varb.insert(make_pair(var, bins));
}

void filler::init_vars( ) {
  const int nbin = 200;
  init_var("mass", "Invariant Mass", "M^{inv}", "[GeV/c^{2}]", axis(nbin, 0.0, 10.0));
  init_var("mass_sq", "Invariant Mass Squared", "(M^{inv})^{2}", "[(GeV/c^{2})^2]", axis(nbin, 0.0, 10.0));
  init_var("mom", "Momentum", "p", "[GeV/c]", axis(nbin, 0.0, 10.0));
  init_var("mom_fwd", "Fwd going momentum", "p", "[GeV/c]", axis(nbin, 0.0, 10.0));
  init_var("mom_bwd", "Bwd going momentum", "p", "[GeV/c]", axis(nbin, 0.0, 10.0));
  init_var("pt", "p_{T}", "p_{T}", "[GeV/c]", axis(nbin, 0.0, 10.0));
  init_var("e", "Energy", "E", "[GeV]", axis(nbin, 0.0, 10.0 ));
  init_var("the", "#theta", "#theta", "[#circ]", axis(nbin, -0.1, m_pi+0.1 ));
  init_var("the_fwd", "Fwd going #theta", "#theta^{fwd}", "[#circ]", axis(nbin, 0.0, 10.0));
  init_var("the_bwd", "Bwd going #theta", "#theta^{bwd}", "[#circ]", axis(nbin, 0.0, 10.0));
  init_var("cost", "cos(#theta)", "cos(#theta)", "", axis(nbin, -1.1, 1.1 ));
  init_var("phi", "#phi", "#phi", "[#circ]", axis(nbin, -m_pi, m_pi ));
  init_var("oa", "Opening Angle", "OA", "[#circ]", axis(nbin, -0.2, m_pi+0.2 ));
  init_var("u", "Mandelstam u", "u", "[(GeV/c^{2})^{2}]", axis(nbin, -20.0, 20.0 ));
  init_var("s", "Mandelstam s", "s", "[(GeV/c^{2})^{2}]", axis(nbin, -20.0, 20.0 ));
  init_var("t", "Mandelstam t", "t", "[(GeV/c^{2})^{2}]", axis(nbin, -20.0, 20.0 ));
}

void filler::init_pair_func_dict() {
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
  pair_func_dict.insert(make_pair("mom_fwd", &filler::mom_p_fwd));
  pair_func_dict.insert(make_pair("mom_bwd", &filler::mom_p_bwd));
  pair_func_dict.insert(make_pair("pt", &filler::pt_p));
  pair_func_dict.insert(make_pair("e", &filler::ene_p));
  pair_func_dict.insert(make_pair("the", &filler::the_p));
  pair_func_dict.insert(make_pair("the_fwd", &filler::the_p_fwd));
  pair_func_dict.insert(make_pair("the_bwd", &filler::the_p_bwd));

  pair_func_dict.insert(make_pair("cost", &filler::cost_p));

  // This makes sense only for pairs and has boosted counterpart
  pair_func_dict.insert(make_pair("oa", &filler::oa_p));
}

void filler::init_pair_func_boost_dict() {
  pair_func_boost_dict.insert(make_pair("mom", &filler::mom_p_b));
  pair_func_boost_dict.insert(make_pair("mom_fwd", &filler::mom_p_fwd_b));
  pair_func_boost_dict.insert(make_pair("mom_bwd", &filler::mom_p_bwd_b));
  pair_func_boost_dict.insert(make_pair("pt", &filler::pt_p_b));
  pair_func_boost_dict.insert(make_pair("e", &filler::ene_p_b));
  pair_func_boost_dict.insert(make_pair("the", &filler::the_p_b));
  pair_func_boost_dict.insert(make_pair("the_fwd", &filler::the_p_fwd_b));
  pair_func_boost_dict.insert(make_pair("the_bwd", &filler::the_p_bwd_b));
  pair_func_boost_dict.insert(make_pair("cost", &filler::cost_p_b));
  pair_func_boost_dict.insert(make_pair("oa", &filler::oa_p_b));
}

void filler::init_func_dict() {
  func_dict.insert(make_pair("mom", &filler::mom));
  func_dict.insert(make_pair("pt", &filler::pt));
  func_dict.insert(make_pair("e", &filler::ene));
  func_dict.insert(make_pair("the", &filler::the));
  func_dict.insert(make_pair("cost", &filler::cost));
  func_dict.insert(make_pair("phi", &filler::phi));
}

void filler::init_func_boost_dict() {
  func_boost_dict.insert(make_pair("mom", &filler::mom_b));
  func_boost_dict.insert(make_pair("pt", &filler::pt_b));
  func_boost_dict.insert(make_pair("e", &filler::ene_b));
  func_boost_dict.insert(make_pair("the", &filler::the_b));
  func_boost_dict.insert(make_pair("cost", &filler::cost_b));
}

void filler::init_all() {
  init_vars();
  init_func_dict();
  init_func_boost_dict();
  init_pair_func_dict();
  init_pair_func_boost_dict();
}

TLorentzVector filler::boost_transf(cr_tlv vect_in, cr_tv3 boost) const {
  TLorentzVector vect_out(vect_in);
  vect_out.Boost(boost);
  return vect_out;
}

void filler::set_funcs(const int &iaxis, const char* var) {
  if (func_dict.find(var) == func_dict.end()) { cout << warning_msg(var, "func_dict") << endl; }
  if (func_boost_dict.find(var) == func_boost_dict.end()) { cout << warning_msg(var, "func_boost_dict") << endl; }
  _func[iaxis] = func_dict[var];
  _func_b[iaxis] = func_boost_dict[var];
}

void filler::set_pair_funcs(const int &iaxis, const char* var) {
  if (pair_func_dict.find(var) == pair_func_dict.end()) { cout << warning_msg(var, "pair_func_dict") << endl; }
  if (pair_func_boost_dict.find(var) == pair_func_boost_dict.end()) { cout << warning_msg(var, "pair_func_boost_dict") << endl; }
  _func_p[iaxis] = pair_func_dict[var];
  _func_p_b[iaxis] = pair_func_boost_dict[var];
}

double filler::value(func __func, cr_tlv v) {
  BOOST_ASSERT_MSG(__func != nullptr, "function pointer not set");
  return (this->*__func)(v);
}

double filler::value(func_boost __func, cr_tlv v, cr_tv3 b) {
  BOOST_ASSERT_MSG(__func != nullptr, "function pointer not set");
  return (this->*__func)(v,b);
}

double filler::value(pair_func __func, cr_tlv v1, cr_tlv v2) {
  BOOST_ASSERT_MSG(__func != nullptr, "function pointer not set");
  return (this->*__func)(v1,v2);
}

double filler::value(pair_func_boost __func, cr_tlv v1, cr_tlv v2, cr_tv3 b) {
  BOOST_ASSERT_MSG(__func != nullptr, "function pointer not set");
  return (this->*__func)(v1,v2,b);
}

double filler::mass(cr_tlv v1, cr_tlv v2) const { return (v1+v2).M(); }
double filler::mass_sq(cr_tlv v1, cr_tlv v2) const { return (v1+v2).M2(); }
double filler::mand_s(cr_tlv v1, cr_tlv v2) const { return (v1+v2).M2(); }
double filler::mand_u(cr_tlv v1, cr_tlv v2) const { return (v2-v1).M2(); }
double filler::mand_t(cr_tlv v1, cr_tlv v2) const { return (v2-v1).M2(); }

double filler::mom(cr_tlv v) const { return v.Vect().Mag(); }
double filler::mom_b(cr_tlv v, cr_tv3 b) const { return mom( boost_transf(v,b) ); }
double filler::mom_p(cr_tlv v1, cr_tlv v2) const { return mom(v1+v2); }
double filler::mom_p_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const { return mom_b(v1+v2,b); }
double filler::mom_p_fwd(cr_tlv v1, cr_tlv v2) const {
  return v1.Vect().Theta()<v2.Vect().Theta()?mom(v1):mom(v2);
}
double filler::mom_p_bwd(cr_tlv v1, cr_tlv v2) const {
  return v1.Vect().Theta()>v2.Vect().Theta()?mom(v1):mom(v2);
}
double filler::mom_p_fwd_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const {
  return v1.Vect().Theta()<v2.Vect().Theta()?mom_b(v1,b):mom_b(v2,b);
}
double filler::mom_p_bwd_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const {
  return v1.Vect().Theta()>v2.Vect().Theta()?mom_b(v1,b):mom_b(v2,b);
}

double filler::pt(cr_tlv v) const { return v.Vect().Pt(); }
double filler::pt_b(cr_tlv v, cr_tv3 boost) const { return pt( boost_transf(v,boost) ); }
double filler::pt_p(cr_tlv v1, cr_tlv v2) const { return pt(v1+v2); }
double filler::pt_p_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const { return pt_b(v1+v2,b); }

double filler::ene(cr_tlv v) const { return v.E(); }
double filler::ene_b(cr_tlv v, cr_tv3 boost) const { return ene(boost_transf(v,boost)); }
double filler::ene_p(cr_tlv v1, cr_tlv v2) const { return ene(v1+v2); }
double filler::ene_p_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const { return ene_b(v1+v2,b); }

double filler::the(cr_tlv v) const { return rtd*(v.Vect().Theta()); }
double filler::the_b(cr_tlv v, cr_tv3 boost) const {
  return rtd*(boost_transf(v,boost).Angle(boost));
}
double filler::the_p(cr_tlv v1, cr_tlv v2) const { return the(v1+v2); }
double filler::the_p_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const { return the_b(v1+v2,b); }
double filler::the_p_fwd(cr_tlv v1, cr_tlv v2) const {
  return v1.Vect().Theta()<v2.Vect().Theta()?the(v1):the(v2);
}
double filler::the_p_bwd(cr_tlv v1, cr_tlv v2) const {
  return v1.Vect().Theta()>v2.Vect().Theta()?the(v1):the(v2);
}
double filler::the_p_fwd_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const {
  return v1.Vect().Theta()<v2.Vect().Theta()?the_b(v1,b):the_b(v2,b);
}
double filler::the_p_bwd_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const {
  return v1.Vect().Theta()>v2.Vect().Theta()?the_b(v1,b):the_b(v2,b);
}

double filler::phi(cr_tlv v) const { return rtd*(v.Vect().Phi()); }
double filler::phi_p(cr_tlv v1, cr_tlv v2) const { return phi(v1+v2); }

double filler::cost(cr_tlv v) const { return v.Vect().CosTheta(); }
double filler::cost_b(cr_tlv v, cr_tv3 boost) const {
  return TMath::Cos(boost_transf(v,boost).Vect().Angle(boost));
}
double filler::cost_p(cr_tlv v1, cr_tlv v2) const { return cost(v1+v2); }
double filler::cost_p_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const { return cost_b(v1+v2,b); }

double filler::oa_p(cr_tlv v1, cr_tlv v2) const {
  return rtd*(v1.Vect().Angle(v2.Vect()));
}
double filler::oa_p_b(cr_tlv v1, cr_tlv v2, cr_tv3 b) const {
  return rtd*(boost_transf(v1,b).Vect().Angle( boost_transf(v2,b).Vect()));
}


const char* filler::_name(const char* var,const char* fname,
  const char* pname){
  return Form("%s_%s_%s",fname, var, pname);
}

const char* filler::_name(const char* var, const char* pname) {
  return Form("%s_%s", var, pname);
}

const char* filler::_name_p(const char* var,const char* fname,
  const char* pname1, const char* pname2) {
  return Form("%s_p_%s_%s_%s",fname, var, pname1, pname2);
}

const char* filler::_name_p(const char* var, const char* pname1,
  const char* pname2) {
  return Form("p_%s_%s_%s", var, pname1, pname2);
}

const char* filler::_title(const char *var, const char *ftitle,
  const char *ptitle) {
  return Form("%s %s%s", ptitle, vart[var], ftitle);
}

const char* filler::_title(const char *var, const char *ptitle) {
  return Form("%s %s", ptitle, vart[var]);
}

const char* filler::_title_p(const char *var, const char *ftitle,
  const char *ptitle1, const char *ptitle2) {
  return Form("%s-%s %s%s", ptitle1, ptitle2, vart[var], ftitle);
}

const char* filler::_title_p(const char *var, const char *ptitle1,
  const char *ptitle2) {
  return Form("%s-%s %s", ptitle1, ptitle2, vart[var]);
}

const char* filler::_atitle(const char *var, const char *ptitle) {
  return Form("%s_{%s}%s", varst[var],ptitle, varu[var]);
}

const char* filler::_atitle_p(const char *var, const char *ptitle1,
  const char *ptitle2) {
  return Form("%s_{%s-%s}%s", varst[var],ptitle1, ptitle2, varu[var]);
}

filler1d::filler1d(const char* h_name, const char* h_title, const int &npart,
		   const int &nbins, const float &min, const float &max ):
  filler(npart) {
  hist = new TH1F(h_name, h_title, nbins, min, max);
}

filler1d::filler1d(const int &npart): filler(new TH1F(),npart){}

filler1d::filler1d(const vector<int> &_ix, const char *var):
  filler(new TH1F(),_ix),func_tag(var),filler_tag(Form("Filler func tag= %s", func_tag))
{ }

filler1d::~filler1d() { delete hist; }

void filler1d::set_bins(const int &nbins, const float &min, const float &max) {
  TH1* htemp = dynamic_cast<TH1*>(hist);
  htemp->SetBins(nbins, min, max);
}

//virtual
void filler1d::operator()(cr_vec_tlv p4s) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
  TH1* htemp = dynamic_cast<TH1*>(hist);
  if (ix.size() == 1) {
    BOOST_ASSERT_MSG(_func[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func[0])(p4s[ix[0]]));
  } else if (ix.size() == 2) {
    BOOST_ASSERT_MSG(_func_p[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func_p[0])(p4s[ix[0]], p4s[ix[1]]));
  }
}

//virtual
void filler1d::operator()(cr_vec_tlv p4s, cr_tv3 boost) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
  TH1* htemp = dynamic_cast<TH1*>(hist);
  if (ix.size() == 1 ) {
    BOOST_ASSERT_MSG(_func_b[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func_b[0])( p4s[ix[0]], boost));
  } else if (ix.size() == 2 ) {
    BOOST_ASSERT_MSG(_func_p_b[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func_p_b[0])( p4s[ix[0]] , p4s[ix[1]], boost ) );
  }
}

//virtual
void filler1d::operator()(cr_vec_tlv, cr_tv3, cr_tv3) {
  cout << "Calling operator() on an 1d filler with two boost vectors, doesn't make any sense, check your code!" << endl;
}



//virtual
void filler1d::operator()(cr_vec_tlv p4s, const double &w) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
  TH1* htemp = dynamic_cast<TH1*>(hist);
  if (ix.size() == 1) {
    BOOST_ASSERT_MSG(_func[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func[0])(p4s[ix[0]]), w);
  } else if (ix.size() == 2) {
    BOOST_ASSERT_MSG(_func_p[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func_p[0])(p4s[ix[0]], p4s[ix[1]]), w);
  }
}

//virtual
void filler1d::operator()(cr_vec_tlv p4s, cr_tv3 boost, const double &w) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
  TH1* htemp = dynamic_cast<TH1*>(hist);
  if (ix.size() == 1 ) {
    BOOST_ASSERT_MSG(_func_b[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func_b[0])( p4s[ix[0]], boost), w);
  } else if (ix.size() == 2 ) {
    BOOST_ASSERT_MSG(_func_p_b[0] != nullptr, filler_tag);
    htemp->Fill( (this->*_func_p_b[0])( p4s[ix[0]] , p4s[ix[1]], boost ), w );
  }
}

//virtual
void filler1d::operator()(cr_vec_tlv, cr_tv3, cr_tv3, const double &w) {
  cout << "Calling operator() on an 1d filler with two boost";
  cout << " vectors, doesn't make any sense, check your code!" << endl;
}

var1d::var1d(const int &_ix, const char *var, const axis &x, const char* frame[],
	     const char* names[][2]): filler1d(vector<int>{_ix},var) {
  set_funcs(0,var);
  set_bins(x.nbins, x.min, x.max);
  set_name( Form("%s", _name(var, frame[0], names[_ix][0]) ));
  set_title( Form("%s;%s", _title(var, frame[1], names[_ix][1]), _atitle(var, names[_ix][1]) ));
}

var1d::var1d(const int &_ix0, const int &_ix1, const char *var, const axis &x, const char* frame[],
      const char* names[][2]): filler1d(vector<int>{_ix0,_ix1},var) {
  set_pair_funcs(0,var);
  set_bins(x.nbins, x.min, x.max);
  set_name( Form("%s", _name_p(var, frame[0], names[_ix0][0], names[_ix1][0]) ));
  set_title( Form("%s;%s", _title_p(var, frame[1], names[_ix0][1], names[_ix1][1]), _atitle_p(var, names[_ix0][1], names[_ix1][1]) ));
}

filler2d::filler2d(const char* h_name, const char* h_title, const int &npart,
		   const int &nbinsx, const float &xmin, const float &xmax,
		   const int &nbinsy, const float &ymin, const float &ymax) :
  filler(npart) {
  hist = new TH2F(h_name, h_title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
}

filler2d::filler2d(const vector<int> &_ix, const vector<int> &_iy,
		   const char *varx, const char *vary):
  filler(new TH2F(),_ix,_iy),func_x(varx),func_y(vary),filler_tag(Form("Filler func tagx = %s tagy = %s", varx, vary)) {
}

void filler2d::set_bins(const int &nbinsx, const float &xmin, const float &xmax,
			const int &nbinsy, const float &ymin, const float &ymax) {
  TH2* htemp = dynamic_cast<TH2*>(hist);
  htemp->SetBins(nbinsx, xmin, xmax, nbinsy, ymin, ymax);
}

void filler2d::set_bins(const axis &x, const axis &y) {
  TH2* htemp = dynamic_cast<TH2*>(hist);
  htemp->SetBins(x.nbins, x.min, y.max, y.nbins, y.min, y.max);
}

double filler2d::get_value(const int &iaxis, const vector<int> &ii, cr_vec_tlv p4s) {
  if (ii.size() == 1) {
    return value(_func[iaxis], p4s[ii[0]]);
  } else if (ii.size() == 2) {
    return value(_func_p[iaxis], p4s[ii[0]], p4s[ii[1]]);
  } else {
    cout << "filler2d::get_value - called with too many particle indices. will fill 1.0 " << endl;
    return 1.0;
  }
}

double filler2d::get_value(const int &iaxis, const vector<int> &ii, cr_vec_tlv p4s, const TVector3 &b) {
  if (ii.size() == 1) {
    return value(_func_b[iaxis], p4s[ii[0]], b);
  } else if (ii.size() == 2) {
    return value(_func_p_b[iaxis], p4s[ii[0]], p4s[ii[1]], b);
  } else {
    cout << "filler2d::get_value - called with too many particle indices. will fill 1.0 " << endl;
    return 1.0;
  }
}

//virtual
void filler2d::operator()(cr_vec_tlv p4s) {
    BOOST_ASSERT_MSG(p4s.size()>=ix.size() , filler_tag);
    TH2* htemp = dynamic_cast<TH2*>(hist);
    htemp->Fill(get_value(0,ix,p4s),get_value(1,iy,p4s));
  }

//virtual
void filler2d::operator()(cr_vec_tlv p4s, cr_tv3 b) {
    BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
    TH2* htemp = dynamic_cast<TH2*>(hist);
    htemp->Fill(get_value(0,ix,p4s,b), get_value(1,iy,p4s,b));
  }

//virtual
void filler2d::operator()(cr_vec_tlv p4s, cr_tv3 b1, cr_tv3 b2) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
  TH2* htemp = dynamic_cast<TH2*>(hist);
  htemp->Fill(get_value(0,ix,p4s,b1), get_value(1,iy,p4s,b2));
}

//virtual
void filler2d::operator()(cr_vec_tlv p4s, const double &w) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size() , filler_tag);
  TH2* htemp = dynamic_cast<TH2*>(hist);
  htemp->Fill(get_value(0,ix,p4s),get_value(1,iy,p4s), w);
}

//virtual
void filler2d::operator()(cr_vec_tlv p4s, cr_tv3 b, const double &w) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
  TH2* htemp = dynamic_cast<TH2*>(hist);
  htemp->Fill(get_value(0,ix,p4s,b), get_value(1,iy,p4s,b), w);
}

//virtual
void filler2d::operator()(cr_vec_tlv p4s, cr_tv3 b1, cr_tv3 b2, const double &w) {
  BOOST_ASSERT_MSG(p4s.size()>=ix.size(), filler_tag);
  TH2* htemp = dynamic_cast<TH2*>(hist);
  htemp->Fill(get_value(0,ix,p4s,b1), get_value(1,iy,p4s,b2), w);
}

// Single-Single
var2d::var2d(const int &_ix, const char *varx, const axis &x, const char* framex[],
	     const int &_iy, const char *vary, const axis &y, const char* framey[],
	     const char* names[][2]): filler2d(vector<int>{_ix},vector<int>{_iy},varx,vary) {
  set_bins(x.nbins,x.min,x.max,y.nbins,y.min,y.max);
  set_funcs(0,varx); // make sure x is called before y (otherwise axes will be inversed)
  set_funcs(1,vary); // make sure x is called before y (otherwise axes will be inversed)
  if (strcmp(framex[0],framey[0])!=0) {
    set_name( Form("%s_%s", _name(varx, framex[0], names[_ix][0]), _name(vary, framey[0], names[_iy][0] ) ));
    set_title( Form("%s vs %s;%s;%s", _title(vary, framey[1], names[_iy][1]), _title(varx, framex[1], names[_ix][1] ),
		    _atitle(varx, names[_ix][1]), _atitle(vary, names[_iy][1]) ));
  } else {
    set_name( Form("%s_%s_%s", framex[0], _name(varx, names[_ix][0]), _name(vary, names[_iy][0] ) ));
    set_title( Form("%s vs %s%s;%s;%s", _title(vary, names[_iy][1]), _title(varx, names[_ix][1]),
		    framex[1], _title(varx, names[_ix][1]), _title(vary, names[_iy][1]) ));
  }
}

// Pair-Pair
var2d::var2d(const int &_ix0, const int &_ix1, const char *varx, const axis &x, const char* framex[],
      const int &_iy0, const int &_iy1, const char *vary, const axis &y, const char* framey[],
      const char* names[][2]): filler2d(vector<int>{_ix0,_ix1},vector<int>{_iy0,_iy1},varx,vary) {
  set_bins(x.nbins,x.min,x.max,y.nbins,y.min,y.max);
  set_pair_funcs(0,varx); // make sure x is called before y (otherwise axes will be inversed)
  set_pair_funcs(1,vary); // make sure x is called before y (otherwise axes will be inversed)
  if (strcmp(framex[0],framey[0])!=0) {
    set_name( Form("%s_%s", _name_p(varx, framex[0], names[_ix0][0], names[_ix1][0]), _name_p(vary, framey[0], names[_iy0][0], names[_iy1][0]) ));
    set_title( Form("%s vs %s;%s;%s", _title_p(vary, framey[1], names[_iy0][1], names[_iy1][1]), _title_p(varx, framex[1], names[_ix0][1], names[_ix1][1]),
		    _atitle_p(varx, names[_ix0][1], names[_ix1][1]), _atitle_p(vary, names[_iy0][1], names[_iy1][1]) ));
  } else {
    set_name( Form("%s_%s_%s", framex[0], _name_p(varx, names[_ix0][0], names[_ix1][0]), _name_p(vary, names[_iy0][0], names[_iy1][0]) ));
    set_title( Form("%s vs %s%s;%s;%s", _title_p(vary, names[_iy0][1], names[_iy1][1]), _title_p(varx, names[_ix0][1], names[_ix1][1]), framex[1],
		    _atitle_p(varx, names[_ix0][1], names[_ix1][1]), _atitle_p(vary, names[_iy0][1], names[_iy1][1]) ));
  }
}

// Singe-Pair
var2d::var2d(const int &_ix, const char *varx, const axis &x, const char* framex[],
	     const int &_iy0, const int &_iy1, const char *vary, const axis &y, const char* framey[],
	     const char* names[][2]): filler2d(vector<int>{_ix},vector<int>{_iy0,_iy1},varx,vary) {
  set_bins(x.nbins,x.min,x.max,y.nbins,y.min,y.max);
  set_funcs(0,varx); // make sure x is called before y (otherwise axes will be inversed)
  set_pair_funcs(1,vary); // make sure x is called before y (otherwise axes will be inversed)
  if (strcmp(framex[0],framey[0])!=0) {
    set_name( Form("%s_%s", _name(varx, framex[0], names[_ix][0]), _name_p(vary, framey[0], names[_iy0][0], names[_iy1][0]) ));
    set_title( Form("%s vs %s;%s;%s", _title_p(vary, framey[1], names[_iy0][1], names[_iy1][1]), _title(varx, framex[1], names[_ix][1]),
		    _atitle(varx, names[_ix][1]), _atitle_p(vary, names[_iy0][1], names[_iy1][1]) ));
  } else {
    set_name( Form("%s_%s_%s", framex[0], _name(varx, names[_ix][0]), _name_p(vary, names[_iy0][0], names[_iy1][0]) ));
    set_title( Form("%s vs %s%s;%s;%s", _title_p(vary, names[_iy0][1], names[_iy1][1]), _title(varx, names[_ix][1]), framex[1],
		    _atitle(varx, names[_ix][1]), _atitle_p(vary, names[_iy0][1], names[_iy1][1]) ));
  }
}

// Pair-Single
var2d::var2d(const int &_ix0, const int &_ix1, const char *varx, const axis &x, const char* framex[],
	     const int &_iy, const char *vary, const axis &y, const char* framey[], const char* names[][2])
  : filler2d(vector<int>{_ix0,_ix1},vector<int>{_iy},varx,vary) {
  set_bins(x.nbins,x.min,x.max,y.nbins,y.min,y.max);
  set_pair_funcs(0,varx); // make sure x is called before y (otherwise axes will be inversed)
  set_funcs(1,vary); // make sure x is called before y (otherwise axes will be inversed)
  if (strcmp(framex[0],framey[0])!=0) {
    set_name( Form("%s_%s", _name_p(varx, framex[0], names[_ix0][0], names[_ix1][0]), _name(vary, framey[0], names[_iy][0]) ));
    set_title( Form("%s vs %s;%s;%s", _title(vary, framey[1], names[_iy][1]), _title_p(varx, framex[1], names[_ix0][1], names[_ix1][1]),
		    _atitle_p(varx, names[_ix0][1], names[_ix1][1]), _atitle(vary, names[_iy][1]) ));
  } else {
    set_name( Form("%s_%s_%s", framex[0], _name_p(varx, names[_ix0][0], names[_ix1][0]), _name(vary, names[_iy][0]) ));
    set_title( Form("%s vs %s%s;%s;%s", _title(vary, names[_iy][1]), _title_p(varx, names[_ix0][1], names[_ix1][1]), framex[1],
		    _atitle_p(varx, names[_ix0][1], names[_ix1][1]), _atitle(vary, names[_iy][1]) ));
  }
}

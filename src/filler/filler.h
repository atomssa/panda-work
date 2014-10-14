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

using std::vector;

struct axis {
axis():nbins(200),min(0),max(10.){}
axis(int _nbins, double _min, double _max ):nbins(_nbins), min(_min), max(_max) {}
  int nbins;
  double min;
  double max;
};

//____________
class filler {

 public:

  filler(TNamed *_hist, const vector<int> &_ix);
  filler (TNamed *_hist, const vector<int> &_ix, const vector<int> &_iy);

  // LEGACY
  filler (TNamed *_hist, const int &npart);
  filler(const int &npart);

  ~filler();

  virtual void operator()(const vector<TLorentzVector>&)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&)=0;
  // if x axis and y axis need to be boosted by different vecotrs -
  // this function should do nothing and communicate annoyingly if called on a 1d filler
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&, const TVector3&)=0;

  // with weight
  virtual void operator()(const vector<TLorentzVector>&, const double &)=0;
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&, const double &)=0;
  // if x axis and y axis need to be boosted by different vecotrs -
  // this function should do nothing and communicate annoyingly if called on a 1d filler
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&, const TVector3&, const double &)=0;

  void Write();
  void Write(const char* prefix);

 protected:

  typedef double(filler::*func)(const TLorentzVector&) const;
  typedef double(filler::*func_boost)(const TLorentzVector&, const TVector3&) const;
  typedef double(filler::*pair_func)(const TLorentzVector&, const TLorentzVector&) const;
  typedef double(filler::*pair_func_boost)(const TLorentzVector&, const TLorentzVector&, const TVector3&) const;

  std::map<const char*, const char*> vart; // variable tex-form for titles
  std::map<const char*, const char*> varst; // variable short tex-form for axes
  std::map<const char*, const char*> varu; // variable units
  std::map<const char*, axis> varb; // default axis for histgram

  std::map<const char*, func> func_dict;
  std::map<const char*, func_boost> func_boost_dict;
  std::map<const char*, pair_func> pair_func_dict;
  std::map<const char*, pair_func_boost> pair_func_boost_dict;

  TNamed *hist;
  std::vector<int> ix;
  std::vector<int> iy;

  std::vector<func> _func;
  std::vector<func_boost> _func_b;
  std::vector<pair_func> _func_p;
  std::vector<pair_func_boost> _func_p_b;

  void set_name(const char*name);
  void set_title(const char*title);

  const char* warning_msg(const char* var, const char* dict);

  void init_var(const char* var, const char* title, const char *stitle, const char*unit, const axis& bins);
  void init_vars( );
  void init_pair_func_dict();
  void init_pair_func_boost_dict();
  void init_func_dict();
  void init_func_boost_dict();
  void init_all();

  TLorentzVector boost_transf(const TLorentzVector &vect_in, const TVector3 &boost) const;

  void set_funcs(const int &iaxis, const char* var);
  void set_pair_funcs(const int &iaxis, const char* var);

  double value(func __func, const TLorentzVector &v);
  double value(func_boost __func, const TLorentzVector &v, const TVector3 &b);
  double value(pair_func __func, const TLorentzVector &v1, const TLorentzVector &v2);
  double value(pair_func_boost __func, const TLorentzVector &v1, const TLorentzVector &v2, const TVector3 &b);

  const char* _name(const char* var,const char* fname, const char* pname);
  const char* _name(const char* var, const char* pname);
  const char* _name_p(const char* var,const char* fname, const char* pname1, const char* pname2);
  const char* _name_p(const char* var, const char* pname1, const char* pname2);
  const char* _title(const char *var, const char *ftitle, const char *ptitle);
  const char* _title(const char *var, const char *ptitle);
  const char* _title_p(const char *var, const char *ftitle, const char *ptitle1, const char *ptitle2);
  const char* _title_p(const char *var, const char *ptitle1, const char *ptitle2);
  const char* _atitle(const char *var, const char *ptitle);
  const char* _atitle_p(const char *var, const char *ptitle1, const char *ptitle2);

  double mass(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double mass_sq(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double mand_s(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double mand_u(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double mand_t(const TLorentzVector &v1, const TLorentzVector &v2) const;

  double mom(const TLorentzVector &v) const;
  double mom_b(const TLorentzVector &v, const TVector3 &b) const;
  double mom_p(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double mom_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const;
  double mom_p_fwd(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double mom_p_bwd(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double mom_p_fwd_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3 &b) const;
  double mom_p_bwd_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3 &b) const;

  double pt(const TLorentzVector &v) const;
  double pt_b(const TLorentzVector &v, const TVector3 &boost) const;
  double pt_p(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double pt_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const;

  double ene(const TLorentzVector &v) const;
  double ene_b(const TLorentzVector &v, const TVector3 &boost) const;
  double ene_p(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double ene_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const;

  double the(const TLorentzVector &v) const;
  double the_b(const TLorentzVector &v, const TVector3 &boost) const;
  double the_p(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double the_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const;
  double the_p_fwd(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double the_p_bwd(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double the_p_fwd_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3 &b) const;
  double the_p_bwd_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3 &b) const;

  double phi(const TLorentzVector &v) const;
  double phi_p(const TLorentzVector &v1, const TLorentzVector &v2) const;

  double cost(const TLorentzVector &v) const;
  double cost_b(const TLorentzVector &v, const TVector3 &boost) const;
  double cost_p(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double cost_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const;

  double oa_p(const TLorentzVector &v1, const TLorentzVector &v2) const;
  double oa_p_b(const TLorentzVector &v1, const TLorentzVector &v2, const TVector3&b) const;

};

//_____________________________
class filler1d: public filler {

 protected:
  const char* func_tag;
  const char* filler_tag;

  filler1d(const char* h_name, const char* h_title, const int &npart, const int &nbins, const float &min, const float &max );
  filler1d(const int &npart);
  filler1d(const vector<int> &_ix, const char *var);
  ~filler1d();

  void set_bins(const int &nbins, const float &min, const float &max);

 public:
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector> &p4s);
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3& boost);
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&, const TVector3&);
  // Virtual overlaod function call operators to pass the 4momenta to be filled
  virtual void operator()(const vector<TLorentzVector> &p4s, const double &w);
  // if needs to be boosted by some
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3& boost, const double &w);
  virtual void operator()(const vector<TLorentzVector>&, const TVector3&, const TVector3&, const double &w);

  TH1* getHist() {return dynamic_cast<TH1*>(hist); }

};

//____________________________________
class var1d: public filler1d {
 public:
 var1d(const int &_ix, const char *var, const axis &x, const char* frame[], const char* names[][2]);
 var1d(const int &_ix0, const int &_ix1, const char *var, const axis &x, const char* frame[], const char* names[][2]);
};

//_____________________________
class filler2d: public filler {

 protected:
  const char *func_x, *func_y;
  const char *filler_tag;

 filler2d(const char* h_name, const char* h_title, const int &npart,
	  const int &nbinsx, const float &xmin, const float &xmax,
	  const int &nbinsy, const float &ymin, const float &ymax);

 filler2d(const vector<int> &_ix, const vector<int> &_iy,
	  const char *varx, const char *vary);

  void set_bins(const int &nbinsx, const float &xmin, const float &xmax,
		const int &nbinsy, const float &ymin, const float &ymax);

  void set_bins(const axis &x, const axis &y);

  ~filler2d() { delete hist; }

 private:

  double get_value(const int &iaxis, const vector<int> &ii, const vector<TLorentzVector> &p4s);
  double get_value(const int &iaxis, const vector<int> &ii, const vector<TLorentzVector> &p4s, const TVector3 &b);

 public:
  virtual void operator()(const vector<TLorentzVector> &p4s);
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3& b);
  virtual void operator()(const vector<TLorentzVector>&p4s, const TVector3&b1, const TVector3&b2);

  virtual void operator()(const vector<TLorentzVector> &p4s, const double &w);
  virtual void operator()(const vector<TLorentzVector> &p4s, const TVector3& b, const double &w);
  virtual void operator()(const vector<TLorentzVector>&p4s, const TVector3&b1, const TVector3&b2, const double &w);

  TH2* getHist() {return dynamic_cast<TH2*>(hist); }

};

//____________________________________
class var2d: public filler2d {
 public:
  // Example consructions
  //var2d(ep,em,"mass",mass_bins,lab_frame,jpsi,"pt",pt_bins,lab_frame,parts,part_titles)
  //var2d(jpsi,"energy",ene_bins,lab_frame,ep,em,"oa",oa_bins,lab_frame,parts,part_titles)
  //var2d(gamma1,"mom",mom_bins,jpsi_frame,)
  // Constructor chaining is allowed only in c++11

  // Single-Single
 var2d(const int &_ix, const char *varx, const axis &x, const char* framex[],
       const int &_iy, const char *vary, const axis &y, const char* framey[],
       const char* names[][2]);

 // Pair-Pair
 var2d(const int &_ix0, const int &_ix1, const char *varx, const axis &x, const char* framex[],
       const int &_iy0, const int &_iy1, const char *vary, const axis &y, const char* framey[],
       const char* names[][2]);

 // Singe-Pair
 var2d(const int &_ix, const char *varx, const axis &x, const char* framex[],
       const int &_iy0, const int &_iy1, const char *vary, const axis &y, const char* framey[],
       const char* names[][2]);

 // Pair-Single
 var2d(const int &_ix0, const int &_ix1, const char *varx, const axis &x, const char* framex[],
       const int &_iy, const char *vary, const axis &y, const char* framey[], const char* names[][2]);

};

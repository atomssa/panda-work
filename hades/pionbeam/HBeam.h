#ifndef _HBEAM_H_
#define _HBEAM_H_

/////////////////////////////////////////////////////////////////////
//
// Beam line transport
//
// Author:  T. Hennino (core), J. Biernat , J.Markert
//          Pluto implementation by IF
// Written: 14.11.2013
//
/////////////////////////////////////////////////////////////////////

//_HADES_CLASS_DESCRIPTION
///////////////////////////////////////////////////////////////////////////////
//
// HBeam
//
//    This class simulates a beam line. Beam particles (HBeamParticle) will be produced using
//    a certain beam profile and momentum smearing. The particle will be tranported through
//    the beam line using matrices for the single beam elements (HBeamElement). The particles
//    cand be registered in detectors along the beam line.
//
//    HBeam pionbeam;
//    pionbeam.setBeam           (HPhysicsConstants::pid("pi-"),1.3,60,60,0.0,0.0); // id, total mom [GeV], beamtube x and y size, xoff,yoff
//    pionbeam.SetBeamResolution (0.01,0.05,0.06);                                  // dpx [rad],dpy [rad],dpz [frac]  [+-] white random
//    pionbeam.setBeamProfile    (.05,0.0);                                         // sigma [mm], flatradius [mm]
//                                                                                  // sigma>0 will select the beam profile (gausian + flat top)
//                                                                                  // sigma==0, flatradius>0 will give an extended beam spot without gaussian borders
//                                                                                  // sigma==0,flatradius==0 spot like beam
//    if(!pionbeam.initBeamLine  ("par_files/pibeam_set6_mod.data",32)) return;               // transform input file and target element number
//    pionbeam.addDetector("det1", -17092.6,2,50.,50.);                             // [mm] relative to HADES 0,0,0    cutype (0,1,2), xcut[mm],ycut[mm]
//    pionbeam.addDetector("det2",  -5400.0,2,50.,50.);                             // [mm] relative to HADES 0,0,0    cutype (0,1,2), xcut[mm],ycut[mm]
//    pionbeam.addDetector("plane", -1300.0,1,60.,60.);                             // [mm] relative to HADES 0,0,0    cutype (0,1,2), xcut[mm],ycut[mm]
//
//    pionbeam.printBeamLine(kTRUE);          // kTRUE : print transform matrices in addition to name and distance
//    pionbeam.printDetectors();
//    pionbeam.printBeamProfile();
//
//
//    //-----------------------------------------------------------
//    // create particles
//    vector<HBeamParticle>& vhistory = pionbeam.newParticle();  // beam particle at : beam,det1,det2,plane.target
//    // check if particle was in acceptance (using the particle at the end of transport)
//    if(vhistory[vhistory.size()-1].fAccepted) {   ......   }
//
//    // get access to all beam elements and detectors (positions,accptance,statistics ...)
//    vector<HBeamElement>& elements  = pionbeam.getElements();
//    vector<HBeamElement>& detectors = pionbeam.getDetectors();
//
///////////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TLorentzVector.h"

#include "TF1.h"

#include "TCanvas.h" // temporary for debug

#include "Math/Functor.h"
#include "Math/RootFinder.h"

#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

#include "TMinuit.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

double Tx[3][8] = {{0.}};
double Ty[3][7] = {{0.}};
double x_mes[3]={0.}, x_mes_e[3]={0.};
double y_mes[3]={0.}, y_mes_e[3]={0.};

//______________________________________________________________________________
Double_t x_exp(int i, Double_t *par)
{
  // Tx[0]->T11, Tx[1]->T12, Tx[2]->T14, Tx[3]->T16, Tx[4]->T116, Tx[5]->T126, Tx[6]->T146, Tx[7]->T166
  // par[0]->x, par[1]->theta, par[2]->y, par[3]->phi, par[4]->del
  return Tx[i][0]*par[0] + Tx[i][1]*par[1] + Tx[i][2]*par[3] + Tx[i][3]*par[4] + Tx[i][4]*par[0]*par[4] + Tx[i][5]*par[1]*par[4] + Tx[i][6]*par[3]*par[4] + Tx[i][7]*par[4]*par[4];
}

//______________________________________________________________________________
Double_t y_exp(int i, Double_t *par)
{
  // Ty[0]->T32, Ty[1]->T33, Ty[2]->T34, Ty[3]->T36, Ty[4]->T336, Ty[5]->T346, Ty[6]->T366
  // par[0]->x, par[1]->theta, par[2]->y, par[3]->phi, par[4]->del
  return Ty[i][0]*par[1] + Ty[i][1]*par[2] + Ty[i][2]*par[3] + Ty[i][3]*par[4] + Ty[i][4]*par[2]*par[4] + Ty[i][5]*par[3]*par[4] + Ty[i][6]*par[4]*par[4];
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   const Int_t ndet = 3;
   //calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (int idet=0; idet<ndet; ++idet) {
     delta = ( x_mes[idet] - x_exp(idet, par) ) / x_mes_e[idet];
     chisq += delta*delta;
     delta = ( y_mes[idet] - y_exp(idet, par) ) / y_mes_e[idet];
     chisq += delta*delta;
   }
   f = chisq;
   //cout << "chisq = " << chisq << endl;
}

// struct to hold reconstruction parameters
struct RecoPars {
  // solve_delta and solve_theta
  double T12d1, T12d2, T16d1, T16d2, T126d1, T126d2, T166d1, T166d2, T14d1, T14d2, T146d1, T146d2;
  // solve_phi_and_y
  double T33d1, T33d2, T34d1, T34d2, T32d1, T32d2, T36d1, T36d2, T336d1, T336d2, T346d1, T346d2, T366d1, T366d2;
};

class HBeamElement {

public:

    TString  fName;          // name of the lement : default  = Element[i]
    Double_t fDistance;      // position along beam line
    Int_t    fId;            // element number
    Double_t Tij [5][5];     // first order transform
    Double_t Tijk[5][5][5];  // second order transform
    Double_t fout[5];        // particle state at element [x,tx,y,ty,dp]
    Double_t fxCut;          // cut ond x for apperture
    Double_t fyCut;          // cut ond y for apperture
    Int_t    fCutType;       // 0: no cut ,1 : radial cut (default), 2: square cut
    // square cut uses x,y , radial cut only x!
    Long64_t fCtAll ;        // count all beam particle transport
    Long64_t fCtFail;        // count if particle is out of apperture
    Bool_t   fAccepted;      // remembers the last result of check()
    Double_t fXoX0;          // MS angular smearing width (Sigma of gaussian)
    Double_t fPitch;         // Detector pitch size (mm)
    Bool_t   fDigi;          // If the hit position on this detector needs to be discretized

    HBeamElement(TString n="not_set",Double_t dist = 0, Int_t id=-1)
    {
	fName     = n;
	fDistance = dist;
	fId       = id;
	fCutType  = 1;
	fxCut     = 60;
	fyCut     = 60;
	fCtAll    = 0;
	fCtFail   = 0;
	fAccepted = kTRUE;
	clear();
    }

    void setElement(TString name,Int_t cuttype,Double_t xcut,Double_t ycut){
	if(name.CompareTo("")!=0) fName = name;
	fCutType = cuttype;
	fxCut = xcut;
        fyCut = ycut;
    }

    bool apply_ms() {
      return ( fXoX0 > 1e-9 );
    }

    void clear()
    {
	for(Int_t i = 0; i < 5; i ++) {
	    fout[i] = 0;
	    for(Int_t j = 0; j < 5; j ++){
		Tij[i][j] =0;
		for(Int_t k = 0; k < 5; k ++){
		    Tijk[i][j][k] = 0;
		}
	    }
	}
    }

    void printTij()
    {
	for(Int_t i = 0; i < 5; i ++) {     // row
	    for(Int_t j = 0; j < 5; j ++){  // column
		cout<<setw(12)<<Tij[i][j]<<" "<<flush;

	    }
	    cout<<endl;
	}
    }

    void printTijk()
    {
	for(Int_t i = 0; i < 5; i ++) { // block
	    for(Int_t k = 0; k < 5; k ++){ // line
		for(Int_t j = 0; j <= k; j ++){    // values
		    cout<<setw(12)<<Tijk[i][j][k]<<" "<<flush;
		}
		cout<<endl;
	    }
	    cout<<endl;
	}
    }


    void print(Bool_t printAll = kFALSE)
    {
	// printAll = kTRUE will print transforms in addition
	cout<<"---------------------------------"<<endl;
	cout<<"ID  : "<<setw(3)<<fId<<" , "<<fName.Data()<<endl;
	cout<<"distance   : "<<setw(12)<<fDistance<<endl;
	cout<<"CutType    : "<<setw(3)<<fCutType<<" , x : "<<setw(12)<<fxCut<<", y : "<<setw(12)<<fyCut<<endl;
	cout<<"Acceptance : "<<setw(8)<<fCtAll<<" , Failed : "<<setw(8)<<fCtFail<<" ,  "<<setw(8)<<(fCtAll > 0?  (fCtFail/(Double_t)fCtAll)*100: 100)<<" %"<<endl;
        if(printAll){
	    cout<<"first  order transform:"<<endl;
	    printTij();
	    cout<<"second order transform:"<<endl;
	    printTijk();
	}
    }

    Bool_t isInAcceptance()
    {
	// checks if the particle is inside acceptance
	// for this element. to be called after beam transport (for elements)
	// and after createDetectorHits() for detectors
	if(fCutType == 0 || fCutType > 2) return kTRUE;
	if(fCutType == 1 ){
	    Double_t radius = TMath::Sqrt(fout[0] * fout[0] + fout[2] * fout[2]) ;
	    if(radius <  fxCut )   return kTRUE;
	    else                   return kFALSE;
	}
	if(fCutType == 2 ){
	    if(TMath::Abs(fout[0]) < fxCut && TMath::Abs(fout[2]) < fyCut ) return kTRUE;
	    else                                                            return kFALSE;
	}
	return kTRUE;
    }

    void check()
    {
	// checks the acceptance and fill the stats counters
	// and the fAccepted flag
	if(!isInAcceptance()) {
	    fCtFail++;
	    fAccepted = kFALSE;
	} else  fAccepted = kTRUE;

	fCtAll++;
    }
};


class HBeamParticle {

public:
    TVector3 fP;                 // momentum vector
    TVector3 fPos;               // position vector
    TVector3 fOffset;            // beam could be displaced : profile function makes use of it
    TVector3 fMomRes;            // momentum resolution

    Double_t fBeamMom;           // nominal beam momentum
    Double_t fBeamMomSmeared;    // smeared beam momentum
    Int_t    fPid;               // type of beam particle
    Int_t    fDetnum;            // detector number where the transported particle was registered
    TString  fName;              // name of detecor where particle was registered
    Double_t fbeamtube_size_x;   // x general acceptance of the beam tube
    Double_t fbeamtube_size_y;   // y general acceptance of the beam tube
    Bool_t   fAccepted;          // whether the particle was reaching the detector without violation acceptance
    Int_t    fStatus;            // Status of TDR solution
    HBeamParticle(){
	clear();
    }

    void clear (){
	fP     .SetXYZ(0.,0.,0.);
	fPos   .SetXYZ(0.,0.,0.);
	fOffset.SetXYZ(0.,0.,0.);
	fMomRes.SetXYZ(0.,0.,0.);
	fBeamMom         = -1.;
	fBeamMomSmeared  = -1.;
	fPid             = -10;
	fDetnum          = -1;
	fbeamtube_size_x = -1.;
	fbeamtube_size_y = -1.;
	fAccepted        = kTRUE;
	fStatus          = -10000;
    }

    // access methods to return "state vector" elements concisely
    Double_t x() const { return fPos.X(); }
    Double_t th() const { return TMath::ATan2(fP.X(), fP.Z()); }
    Double_t y() const { return fPos.Y(); }
    Double_t ph() const { return TMath::ATan2(fP.Y(), fP.Z()); }
    Double_t dp() const { return (fBeamMomSmeared-fBeamMom)/fBeamMom; }

    // method to set fPos and fP from "state vector",
    // additional z position and beam momentum required because state vector doesn't hold this info
    // Units of state vector have to be pluto units
    // assumes state array has 5 elements
    void set_state(double state[], double beam_mom, double z) {
      fPos.SetXYZ(state[0], state[2], z );
      const double p  = beam_mom * (1. + state[4]);
      const double pz = p / sqrt(1.+ tan(state[1])*tan(state[1]) + tan(state[3])*tan(state[3]));
      const double px = pz * tan(state[1]);
      const double py = pz * tan(state[3]);
      fP.SetXYZ(px,py,pz);
    }

    void print(Bool_t printAll = kFALSE){
	// printAll=kTRUE print in addition momentum vector
	cout<<"---------------------------------"<<endl;
	cout<<setw(12)<<fName.Data()<<" acc : "<<fAccepted<<" , ID "<<fPid<<" , Beam Mom : "<<setw(12)<<fBeamMom<<endl;
	cout<<"Beam tube [x,y]       : "<<setw(12)<<fbeamtube_size_x <<" , "<<setw(12)<<fbeamtube_size_y <<endl;
	cout<<"P smear   [x,y,z]     : "<<setw(12)<<fMomRes.X()<<" , "<<setw(12)<<fMomRes .Y()<<" , "<<setw(12)<<fMomRes.Z()<<endl;
	cout<<"Position  [x,y,z]     : "<<setw(12)<<fPos   .X()<<" , "<<setw(12)<<fPos    .Y()<<" , "<<setw(12)<<fPos   .Z()<<endl;
	cout<<"Offset    [x,y,z]     : "<<setw(12)<<fOffset.X()<<" , "<<setw(12)<<fOffset .Y()<<" , "<<setw(12)<<fOffset.Z()<<endl;

	if(printAll){
	    cout<<"P         [x,y,z,tot] : "<<setw(12)<<fP     .X()<<" , "<<setw(12)<<fP      .Y()<<" , "<<setw(12)<<fP     .Z()<<" , "<<setw(12)<<fP  .Mag()<<endl;
	    cout<<"Detector ID           : "<<setw(12)<<fDetnum    <<endl;
	}


    }

};




class HBeam {



private:

    static Double_t profile(Double_t *x, Double_t *par)
    {
	Double_t xx    = x[0];
	Double_t c     = 1;
	Double_t mean  = 0;
	Double_t sigma = par[0];
	Double_t flat  = par[1];
	Double_t result = 0;

	if(xx < mean - flat) {         // left gaus
	    result = c * TMath::Gaus(xx,mean - flat,sigma);
	} else if(xx > mean + flat) {  // right gaus
	    result = c * TMath::Gaus(xx,mean + flat,sigma);
	} else {                       // flat top
	    result = c;
	}

	return result;
    }


    HBeamParticle fBeam;                 // beam particle
    TF1*          fProfile;              // beam profile

    vector <HBeamElement> fdetectors;    // vector of all detectors
    vector <HBeamElement>  felements;    // vector of all beam elements
    vector <HBeamParticle>  fhistory;    // vector of beam particle at all detectors
    vector <HBeamParticle>  fms_history;  // vector of beam particle state after each MS generation step
    vector <HBeamParticle>  fsolution;   // vector to hold solution states
    Int_t                fnum_target;    // index of elemnt target
    Double_t            ftoPluto  [5];   // units -> transport ->Pluto
    Double_t            ffromPluto[5];   // units -> Pluto     ->transport
    Bool_t setTargetElement(UInt_t n);

    TMinuit *gMinuit; // for fitting
    Double_t arglist[10];
    Int_t ierflg;

    vector<Int_t> det_idx;
    RecoPars par;

public:

    HBeam();
    ~HBeam();
    Int_t                    addDetector (TString name, Double_t thikness, Double_t rad_len, bool digi, Double_t pitch, Double_t distance,Int_t cutType=2,Double_t xcut=50,Double_t ycut=50);
    vector<HBeamElement>&    getDetectors()     { return fdetectors ;  }
    vector<HBeamElement>&    getElements()      { return felements ;   }
    vector<HBeamParticle>&   getHistory()       { return fhistory ;    }

    void                     setBeam           (Int_t id,Double_t p, Double_t _beamtube_size_x,Double_t _beamtube_size_y,Double_t xoff = 0.0,Double_t yoff = 0.0);
    void                     setBeamResolution (Double_t dpx = 0.001,Double_t dpy = 0.005,Double_t dpz = 0.06);
    void                     setBeamProfile    (Double_t sigma_beam = 0.05,Double_t flat_radius = 0);
    Bool_t                   initBeamLine      (TString filename,Int_t targetElementNum ,Bool_t debug = kFALSE);
    vector<HBeamParticle>&   newParticle       ();
    void                     printBeamLine     (Bool_t printAll=kFALSE);
    void                     printBeamProfile  ();
    void                     printDetectors    ();

    void                     createBeam        (HBeamParticle& part);
    Bool_t                   transport         (HBeamParticle& part);
    Bool_t                   transport         (HBeamParticle& part, int from);
    Bool_t                   transport_with_ms (HBeamParticle& part);
    Bool_t                   createDetectorHits(HBeamParticle& part);

    void                     setRecoParams(vector<Int_t> idx);

    Double_t                 digitize_pos(const Double_t &pos, const Double_t &pitch);
    Double_t                 solve_delta(const Double_t &xd1, const Double_t &xd2, const Double_t &phi, int &status);
    Double_t                 solve_theta(const Double_t &xd1, const Double_t &xd2, const Double_t &phi, const Double_t &delta);
    void                     solve_phi_and_y(const Double_t &yd1, const Double_t &yd2, const Double_t &th, const Double_t &del, Double_t &ph, Double_t &y);
    void                     transport_to(const HBeamElement &det, const vector<Double_t> &state_in, vector<Double_t> &state_out);
    void                     solve_state();

    void                     store_solution(const char *tag, vector<Double_t> &state_targ, HBeamParticle& part);
    int                     solution_tdr(Double_t *x, Double_t *y, vector<Double_t> &state);
    int                     solution_minuit(vector<Double_t> &state);

    vector<HBeamParticle>&   get_ms_history() { return fms_history; }
    vector<HBeamParticle>&   get_solution() { return fsolution; }

    TF1 *gen_func;

};

HBeam::HBeam()
{
  fnum_target = -1;
  fProfile = new TF1("fprofile",profile,-20,20,2);
  fProfile->SetParNames("sigma","flatradius");
  setBeamProfile(0.05,0);
  setBeamResolution(0.01,0.05,0.06); //
  ftoPluto[0] = 10.;      // x  cm   -> mm
  ftoPluto[1] =  0.001;   // tx mrad -> rad
  ftoPluto[2] = 10.;      // y  cm   -> mm
  ftoPluto[3] =  0.001;   // ty mrad -> rad
  ftoPluto[4] =  0.01;    // dp %    -> frac

  for(Int_t i = 0; i < 5; i ++) ffromPluto[i] = 1./ftoPluto[i];

  // Initialize gMinuit for track fitting
  gMinuit = new TMinuit(5);
  gMinuit->SetFCN(fcn);
  gMinuit->SetPrintLevel(-1);

  gen_func = new TF1("pol2","[0]+x*x*[1]",-1.1999999,1.19999999);
  gen_func->SetParameter(0,1);
  gen_func->SetParameter(1,100);

}

HBeam::~HBeam() {
    if(fProfile) delete fProfile;
}

Int_t HBeam::addDetector(TString name, Double_t thikness, Double_t rad_len, bool digi, Double_t pitch, Double_t distance, Int_t cutType, Double_t xcut, Double_t ycut)
{
    // adds a detector "name" at place "distance" [in mm]
    // (distance must be negative if placed before the target)
    // cutType (0,1,2)(no,radial,square) xcut/ycut [mm]

    HBeamElement det(name,distance,fdetectors.size()+1);
    det.fxCut    = TMath::Abs(xcut);
    det.fyCut    = TMath::Abs(ycut);
    det.fCutType = cutType;

    det.fXoX0 = thikness / rad_len;

    det.fPitch = pitch;
    det.fDigi = digi;

    cout << "detector: " << name << " dist= " << distance << endl;
    // Intrapolate detector transport matrix elements from nearest neighboring entries
    for (unsigned int iel=0; iel<felements.size(); ++iel) {
      if ( felements[iel].fDistance <= distance && distance <= felements[iel+1].fDistance ) {
	Double_t dist_num = (distance - felements[iel].fDistance);
	Double_t dist_den = (felements[iel+1].fDistance - felements[iel].fDistance);
	if (dist_den != 0 ) {
	  Double_t frac = dist_den!=0? dist_num/dist_den: 1.0;
	  for (int j=0; j<5; ++j) {
	    det.Tij[0][j] = felements[iel].Tij[0][j] + frac* (felements[iel+1].Tij[0][j] - felements[iel].Tij[0][j])  ;
	    det.Tij[1][j] = felements[iel].Tij[1][j];
	    det.Tij[2][j] = felements[iel].Tij[2][j] + frac* (felements[iel+1].Tij[2][j] - felements[iel].Tij[2][j])  ;
	    det.Tij[3][j] = felements[iel].Tij[3][j];
	    det.Tij[4][j] = felements[iel].Tij[4][j];
	  }
	  for (int j=0; j<5; ++j) {
	    for (int k=0; k<5; ++k) {
	      det.Tijk[0][j][k] = felements[iel].Tijk[0][j][k] + frac * (felements[iel+1].Tijk[0][j][k] - felements[iel].Tijk[0][j][k] );
	      det.Tijk[1][j][k] = felements[iel].Tijk[1][j][k];
	      det.Tijk[2][j][k] = felements[iel].Tijk[2][j][k] + frac * (felements[iel+1].Tijk[2][j][k] - felements[iel].Tijk[2][j][k] );
	      det.Tijk[3][j][k] = felements[iel].Tijk[3][j][k];
	      det.Tijk[4][j][k] = felements[iel].Tijk[4][j][k];
	    }
	  }
	} else {
	  //special case for hades target added as "detector". Interpolating between the last two elements doesn't work
	  // Therefore, just copy the matrix elements from either of two elements
	  for (int ii=0; ii<5; ++ii)
	    for (int jj=0; jj<5; ++jj)
	      det.Tij[ii][jj] = felements[iel].Tij[ii][jj];

	  for (int ii=0; ii<5; ++ii)
	    for (int jj=0; jj<5; ++jj)
	      for (int kk=0; kk<5; ++kk)
		det.Tijk[ii][jj][kk] = felements[iel].Tijk[ii][jj][kk];
	}
      }
    }

    fdetectors.push_back(det);
    return fdetectors.size();

}

void HBeam::setRecoParams(vector<Int_t> idx) {

  for (unsigned int i=0; i<idx.size(); ++i) det_idx.push_back(idx[i]);

  // solve_delta
  par.T12d1 = fdetectors[det_idx[0]].Tij[0][1];
  par.T12d2 = fdetectors[det_idx[1]].Tij[0][1];
  par.T16d1 = fdetectors[det_idx[0]].Tij[0][4];
  par.T16d2 = fdetectors[det_idx[1]].Tij[0][4];
  par.T126d1 = fdetectors[det_idx[0]].Tijk[0][1][4];
  par.T126d2 = fdetectors[det_idx[1]].Tijk[0][1][4];
  par.T166d1 = fdetectors[det_idx[0]].Tijk[0][4][4];
  par.T166d2 = fdetectors[det_idx[1]].Tijk[0][4][4];
  par.T14d1 = fdetectors[det_idx[0]].Tij[0][3];
  par.T14d2 = fdetectors[det_idx[1]].Tij[0][3];
  par.T146d1 = fdetectors[det_idx[0]].Tijk[0][3][4];
  par.T146d2 = fdetectors[det_idx[1]].Tijk[0][3][4];

  // solve_phi_and_y
  par.T33d1 = fdetectors[det_idx[0]].Tij[2][2];
  par.T33d2 = fdetectors[det_idx[1]].Tij[2][2];
  par.T34d1 = fdetectors[det_idx[0]].Tij[2][3];
  par.T34d2 = fdetectors[det_idx[1]].Tij[2][3];
  par.T32d1 = fdetectors[det_idx[0]].Tij[2][1];
  par.T32d2 = fdetectors[det_idx[1]].Tij[2][1];
  par.T36d1 = fdetectors[det_idx[0]].Tij[2][4];
  par.T36d2 = fdetectors[det_idx[1]].Tij[2][4];
  par.T336d1 = fdetectors[det_idx[0]].Tijk[2][2][4];
  par.T336d2 = fdetectors[det_idx[1]].Tijk[2][2][4];
  par.T346d1 = fdetectors[det_idx[0]].Tijk[2][3][4];
  par.T346d2 = fdetectors[det_idx[1]].Tijk[2][3][4];
  par.T366d1 = fdetectors[det_idx[0]].Tijk[2][4][4];
  par.T366d2 = fdetectors[det_idx[1]].Tijk[2][4][4];

  // Static parameters for Minuit minimization function
  // Tx[0]->T11, Tx[1]->T12, Tx[2]->T14, Tx[3]->T16, Tx[4]->T116, Tx[5]->T126, Tx[6]->T146, Tx[7]->T166
  // Ty[0]->T32, Ty[1]->T33, Ty[2]->T34, Ty[3]->T36, Ty[4]->T336, Ty[5]->T346, Ty[6]->T366
  // int det_idx[3] = {0 /*det1*/, 1 /*det2*/, 3 /*diam*/};
  for (int idet=0; idet<3; ++idet) {
    cout << "det_idx[" << idet << "] = " << det_idx[idet] << " name= " << fdetectors[det_idx[idet]].fName << endl;
    Tx[idet][0] /*T11*/ = fdetectors[det_idx[idet]].Tij[0][0];
    Tx[idet][1] /*T12*/ = fdetectors[det_idx[idet]].Tij[0][1];
    Tx[idet][2] /*T14*/ = fdetectors[det_idx[idet]].Tij[0][3];
    Tx[idet][3] /*T16*/ = fdetectors[det_idx[idet]].Tij[0][4];
    Tx[idet][4] /*T166*/ = fdetectors[det_idx[idet]].Tijk[0][4][4];
    Tx[idet][5] /*T126*/ = fdetectors[det_idx[idet]].Tijk[0][1][4];
    Tx[idet][6] /*T146*/ = fdetectors[det_idx[idet]].Tijk[0][3][4];
    Tx[idet][7] /*T166*/ = fdetectors[det_idx[idet]].Tijk[0][4][4];

    Ty[idet][0] /*T32*/ = fdetectors[det_idx[idet]].Tij[2][1];
    Ty[idet][1] /*T33*/ = fdetectors[det_idx[idet]].Tij[2][2];
    Ty[idet][2] /*T34*/ = fdetectors[det_idx[idet]].Tij[2][3];
    Ty[idet][3] /*T36*/ = fdetectors[det_idx[idet]].Tij[2][4];
    Ty[idet][4] /*T336*/ = fdetectors[det_idx[idet]].Tijk[2][2][4];
    Ty[idet][5] /*T346*/ = fdetectors[det_idx[idet]].Tijk[2][3][4];
    Ty[idet][6] /*T366*/ = fdetectors[det_idx[idet]].Tijk[2][4][4];
  }

}

Bool_t HBeam::setTargetElement(UInt_t n){
    if (n < 1 || n >= felements.size()){
	cout<<"HBeam::setTargetElement() : specifified target number outside range!"<<endl;
	return kFALSE;
    }
    fnum_target = n;
    Double_t shift = felements[fnum_target].fDistance;
    for (UInt_t i = 0; i < felements.size(); i ++){
	felements[i].fDistance -= shift;
    }
    return kTRUE;
}


void HBeam::setBeam(Int_t id,
			Double_t p,
			Double_t _beamtube_size_x,Double_t _beamtube_size_y,
			Double_t xoff,Double_t yoff)
{
    fBeam.fBeamMom = TMath::Abs(p);
    fBeam.fPid     = id;
    fBeam.fOffset.SetXYZ(xoff,yoff,0.);
    fBeam.fbeamtube_size_x = TMath::Abs(_beamtube_size_x);
    fBeam.fbeamtube_size_y = TMath::Abs(_beamtube_size_y);
}

void HBeam::setBeamResolution(Double_t dpx,Double_t dpy,Double_t dpz)
{
    fBeam.fMomRes.SetXYZ(TMath::Abs(dpx),TMath::Abs(dpy),TMath::Abs(dpz));
}

void HBeam::setBeamProfile(Double_t sigma_beam,Double_t flat_radius)
{
    fProfile->SetParameter(0,TMath::Abs(sigma_beam));
    fProfile->SetParameter(1,TMath::Abs(flat_radius));
}


void  HBeam::createBeam(HBeamParticle& part)
{

    //----------------------------------------
    // beam profile
    // x-y symmetric profile
    Double_t x = 0;
    if      ( fProfile->GetParameter(0) != 0)                                  x  = fProfile->GetRandom();
    else if ( fProfile->GetParameter(0) == 0  && fProfile->GetParameter(1)!=0) x  = 2 * (gRandom->Rndm() - 0.5) * fProfile->GetParameter(3);
    else                                                                       x  = 0;
    Double_t y = 0;
    if      ( fProfile->GetParameter(0) != 0)                                  y  = fProfile->GetRandom();
    else if ( fProfile->GetParameter(0) == 0  && fProfile->GetParameter(1)!=0) y  = 2 * (gRandom->Rndm() - 0.5) * fProfile->GetParameter(3);
    else                                                                       y  = 0;
    TVector3 p1(x,y,felements[fnum_target].fDistance);
    p1 += fBeam.fOffset;  // beam displacement
    part.fPos += p1;

    Double_t p = 0.0;
    if ( fBeam.fMomRes.Z()==0 ) {
      //double delrnd = (double) ( (int) (gRandom->Rndm() * 13.0) - 6 ) / 100.0;
      double delrnd = (double) ( (int) (gen_func->GetRandom() * 5) ) / 100.0;
      p  = fBeam.fBeamMom * (1.0 + delrnd);
    } else {
      // beam divergence
      p  = fBeam.fBeamMom + fBeam.fBeamMom * (gRandom->Rndm() - 0.5) * 2*fBeam.fMomRes.Z();
    }
    Double_t px = p * (gRandom->Rndm() - 0.5) * 2*fBeam.fMomRes.X();
    Double_t py = p * (gRandom->Rndm() - 0.5) * 2*fBeam.fMomRes.Y();
    Double_t pz = TMath::Sqrt(p*p - px*px - py*py);
    part.fP.SetXYZ(px,py,pz);
    part.fBeamMomSmeared = p;
    part.fName   = "beam";
    part.fDetnum = -1;
    //----------------------------------------

}

Bool_t HBeam::transport(HBeamParticle& part){
  return transport(part, 0);
}

Bool_t HBeam::transport(HBeamParticle& part, int from){

    //----------------------------------------
    // coordinates  x [cm], px/pz [mrad],y [mm], py/pz [mrad], dp to cental p [%]
    //
    //
    //   x_i =  sum (j=0,5) T_ij*x_j + sum (j=0,5)(k=j,5)T_ijk*x_j*x_k
    //
    //
    //

    Double_t state   [5] = { 0,  0, 0,  0, 0};
    Double_t stateold[5] = {
	part.fPos.X()                        *ffromPluto[0],
	TMath::ATan2(part.fP.X(),part.fP.Z())*ffromPluto[1],
	part.fPos.Y()                        *ffromPluto[2],
	TMath::ATan2(part.fP.Y(),part.fP.Z())*ffromPluto[3],
	((part.fBeamMomSmeared-fBeam.fBeamMom)/fBeam.fBeamMom)*ffromPluto[4]
    };

    for(UInt_t el = from; el < felements.size(); el++){
	for(Int_t i = 0; i < 5; i ++){ // state vars
	    // first order Term
	    for(Int_t j = 0; j < 5; j ++){ //
		state[i] += felements[el].Tij[i][j] * stateold[j];
	    }
	    // second order Term
	    for(Int_t j = 0; j < 5; j ++){ //
	    	for(Int_t k = j; k < 5; k ++){ //
	    	    state[i] += felements[el].Tijk[i][j][k] * stateold[j] * stateold[k];
	    	}
	    }
	}

	// propagate and reset
	for(Int_t i = 0; i < 5; i ++) {
	    felements[el].fout[i] = state[i] * ftoPluto[i];
	    state[i] = 0.;
	}
	felements[el].fout[4] = ((part.fBeamMomSmeared-fBeam.fBeamMom)/fBeam.fBeamMom) ; // keep dp
	felements[el].check(); // stats and accepted flag
    }

    return kTRUE;
}

Bool_t HBeam::transport_with_ms(HBeamParticle& part){

  Double_t init_state[5] = { part.x(), part.th(), part.y(), part.ph(), part.dp() };

  part.fName = "beam";
  fhistory.push_back(part);

  part.fName = "ORIG_BEAM_PROFILE";
  transport(part, 0);
  fms_history.push_back(part);

  for (int i=0; i<5; ++i) init_state[i] *= ffromPluto[i]; // convert to transport units before applying transport matrices

  // Propagate to each detector where MS can occur (fXoX0!=0), MS smear and propagate back to the origin
  // Detectors need to be sorted in order of decreasing distance from HADES target. This can be insured by the order in which the corresponding addDetector() functions are called
  // The hits on a detector with MS should be generated and stored before the back propagation!
  for (unsigned int idet=0; idet<fdetectors.size(); ++idet) {

    HBeamElement &det = fdetectors[idet];

    Double_t state[5] = {0.0};
    // propagate state to detector where MS is to be calculated
    for(Int_t i = 0; i < 5; i ++){
      // first order Term
      for(Int_t j = 0; j < 5; j ++){ //
	state[i] += det.Tij[i][j] * init_state[j];
      }
      // second order Term
      for(Int_t j = 0; j < 5; j ++){ //
	for(Int_t k = j; k < 5; k ++){ //
	  state[i] += det.Tijk[i][j][k] * init_state[j] * init_state[k];
	}
      }
    }

    // This state has to be stored as the hit position for this detector (in pluto units!). It is the position before MS happens
    // In createDetectorHits function, the setting of state should be skipped if it has already been done here (test apply_ms()==0)
    for (int i=0; i< 5; ++i) det.fout[i] = state[i] * ftoPluto[i];
    det.check();

    // store this as the state of particle at detector before MS (This is also the position if MS is not requested for a particular detector
    double tmp_state[5] = {0.};
    for (int i=0; i<5; ++i) tmp_state[i] = state[i]*ftoPluto[i];
    part.set_state(tmp_state, fBeam.fBeamMom, det.fDistance);
    part.fName = det.fName;
    part.fAccepted = det.fAccepted;
    fhistory.push_back(part);

    if ( ! det.apply_ms() ) continue; // If MS is not requested for a detector, nothing to do ....

    // Calculate "back propagation" matrix to be used for recalculating initial state after introducing MS smear
    Double_t dp = part.dp() * ffromPluto[4];
    TMatrixD prop(4,4);
    for (int i=0; i<4; ++i) {
      for (int j=0; j<4; ++j) {
	prop[i][j] = det.Tij[i][j] + det.Tijk[i][j][4] * dp;
      }
    }
    prop.Invert();

    // Apply MS smearing to the state vector (position and angles), in cm, mrad
    TLorentzVector tmp;
    tmp.SetXYZM(part.fP.X(),part.fP.Y(),part.fP.Z(), 139.56995*0.001);
    const double sigma_ms = 13.6 * sqrt(det.fXoX0) / ( tmp.P() * tmp.Beta() ); /* mrad */
    //cout << " sigma_ms = " << sigma_ms << endl;
    state[1] += gRandom->Gaus(0.0, sigma_ms);
    state[3] += gRandom->Gaus(0.0, sigma_ms);

    // store this as the state of particle at detector where MS occured. This is the state of the particle after it has MS'ed in the detector
    for (int i=0; i<5; ++i) tmp_state[i] = state[i]*ftoPluto[i];
    part.set_state(tmp_state, fBeam.fBeamMom, det.fDistance);
    part.fName = "POST_" + det.fName + "_MS"; // This is particle state just after MS scattering on a detector, before back propagation
    part.fAccepted = det.fAccepted;
    fhistory.push_back(part);

    // Prepare state vector for multiplication by reverse propation matrix
    state[0] -= (det.Tij[0][4]/*T16*/ * dp + det.Tijk[0][4][4]/*T166*/ * dp * dp );
    state[1] -= (det.Tij[1][4]/*T26*/ * dp + det.Tijk[1][4][4]/*T266*/ * dp * dp );
    state[2] -= (det.Tij[2][4]/*T36*/ * dp + det.Tijk[2][4][4]/*T366*/ * dp * dp );
    state[3] -= (det.Tij[3][4]/*T46*/ * dp + det.Tijk[3][4][4]/*T466*/ * dp * dp );
    //state[0] -= (det.Tij[0][4]/*T16*/ * dp + det.Tijk[0][4][4]/*T166*/ * dp * dp );
    //state[1] -= (det.Tij[1][4]/*T26*/ * dp + det.Tijk[1][4][4]/*T266*/ * dp * dp - gRandom->Gaus(0.0, sigma_ms) );
    //state[2] -= (det.Tij[2][4]/*T36*/ * dp + det.Tijk[2][4][4]/*T366*/ * dp * dp );
    //state[3] -= (det.Tij[3][4]/*T46*/ * dp + det.Tijk[3][4][4]/*T466*/ * dp * dp - gRandom->Gaus(0.0, sigma_ms) );

    TVectorD tvd_state(4, state);
    tvd_state *= prop; // In place multiplication (equiv. to: tvd_state = prop * tvd_state )

    // reinitialize "loop variable" init_state
    for (int i=0; i<4; ++i) init_state[i] = tvd_state[i];

    // Store particle as beam profile after modification to account for MS in current detector
    for (int i=0; i<5; ++i) tmp_state[i] = init_state[i] * ftoPluto[i];
    part.set_state(tmp_state /* in pluto units */, fBeam.fBeamMom, felements[0].fDistance);
    part.fName = Form("POST_%s_MS_BEAM_PROFILE", det.fName.Data());
    fms_history.push_back(part);

    // redo transport for all elements downstream from this detector
    unsigned int next_elment = 0;
    for (; next_elment<felements.size(); ++next_elment) {
      if (det.fDistance <= felements[next_elment].fDistance) break;
    }
    transport(part, next_elment);

  }

  // Accepted if it is within the acceptance of all detectors (at entrance) and all elements
  bool accepted = true;
  for (unsigned int i=0; i<felements.size(); ++i) accepted = accepted && felements[i].fAccepted;
  for (unsigned int i=0; i<fdetectors.size(); ++i) accepted = accepted && fdetectors[i].fAccepted;

  // This is no longer necessary
  //// convert state to pluto units and set HBeamParticle object's initial position and momenta to MS smeared initial condition
  //for (int i=0; i<5; ++i) init_state[i] *= ftoPluto[i];
  //part.set_state(init_state /*in pluto units*/,  fBeam.fBeamMom, felements[fnum_target].fDistance );

  return accepted;

}

Bool_t HBeam::createDetectorHits(HBeamParticle& part)
{
    //calc particles at detector places
    Double_t state   [5] = { 0,  0, 0,  0, 0};
    Double_t p,px,py,pz;
    Bool_t   accepted    = kTRUE;
    Bool_t   acceptedDet = kTRUE;

    for (UInt_t i = 0; i < fdetectors.size(); i ++) {
	accepted = kTRUE;
	for (UInt_t k = 0; k < felements.size() - 1; k ++) {
	    if(!felements[k].fAccepted) {
		accepted = kFALSE;
	    }
	    if (felements[k].fDistance <= fdetectors[i].fDistance &&
		fdetectors[i].fDistance < felements[k+1].fDistance) {

		Double_t frac = (fdetectors[i].fDistance - felements[k].fDistance) / (felements[k+1].fDistance - felements[k].fDistance);

		for(Int_t j = 0; j < 5; j ++){
		    state[j] = felements[k].fout[j] + frac*(felements[k+1].fout[j] - felements[k].fout[j]);
		}

		p  = fBeam.fBeamMom * (1. + state[4]);
		pz = p / sqrt(1.+ tan(state[1])*tan(state[1]) + tan(state[3])*tan(state[3]));
		px = pz * tan(state[1]);
		py = pz * tan(state[3]);

		part.fP  .SetXYZ(px,py,pz);
		part.fPos.SetXYZ(state[0],state[2],fdetectors[i].fDistance);
		part.fDetnum   = fdetectors[i].fId;
		part.fName     = fdetectors[i].fName;
		part.fAccepted = accepted;

		for (int j=0; j<5; ++j) fdetectors[i].fout[j] = state[j];
		fdetectors[i].check();

		if(!fdetectors[i].fAccepted) {
		    part.fAccepted = kFALSE;
		    acceptedDet    = kFALSE;
		}
		fhistory.push_back(part);
	    }
	}
    }
    if(!felements[felements.size()-1].fAccepted) accepted = kFALSE;

    return  accepted&acceptedDet;
}

inline Double_t HBeam::digitize_pos(const Double_t &pos, const Double_t &pitch) {
  //const double pos_digi = ( (double) ( floor( pos/pitch ) ) * pitch ) + pitch/2 + ((gRandom->Rndm()-0.5)*pitch);
  //const double pos_digi = ( (double) ( floor( pos/pitch + 0.5) ) * pitch ) + pitch/2;
  //cout << "pitch= " << pitch << "  pos_in = " << pos << " pos_digi= " << pos_digi << " diff= " << pos_digi-pos << endl;
  double pos_digi = pitch* ( int(pos/fabs(pos))/2.0 + int(pos/pitch) ) ; // + ((gRandom->Rndm()-0.5)*pitch);
  return pos_digi;
}

double pars[4];

double third_order_pol(Double_t x) {
  //cout << "pol= " << pars[0] << "x + " << pars[1] << "x^2 + " << pars[3] << "x^3" << endl;
  return pars[0] + pars[1]*x + pars[2]*x*x + pars[3]*x*x*x;
}

Double_t HBeam::solve_delta(const Double_t &xd1, const Double_t &xd2, const Double_t &phi, int &status) {

  // if phi is non-zero (after first iteration), these transformations will take that into account
  const double _xd1 = xd1 - par.T14d1*phi;
  const double _xd2 = xd2 - par.T14d2*phi;
  const double _T16d1 = par.T16d1 - par.T146d1*phi;
  const double _T16d2 = par.T16d2 - par.T146d2*phi;

  pars[0] = _xd2*par.T12d1 - _xd1*par.T12d2;
  pars[1] = (_xd2*par.T126d1 - par.T12d1*_T16d2) - (_xd1*par.T126d2 - par.T12d2*_T16d1);
  pars[2] = (_T16d1*par.T126d2 + par.T12d2*par.T166d1) - (_T16d2*par.T126d1 + par.T12d1*par.T166d2);
  pars[3] = (par.T126d2*par.T166d1 - par.T126d1*par.T166d2);

  bool debug_3rdO_pol = false;
  if (debug_3rdO_pol) {
    TF1 fp3("pol3","pol3",-20,20);
    for (int i=0; i<3; ++i) fp3.SetParameter(i,pars[i]);
    TCanvas tc("tc","tc",1000,1000);
    tc.cd();
    gPad->SetGridy();
    fp3.Draw();
    tc.Print("func.png");
    string tmp;
    cout << "Hit anything to continue";
    cin >> tmp;
  }

  //cout << "==================================================================" << endl;
  //cout << "Input polynomial = " << pars[0] << "x + " << pars[1] << "x^2 + " << pars[3] << "x^3" << endl;

  // Numerical root finding method from: http://root.cern.ch/drupal/content/root-finder-algorithms

  Double_t del_max = 8.0;
  ROOT::Math::Functor1D f1D(&third_order_pol);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kBRENT);
  //ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
  //ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_FALSE_POS);
  //ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BRENT);

  bool status_solved = false;
  while (!status_solved) {
    del_max += 1.0;
    if (del_max>11) break;
    rf.SetFunction( f1D, -del_max, del_max);
    rf.Solve();
    status_solved = (rf.Status() == 0);
  }

  //cout << "Roots of 3rd order equation: " << rf.Root() << endl;
  if (!status_solved) {
    cout << "Still trouble solving up to del_max= " << del_max << endl;
  }

  status = rf.Status();

  return rf.Root();

  //const double pars[4] = {
  //  _xd2*par.T12d1 - _xd1*par.T12d2,
  //  (_xd2*par.T126d1 - par.T12d1*_T16d2) - (_xd1*par.T126d2 - par.T12d2*_T16d1),
  //  (_T16d1*par.T126d2 + par.T12d2*par.T166d1) - (_T16d2*par.T126d1 + par.T12d1*par.T166d2),
  //  (par.T126d2*par.T166d1 - par.T126d1*par.T166d2)
  //};

  //// Numerical root finding method from: http://root.cern.ch/drupal/content/root-finder-algorithms
  //// SetParameter function of TF1 is incredibly slow. Use TFormula to enter the parameters directly
  //const Double_t del_max = 10.0;
  //TF1 f("ThirdOrdPol",Form("%f+%f*x+%f*x*x+%f*x*x*x",pars[0],pars[1],pars[2],pars[3]),-del_max,del_max);
  //ROOT::Math::WrappedTF1 wf1(f);
  //ROOT::Math::BrentRootFinder brf;
  //brf.SetFunction( wf1, -del_max, del_max );
  //brf.Solve();
  ////cout << "Roots of 3rd order equation: " << brf.Root() << endl;
  //return brf.Root();

}

double HBeam::solve_theta(const Double_t &xd1, const Double_t &xd2, const Double_t &phi, const Double_t &delta) {

  // if phi is non-zero (after first iteration), these transformations will take that into account
  const double _xd1 = xd1 - par.T14d1*phi;
  const double _T16d1 = par.T16d1 - par.T146d1*phi;
  const double theta_d1 = (_xd1 - _T16d1 * delta - par.T166d1 * delta * delta) / ( par.T12d1 + par.T126d1 * delta);

  //double _xd2 = xd2 - par.T14d2*phi;
  //double _T16d2 = par.T16d2 - par.T146d2*phi;
  //double theta_d2 = (_xd2 - _T16d2 * delta - par.T166d2 * delta * delta) / ( par.T12d2 + par.T126d2 * delta);

  //cout << "element[0] th= " << felements[0].fout[1]*ffromPluto[1] ;
  //cout << " Theta_d1= " << theta_d1 << " Theta_d2= " << theta_d2 << endl;
  //assert( fabs( ( theta_d1 - theta_d2)/theta_d2 ) < 1 );

  return theta_d1;
}

void HBeam::solve_phi_and_y(const Double_t &yd1, const Double_t &yd2, const Double_t &th, const Double_t &del, Double_t &ph, Double_t &y) {

  const double a1 = par.T33d1 + par.T336d1*del;
  const double a2 = par.T33d2 + par.T336d2*del;
  const double b1 = par.T34d1 + par.T346d1*del;
  const double b2 = par.T34d2 + par.T346d2*del;
  const double g1 = yd1 - par.T32d1*th - par.T36d1*del - par.T366d1*del*del;
  const double g2 = yd2 - par.T32d2*th - par.T36d2*del - par.T366d2*del*del;
  y = (g1*b2 - g2*b1)/(b2*a1 - b1*a2);
  ph = (g1 - a1*y) / b1;

}

void HBeam::transport_to(const HBeamElement &det, const vector<Double_t> &state_in /* Pluto Units Required*/, vector<Double_t> &state_out /* Answer will be in pluto units too*/) {

  for (int i=0; i<5; ++i) state_out[i] = 0.0; // for safety
  vector<Double_t> _state_in(5);
  for (int i=0; i<5; ++i) _state_in[i] = state_in[i] * ffromPluto[i];

  for(Int_t i = 0; i < 5; i ++){ // state vars
    // first order Term
    for(Int_t j = 0; j < 5; j ++){ //
      state_out[i] += det.Tij[i][j] * _state_in[j];
    }
    // second order Term
    for(Int_t j = 0; j < 5; j ++){ //
      for(Int_t k = j; k < 5; k ++){ //
	state_out[i] += det.Tijk[i][j][k] * _state_in[j] * _state_in[k];
      }
    }
  }
  for (int i=0; i<5; ++i) state_out[i] *= ftoPluto[i];
}

int HBeam::solution_tdr(Double_t *x, Double_t *y, vector<Double_t> &state) {

  // transport units
  // Solve for momentum and position using procedure described in TDR at 1) prod. target 2) Hades 3) Diamond Det, store in fsolution with appropriate name
  // [0]->x [1]->the [2]->y, [3]->phi, [4]->Delta
  // Assumption of the method: Zero object size in horizontal, so state[0] remains at 0.0

  const Double_t tolerance = 1e-4; /// required accuracy
  const Int_t max_niter = 2;
  vector<Double_t> prev(5,0.0);
  int niter = 0;
  bool converged = false;
  int status = 0;
  while (niter<max_niter && !converged ) {
    for (int i=0; i<5; ++i) prev[i] = state[i];

    //state[4] = solve_delta(xd1, xd2, state[3]);
    //state[1] = solve_theta(xd1, xd2, state[3], state[4]);
    //solve_phi_and_y(yd1, yd2, state[1], state[4], state[3], state[2]);

    //cout << "===================== STEP " << niter << "=======================" <<endl;
    //cout << "state0 = " << state[0] << " state1 = " << state[1] << " state2 = " << state[2] << " state3 = " << state[3] << " state4 = " << state[4] << endl;
    //cout <<  "x[0]= " << x[0] << " x[1]= " << x[1] << " phi= " << state[3] << endl;

    state[4] = solve_delta(x[0], x[1], state[3], status);
    state[1] = solve_theta(x[0], x[1], state[3], state[4]);
    solve_phi_and_y(y[0], y[1], state[1], state[4], state[3], state[2]);

    //cout << "============== END OF STEP " << niter << "=======================" <<endl;

    converged = true;
    for (int i=1; i<5; ++i) converged = converged && ( ((state[i] - prev[i])/state[i]) < tolerance );
    ++niter;
  }

  return status;

}

int HBeam::solution_minuit(vector<Double_t> &state /*Everyting in Transport units here*/ ) {

  // Does this have to be done for each fit?
  //arglist[0] = 1;
  //ierflg = 0;
  //gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // Set starting values and step sizes for parameters
  //static Double_t vstart[4] = {3, 1 , 0.1 , 0.01};
  Double_t step[5] = {0.0};
  const Double_t nominal[5] = { 0.5, 10, 0.05, 50, 6.};
  const Double_t grid_size = 1000;
  for (int i=0; i<5; ++i) { step[i]= nominal[i]/grid_size; }

  //cout << "Step0 = " << step[0] << " Step1 = " << step[1] << " Step2 = " << step[2] << " Step3 = " << step[3] << " Step4 = " << step[4] << endl;

  // Solution from manual iteration will be used for initial value
  // Initial step size will be taken a constant fraction of the nominal max value (10th to start with)
  gMinuit->mnparm(0, "x0",  state[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "th0", state[1], step[1], 0, 0, ierflg);
  gMinuit->mnparm(2, "y0",  state[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "ph0", state[3], step[3], 0, 0, ierflg);
  //const Double_t del_max = 10.0; // 8.0
  //gMinuit->mnparm(4, "del", state[4], step[4], -del_max, del_max, ierflg);
  gMinuit->mnparm(4, "del", state[4], step[4], 0, 0, ierflg);

  gMinuit->FixParameter(1);
  gMinuit->FixParameter(3);

  gMinuit->FixParameter(0);
  gMinuit->FixParameter(2);

  //gMinuit->FixParameter(4);

  // Now ready for minimization step
  //arglist[0] = 500;
  //arglist[1] = 1.;
  arglist[0] = 5000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
  //gMinuit->mnexcm("SEE", arglist ,2,ierflg);

  //cout << "erflg = " << ierflg << endl;

  //// Print results what does this serve?
  //Double_t amin,edm,errdef;
  //Int_t nvpar,nparx,icstat;
  //gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  double init[5] = {0.0};
  for (int i=0; i<5; ++i) { init[i] = state[i]; }
  //double diff_Del = state[4];
  bool verb= false;
  if (verb) {
    if (ierflg!=0) {
      cout << "ERROR  state_bef: " << state[0] << " " << state[1] << " " << state[2] << " " << state[3] << " " << state[4] << endl;
    } else {
      cout << "OKOKOK state_bef: " << state[0] << " " << state[1] << " " << state[2] << " " << state[3] << " " << state[4] << endl;
    }
  }

  // here is where the real fitting happens
  Double_t error = 0.0;
  for (int i=0; i<5; ++i) gMinuit->GetParameter(i, state[i], error);

  for (int i=0; i<4; ++i) {
    if (init[i] != state[i]) {cout << "WTF init[" << i << "]("  << init[i] << ") != " << " state[" << i << "](" << state[i] << ")" << endl;  }
  }

  if (verb) {
    if (ierflg!=0) {
      cout << "ERROR  state_aft: " << state[0] << " " << state[1] << " " << state[2] << " " << state[3] << " " << state[4] << endl;
    } else {
      cout << "OKOKOK state_aft: " << state[0] << " " << state[1] << " " << state[2] << " " << state[3] << " " << state[4] << endl;
    }
  }

  return 0;

}

void HBeam::solve_state() {

  //cout << "========================================== " << endl;
  // fetch hit position in each detector whose hit position is used in the reconstrction
  for (unsigned int i=0; i<det_idx.size(); ++i) {

    bool digi = fdetectors[det_idx[i]].fDigi;
    Double_t pitch = fdetectors[det_idx[i]].fPitch*ffromPluto[0];  // cm

    Double_t x_det = fdetectors[det_idx[i]].fout[0]*ffromPluto[0];  // cm
    Double_t y_det = fdetectors[det_idx[i]].fout[2]*ffromPluto[2];  // cm

    x_mes[i] = digi ? digitize_pos( x_det , pitch ): x_det;
    y_mes[i] = digi ? digitize_pos( y_det , pitch ): y_det;
    //cout << "det " << det_idx[i] << "  " << fdetectors[det_idx[i]].fName << " x_det = " << x_det
    //	 << " pitch= " << pitch << " x_det_digi= " << x_mes[i] << endl;
    x_mes_e[i] = y_mes_e[i] = pitch/TMath::Sqrt(12);
  }

  HBeamParticle part(fBeam);
  vector<Double_t> state_targ(5,0.0);

  int state = solution_tdr(&x_mes[0], &y_mes[0], state_targ);
  part.fStatus = state;
  store_solution("TDR", state_targ, part);

  double init[5] = {0.0};
  for (int i=0; i<5; ++i) { init[i] = state_targ[i]; }
  solution_minuit(state_targ);
  store_solution("MIN", state_targ, part);

  // Save the "dgitized" positions in detector 1 and 2
  state_targ[1] = -9999.0;
  state_targ[3] = -9999.0;
  state_targ[4] = -9999.0;
  for (int idet=0; idet<3; idet++) {
    state_targ[0] = x_mes[idet] * ftoPluto[0];
    state_targ[2] = y_mes[idet] * ftoPluto[2];
    part.set_state(&state_targ[0], fBeam.fBeamMom, fdetectors[det_idx[idet]].fDistance);
    part.fName = Form("DIGI_POS_%s",fdetectors[det_idx[idet]].fName.Data());
    fsolution.push_back(part);
  }

}

void HBeam::store_solution(const char *tag, vector<Double_t> &state_targ, HBeamParticle& part) {
  // switch state to pluto units to store state at target and transport to start detector and hades
  for (int i=0; i<5; ++i) state_targ[i]*= ftoPluto[i];

  part.set_state(&state_targ[0], fBeam.fBeamMom, felements[fnum_target].fDistance);
  part.fName = Form("%s_SOLUTION_BEAM_INITIAL",tag);
  fsolution.push_back(part);

  vector<Double_t> state_diam(5,0.0);
  transport_to(fdetectors[det_idx[2]], state_targ, state_diam);
  part.set_state(&state_diam[0], fBeam.fBeamMom, fdetectors[det_idx[2]].fDistance);
  part.fName = Form("%s_SOLUTION_DIAM",tag);
  fsolution.push_back(part);

  vector<Double_t> state_hades(5,0.0);
  transport_to(fdetectors[det_idx[3]], state_targ, state_hades);
  part.set_state(&state_hades[0], fBeam.fBeamMom, fdetectors[det_idx[3]].fDistance);
  part.fName = Form("%s_SOLUTION_HAD",tag);
  fsolution.push_back(part);

  // switch state back to transport units for next step
  for (int i=0; i<5; ++i) state_targ[i] *= ffromPluto[i];
}

vector <HBeamParticle>& HBeam::newParticle()
{
    HBeamParticle part(fBeam);
    fhistory.clear();
    fms_history.clear();
    fsolution.clear();

    createBeam(part);         // beam profile + smearing

    Bool_t accepted = transport_with_ms(part);    // MS effect added for each detector registered with non zero thinkness. Creates and stores detector hits

    //part.fName = "beam";
    //fhistory.push_back(part);
    //transport(part);          // transport particle through beamline
    //Bool_t accepted = createDetectorHits(part); // fill detectors

    // calc particle at 0,0,0 ideal target
    Double_t* out = &felements[fnum_target].fout[0];
    Double_t p,px,py,pz;
    p  = fBeam.fBeamMom * (1. + out[4]);
    pz = p / sqrt(1.+ tan(out[1])*tan(out[1]) + tan(out[3])*tan(out[3]));
    px = pz * tan(out[1]);
    py = pz * tan(out[3]);

    part.fP  .SetXYZ(px,py,pz);
    part.fPos.SetXYZ(out[0],out[2],felements[fnum_target].fDistance);
    part.fDetnum = -1;
    part.fName   = felements[fnum_target].fName;
    part.fAccepted = accepted;
    fhistory.push_back(part);

    solve_state();

    return fhistory;
}


void HBeam::printBeamLine(Bool_t trans)
{
    cout<<"############################################################################"<<endl;
    cout<<"BEAMLINE ELEMENTS :"<<endl;
    for(UInt_t i = 0; i < felements.size(); i ++){
	felements[i].print(trans);
    }
    cout<<"############################################################################"<<endl;
}
void HBeam::printBeamProfile()
{
    cout<<"############################################################################"<<endl;
    cout<<"BEAM PROFILE :"<<endl;
    for(Int_t i = 0 ; i < fProfile->GetNpar();i ++){
	cout<<setw(20)<<fProfile->GetParName(i)<<" : "<<setw(12)<<fProfile->GetParameter(i)<<endl;
    }
    fBeam.print(kFALSE);
    cout<<"############################################################################"<<endl;
}
void HBeam::printDetectors()
{
    cout<<"############################################################################"<<endl;
    cout<<"DETECTORS :"<<endl;
    for(UInt_t i = 0; i < fdetectors.size(); i ++){
	fdetectors[i].print(kTRUE);
    }
    cout<<"############################################################################"<<endl;
}

Bool_t HBeam::initBeamLine(TString filename,Int_t targetElementNum,Bool_t debug)
{
    //opens the ASCII description of the beam line


    cout<<"HBeam:: reading input!"<<endl;
    Int_t nelem = 0;


    std::ifstream input;
    input.open(filename.Data());

    if(input.fail()){
	cout<<"HBeam:: could not open input!"<<endl;
	return kFALSE;
    }


    Char_t str[1000];
    Int_t  nread = 1000;
    Bool_t readfile = kTRUE;
    //----------------------------
    // local dummy vars
    Double_t aa;
    Int_t ind1;
    Int_t ind23;
    Double_t distance;
    //----------------------------


    while (readfile && !input.fail()) {

	input>>distance;
	if(input.peek()!='\n') input.ignore(nread,'\n');
	if(input.fail() && nelem > 0 ) continue;
	if(input.fail() ){ cout<<"HBeam:: "<<nelem<<" failed reading distance"<<endl; return kFALSE;}
	distance *= 1000;  // [m] -> [mm]
	if(debug) cout<<nelem<<" ------------------------"<<distance<<endl;
	for(Int_t i = 0; i < 5; i++) {
	    input.ignore(nread,'\n');
	}
	if(input.fail() ) {
	    cout<<"HBeam:: "<<nelem<<" failed after first block "<<endl;
	    return kFALSE;
	}
	input.getline(str,nread,'\n');  // *TRANSFORM 1*
	if(debug) cout<<"HBeam:: "<<"1st "<<str<<endl;

	HBeamElement e(Form("Element[%i]",nelem),distance,nelem);
	// first order transform
	Int_t ind = 0;
	for(Int_t i = 0; i < 6; i++){ // lines
	    if(i == 4) {
		input>>aa>>aa>>aa>>aa>>aa>>aa;
		continue;
	    }
	    if(i > 4) ind = i-1;
	    else      ind = i;
	    input>>e.Tij[ind][0]>>e.Tij[ind][1]>>e.Tij[ind][2]>>e.Tij[ind][3]>>aa>>e.Tij[ind][4];
	}
	if(debug) e.printTij();

	if(input.fail() ) {
	    cout<<"HBeam:: "<<nelem<<" failed after first order "<<endl;
	    return kFALSE;
	}
	if(input.peek() =='\n') input.ignore(nread,'\n');
	input.getline(str,nread,'\n');   // *2ND ORDER TRANSFORM*

	if(debug) cout<<"HBeam:: "<<"2nd : "<<str<<endl;
	//C
	//C     second order transform
	//C

	for(Int_t i = 0; i < 5; i++){ // blocks i
	    input>>ind1>>ind23>>e.Tijk[i][0][0];
	    input>>ind1>>ind23>>e.Tijk[i][0][1]>>ind1>>ind23>>e.Tijk[i][1][1];
	    input>>ind1>>ind23>>e.Tijk[i][0][2]>>ind1>>ind23>>e.Tijk[i][1][2]>>ind1>>ind23>>e.Tijk[i][2][2];
	    input>>ind1>>ind23>>e.Tijk[i][0][3]>>ind1>>ind23>>e.Tijk[i][1][3]>>ind1>>ind23>>e.Tijk[i][2][3]>>ind1>>ind23>>e.Tijk[i][3][3];
	    input>>ind1>>ind23>>aa             >>ind1>>ind23>>aa             >>ind1>>ind23>>aa             >>ind1>>ind23>>aa             >>ind1>>ind23>>aa;
	    input>>ind1>>ind23>>e.Tijk[i][0][4]>>ind1>>ind23>>e.Tijk[i][1][4]>>ind1>>ind23>>e.Tijk[i][2][4]>>ind1>>ind23>>e.Tijk[i][3][4]>>ind1>>ind23>>aa>>ind1>>ind23>>e.Tijk[i][4][4];
	}

	//// Important fix, 5th block doesn't refer to delta component index
	//// The delta block inices are all zoro, meaning delta doesn't change
	for (int jj=0; jj<5; jj++) {
	  for (int kk=0; kk<5; kk++) {
	    e.Tijk[4][jj][kk] = 0.0;
	  }
	}

	if(debug) e.printTijk();

	if(input.fail() ) {
	    cout<<"HBeam:: "<<nelem<<" failed after second order "<<endl;
	    return kFALSE;
	}

	felements.push_back(e);
	input.ignore(nread,'\n');

	if(input.fail() ) return kFALSE;

	if (!input.eof())
	    nelem++;
	else
	    readfile = kFALSE;

    } //end nelem

    input.close();
    cout<<"HBeam:: finished reading input!"<<endl;
    return setTargetElement(targetElementNum);

}


#endif

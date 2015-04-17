#ifndef BREMPIDREADER_H
#define BREMPIDREADER_H

#include <vector>

#include "FairTask.h"
#include "TObject.h"


#include "TString.h"
#include "TClonesArray.h"

class PndPidCandidate;

class PndPidBremCorrected4Mom;
class PndEmcCluster;
class PndEmcBump;

class TFile;
class TTree;
class TString;

class FairMCEventHeader;

// Path of file:
//  ----- $pandaroot/pid/PidCorr

class BremPidReader: public FairTask
{

 public:

        BremPidReader();

	// Destructor
	virtual ~BremPidReader();

	// Methods
	/** Virtual method Init **/
	virtual InitStatus Init();

	/** Virtual method Exec **/
	virtual void Exec(Option_t* opt);

        void SetStorageOfData(Bool_t p = kTRUE) {fPersistance=p;};

        virtual void FinishTask();

	void set_output_name(TString arg) { output_name = arg; }

 private:

        PndPidBremCorrected4Mom* AddBremCorrected4Mom();

	void print_mc();
	void print_cands();
	void print_vect(std::vector<int> &);

        //double GetSepPhotonE(PndPidCandidate *, int&);
        //double GetMergPhotonE(PndPidCandidate *, int&);

	void fill_bump_list(std::vector<std::vector<int> >&);
	void brem_matching(std::vector<std::vector<int> >&, std::vector<std::vector<int> >&, std::vector<int>&, std::vector<int>&, std::vector<int>&, std::vector<int>&);
	void find_ancestory(const int&, std::vector<int>&);
	//void get_mc_brem_photons(const int&, std::vector<double>&, std::vector<double>&, std::vector<double>&,
	//			 std::vector<double>&, std::vector<double>&, std::vector<std::vector<int> >&);
	//double GetSepPhotonE_fromBumps(PndPidCandidate*, double&, double&, std::vector<int>&,
	//			       std::vector<double>&, std::vector<double>&, std::vector<double>&,
	//			       std::vector<double>&, std::vector<std::vector<int> >&);
	void get_mc_brem_photons(const int&, std::vector<std::vector<int> >&);
	void GetSepPhotonE_fromBumps(PndPidCandidate*, double&, double&, double&, std::vector<std::vector<int> >&);
	void GetMergPhotonE(PndPidCandidate *, double&, double&, double&, double&, double&, double&);
	double corrected_mom(const double&);

	int nEvt;

	/** Input array of PndEmcClusters **/
	TClonesArray* fBumpArray;
	TClonesArray* fClusterArray;
        TClonesArray* fDigiArray;
        TClonesArray* fHitArray;

        TClonesArray* fPhiBumpArray;

        TClonesArray* fChargedCandidateArray;
        TClonesArray* fNeutralCandidateArray;

        TClonesArray* fBremCorrected4MomArray;

	TClonesArray* fMcArray;

	FairMCEventHeader* fMCHeader;

        double fRecMomOfEle;
        double fRecThetaOfEle;
        double fRecPhiOfEle;
        int fCharge;

	Bool_t fPersistance; // switch to turn on/off storing the arrays to a file
	// Data members

	Double_t radMaxTracking_cm;

        BremPidReader(const BremPidReader& L);
        BremPidReader& operator= (const BremPidReader&) {return *this;};

	TString output_name;
	TFile* f;
	TTree* t;
	static const int nch_max = 10;
	static const int nNeutMax = 100;
	int nch;
	int charge[nch_max];
	float mom_mc[nch_max];
	float mom_rec[nch_max];
	float mom_cor[nch_max];
	float mom_wcor[nch_max];
	float mom_sep[nch_max];
	float mom_mrg[nch_max];
	float mom_stored[nch_max];

	float mom_sep_w[nch_max];
	float mom_sep_w_bf[nch_max];
	float mom_mrg_w[nch_max];
	float mom_mrg_w_bf[nch_max];
	float mom_mrg_pc[nch_max];
	float mom_mrg_w_pc[nch_max];
	float mom_mrg_w_bf_pc[nch_max];
	float mom_wbfcor[nch_max];
	float mom_wbfcor_mw_bf[nch_max];
	float mom_wbfcor_mw_bf_pc[nch_max];

	float phi[nch_max];
	float the[nch_max];
	float phi_mc[nch_max];
	float the_mc[nch_max];
	int pdg_mc[nch_max];
	int nphot_sep[nch_max];
	int nphot_mrg[nch_max];
	int nphot_mrg_pc[nch_max];
	int is_prim[nch_max];

	int _nmcb[nch_max];
	int imcb_s[nch_max];
	int imcb_e[nch_max];

	static const int nmcb_max = 200;
	int nmcb;
	float mcb_phi[nmcb_max];
	float mcb_the[nmcb_max];
	float mcb_ene[nmcb_max];
	float mcb_rad[nmcb_max];
	float mcb_zed[nmcb_max];
	int mcb_match[nmcb_max];
	int mcb_score[nmcb_max];
	int mcb_match_ab[nmcb_max];
	int mcb_score_ab[nmcb_max];

	int _nsb[nch_max];
	int isb_s[nch_max];
	int isb_e[nch_max];

	static const int nsb_max = 100;
	int nsb;
	int sb_idx[nsb_max];
	float sb_phi[nsb_max];
	float sb_the[nsb_max];
	float sb_ene[nsb_max];
	float sb_rcalc[nsb_max];
	int sb_match[nsb_max];
	int sb_score[nsb_max];

	// All bumps list
	static const int nab_max = 100;
	int nab;
	float ab_phi[nab_max];
	float ab_the[nab_max];
	float ab_ene[nab_max];
	int ab_isb[nab_max];
	int ab_ich[nab_max];
	int ab_match[nsb_max];
	int ab_score[nsb_max];

	int nNeutCand;
	int nChCand;
	int nMcTrack;

	ClassDef(BremPidReader,1);

};

#endif

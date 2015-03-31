#ifndef BREMPIDREADER_H
#define BREMPIDREADER_H

#include <vector>
#include <map>

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

	void print_cands();

        double GetSepPhotonE(PndPidCandidate *, int&);
        double GetSepPhotonE_fromBumps(PndPidCandidate *, int&);
        double GetMergPhotonE(PndPidCandidate *, int&);

        void GetEmcPhiBumpList(int iClust);

	void get_earliest_brem(int &n, double &r, double &e);

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

        Double_t fSepPhotonE;
        Double_t fMergPhotonE;

	std::vector<int> fSepClustId;

	std::vector<PndEmcBump*> fEmcPhiBumpList;

	Bool_t fPersistance; // switch to turn on/off storing the arrays to a file
	// Data members

        BremPidReader(const BremPidReader& L);
        BremPidReader& operator= (const BremPidReader&) {return *this;};

	TString output_name;
	TFile* f;
	TTree* t;
	static const int nChMax = 10;
	static const int nNeutMax = 100;
	int nChCand;
	int ch_charge[nChMax];
	float ch_mom_mc[nChMax];
	float ch_mom_rec[nChMax];
	float ch_mom_cor[nChMax];
	float ch_mom_sep[nChMax];
	float ch_mom_mrg[nChMax];
	float ch_mom_stored[nChMax];
	float ch_phi[nChMax];
	float ch_the[nChMax];

	float brem_rad[nChMax];
	float brem_e[nChMax];
	int nbrem_trk[nChMax];

	int ch_nphot_sep[nChMax];
	int ch_nphot_mrg[nChMax];
	int ch_is_prim[nChMax];

	int nNeutCand;

	ClassDef(BremPidReader,1);

};

#endif

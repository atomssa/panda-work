#ifndef PID_BREMCORRECTORNT_H
#define PID_BREMCORRECTORNT_H

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

class FairMCEventHeader;

// Path of file:
//  ----- $pandaroot/pid/PidCorr

class PndPidBremCorrectorNT: public FairTask
{

 public:

        PndPidBremCorrectorNT();

	// Destructor
	virtual ~PndPidBremCorrectorNT();

	// Methods
	/** Virtual method Init **/
	virtual InitStatus Init();

	/** Virtual method Exec **/
	virtual void Exec(Option_t* opt);

        void SetStorageOfData(Bool_t p = kTRUE) {fPersistance=p;};

        virtual void FinishTask();

 private:

        PndPidBremCorrected4Mom* AddBremCorrected4Mom();

	void print_cands();

        double GetSepPhotonE(PndPidCandidate *, int&);
        double GetMergPhotonE(PndPidCandidate *, int&);

        void GetEmcPhiBumpList(int iClust);

	int fMode;

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

        PndPidBremCorrectorNT(const PndPidBremCorrectorNT& L);
        PndPidBremCorrectorNT& operator= (const PndPidBremCorrectorNT&) {return *this;};

	int nChCand;
	int nNeutCand;

	ClassDef(PndPidBremCorrectorNT,1);

};

#endif

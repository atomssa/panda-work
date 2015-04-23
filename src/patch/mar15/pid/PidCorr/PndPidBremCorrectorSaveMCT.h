#ifndef PID_BREMCORRECTORMCT_H
#define PID_BREMCORRECTORMCT_H

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

// Path of file:
//  ----- $pandaroot/pid/PidCorr

class PndPidBremCorrectorSaveMCT: public FairTask
{

 public:

        PndPidBremCorrectorSaveMCT();

	// Destructor
	virtual ~PndPidBremCorrectorSaveMCT();

	// Methods
	/** Virtual method Init **/
	virtual InitStatus Init();

	/** Virtual method Exec **/
	virtual void Exec(Option_t* opt);

        void SetStorageOfData(Bool_t p = kTRUE) {fPersistance=p;};


        virtual void FinishTask() {};

 private:

        PndPidBremCorrected4Mom* AddBremCorrected4Mom();

        double GetSepPhotonE(PndPidCandidate *, std::vector<Int_t>&);
        double GetMergPhotonE(PndPidCandidate *, std::vector<Int_t>&);

        void GetEmcPhiBumpList(int iClust);

	/** Input array of PndEmcClusters **/
	TClonesArray* fBumpArray;
	TClonesArray* fClusterArray;


        TClonesArray* fPhiBumpArray;

        TClonesArray* fChargedCandidateArray;
        TClonesArray* fNeutralCandidateArray;

        TClonesArray* fBremCorrected4MomArray;

        double fRecMomOfEle;
        double fRecThetaOfEle;
        double fRecPhiOfEle;
        int fCharge;

        Double_t fSepPhotonE;
        Double_t fMergPhotonE;

	std::vector<PndEmcBump*> fEmcPhiBumpList;

	Bool_t fPersistance; // switch to turn on/off storing the arrays to a file
	// Data members

        PndPidBremCorrectorSaveMCT(const PndPidBremCorrectorSaveMCT& L);
        PndPidBremCorrectorSaveMCT& operator= (const PndPidBremCorrectorSaveMCT&) {return *this;};

	ClassDef(PndPidBremCorrectorSaveMCT,1);

};

#endif

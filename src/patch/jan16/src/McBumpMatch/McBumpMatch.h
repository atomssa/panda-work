#ifndef MCBUMPMATCH_H
#define MCBUMPMATCH_H

#include <vector>

#include "FairTask.h"
#include "TObject.h"


#include "TString.h"
#include "TClonesArray.h"

class PndEmcCluster;
class PndEmcBump;

class TString;

class McBumpMatch: public FairTask
{

 public:

        McBumpMatch();

	// Destructor
	virtual ~McBumpMatch() {}

	// Methods
	/** Virtual method Init **/
	virtual InitStatus Init();

	/** Virtual method Exec **/
	virtual void Exec(Option_t* opt);

        virtual void FinishTask() {}

 private:

	void print_mc();
	void print_cands();
	void print_vect(std::vector<int> &);

	bool split_detector();
	void find_ancestory(const int&, std::vector<int>&);
	void find_contributors(PndEmcBump *bump, std::vector<int>&);

	void mc_match_bumps_main_contrib(std::vector<PndEmcBump*>, std::vector<int>&);
	void get_mc_brem_photons(const int&, std::vector<int>&);

	int nEvt;
	int nMcTrack;
	int nChCand;
	int nNeutCand;

	double radMaxTracking_cm;

	/** Input array of PndEmcClusters **/
	TClonesArray* fBumpArray;
	TClonesArray* fClusterArray;
        TClonesArray* fDigiArray;
        TClonesArray* fHitArray;

        TClonesArray* fChargedCandidateArray;
        TClonesArray* fNeutralCandidateArray;

	TClonesArray* fMcArray;

        McBumpMatch(const McBumpMatch& L);
        McBumpMatch& operator= (const McBumpMatch&) {return *this;};

	ClassDef(McBumpMatch,1);

};

#endif

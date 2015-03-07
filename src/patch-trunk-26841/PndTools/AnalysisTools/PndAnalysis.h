#ifndef PNDANALYSIS_H
#define PNDANALYSIS_H 1

#include <string>
#include <vector>

#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

#include "RhoCandList.h"
#include "FairRootManager.h"
#include "FairTrackParP.h"
class TChain;
class TClonesArray;
class TBranch;

class TDatabasePDG;
class FairMCEventHeader;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndMCTrack;
//class PndEventInfo;

class TGeant3;

class PndAnalysis
{
  public:
    PndAnalysis(TString tname1="", TString tname2="",TString algnamec="PidAlgoIdealCharged", TString algnamen="PidAlgoIdealNeutral");
    virtual ~PndAnalysis();

    void Rewind();
    void Reset() { Rewind(); }
    Int_t  GetEvent(Int_t n=-1);
    void GetEventInTask();
    FairMCEventHeader* GetEventHeader();
    Int_t  GetEntries();
    Bool_t FillList(RhoCandList& l, TString listkey="All", TString pidTcaNames="");
    void SetVerbose(Int_t level) {fVerbose = level;}
    //void SetPidChargedName(TString s) {fChargedPidName = s;}
    //void SetPidNeutralName(TString s) {fNeutralPidName = s;}
    //void SetTracksName(TString s) {fTracksName = s;}
    //void SetTracksName2(TString s) {fTracksName2 = s;}
    Bool_t PropagateToIp(RhoCandidate* cand);
    Bool_t PropagateToZAxis(RhoCandidate* cand);
    Bool_t PropagateToPoint(RhoCandidate* cand, TVector3 mypoint);
    Bool_t PropagateToPlane(RhoCandidate* cand, TVector3 origin, TVector3 dj, TVector3 dk);

    FairTrackParP GetFirstPar(RhoCandidate* cand);
    Bool_t ResetCandidate(RhoCandidate* cand);
    Bool_t ResetDaughters(RhoCandidate* cand);
    Bool_t McTruthMatch(RhoCandidate* cand, Int_t level = 2, bool verbose=false);
    Int_t McTruthMatch(RhoCandList& list, Int_t level = 2, bool verbose=false);
	// this also allows match, when some soft photons are missing in decay
	void McMatchAllowPhotos(int maxn=1, double thresh=0.05){fPhotosMax=maxn;fPhotosThresh=thresh;}

    //FIXME: This is an aweful solution to access the correct
    //track array from a fitter object. [R.K.03'11]
    TClonesArray* GetTrackArrayPointer() const {return fTracks;};

    Bool_t Propagator(int mode, FairTrackParP& tStart, RhoCandidate* cand,
                      TVector3 point=TVector3(0,0,0), Bool_t skipcov=kFALSE, Bool_t overwrite=kFALSE, TVector3 planej=TVector3(1,0,0), TVector3 planek=TVector3(0,1,0));

  private:

    void Init();
    void BuildMcCands();
    Bool_t GetMcCandList(RhoCandList& l);
    void ReadRecoCandidates();
	void Cleanup();
	void ReadCandidates();

    TClonesArray* ReadTCA(TString tcaname);

    Bool_t MctMatch(RhoCandidate& c, RhoCandList& mct, Int_t level = 2, bool verbose=false){ return MctMatch(&c,mct,level,verbose); };
    Bool_t MctMatch(RhoCandidate* c, RhoCandList& mct, Int_t level = 2, bool verbose=false);
    // Private Member Variables

    FairRootManager*   fRootManager;
    PndAnaPidSelector* fPidSelector;
    PndAnaPidCombiner* fPidCombiner;
    Int_t              fEvtCount;
    Int_t              fChainEntries;
    Bool_t             fEventRead;
    Bool_t             fBuildMcCands;
    Int_t              fVerbose;
    Int_t              fPhotosMax;
    Double_t           fPhotosThresh;

    std::map<int,RhoCandidate*> fMcPresenceMap;

    TDatabasePDG* fPdg;

    TClonesArray* fChargedCands;
    TClonesArray* fBremCorr;
    TClonesArray* fNeutralCands;
    TClonesArray* fChargedProbability;
    TClonesArray* fNeutralProbability;
    TClonesArray* fTracks;
    TClonesArray* fTracks2;
    TClonesArray* fMcCands;
    TClonesArray* fMcTracks;

    RhoCandList fAllCandList;
    RhoCandList fChargedCandList;
    RhoCandList fNeutralCandList;
    RhoCandList fMcCandList;

    TString fChargedPidName;
    TString fNeutralPidName;
    TString fTracksName;
    TString fTracksName2;
    ClassDef(PndAnalysis,0);
};


#endif

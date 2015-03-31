#ifndef PndTrkTracking2_H
#define PndTrkTracking2_H 1
#include <vector>


// #include "FairRootManager.h"
#include "FairTask.h"
#include "PndGeoSttPar.h"
#include "PndMCTrack.h"
#include "PndTrkPrintouts.h"
#include "PndSttTrack.h"
#include "PndSttTrackFinder.h"
#include "PndSttTube.h"
#include "PndTrkVectors.h"

#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TList.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

static const Double_t
THETAMIN		= 0.,
  THETAMAX		= 2.*3.141592654;


class FairMCPoint;

class PndTrkTracking2 : public FairTask
{

 public:

  /** Default constructor **/
  PndTrkTracking2();

  /** Second constructor **/
  PndTrkTracking2(Int_t verbose);

  /** Third constructor **/
  PndTrkTracking2(int istamp, bool  iplot, bool doMcComparison);

  /** Fourth constructor **/
  PndTrkTracking2(int istamp, bool  iplot, bool doMcComparison, bool doSciTil);

  /** Destructor **/
  ~PndTrkTracking2();


  void Cleanup( ){fYesCleanStt=true; fYesCleanMvd=true; return;};
  void NOCleanup( ){fYesCleanStt=false;  fYesCleanMvd=false;   return;};

  void CleanupStt( ){fYesCleanStt=true; return;};
  void NOCleanupStt( ){fYesCleanStt=false; return;};
  void CleanupMvd( ){fYesCleanMvd=true; return;};
  void NOCleanupMvd( ){fYesCleanMvd=false;  return;};


  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt);

  /** Virtual method Init **/
  virtual InitStatus Init();


  void NoMvdAloneTracking( ){ fMvdAloneTracking=false; return;};


  void PrintTime()
  {
  	cout<<"\nMy calculation of the time is :";
  	cout << " real time " << frtime2 << " sec., CPU time " << fctime2
  	<< " seconds." << endl << endl;
  	return;

  };

  void SetInputBranchName(
	char* string1,
	char* string2,
	char* string3
	)
  {
	sprintf(fSttBranch,"%s", string1);
	sprintf(fMvdPixelBranch,"%s", string2);
	sprintf(fMvdStripBranch,"%s", string3);
	return;
  };

  void SetEventsToPlot(int nev ){fNevents_to_plot = nev;};

  void SetParContainers();

  /** set persistence flag **/
  void SetPersistence(Bool_t persistence) { fPersistence = persistence; }

  void YesMvdAloneTracking( ){ fMvdAloneTracking=true; return;};


  //  write out the histograms
  void WriteHistograms();



 private:

// divisions of the Theta angle range [0, 90  degrees );
#define	BOCA_90DEGREES_DIVISIONS	 180

  static const Short_t
	LEGIANDRE_NRADIUSDIV	= 100,
	LEGIANDRE_NTHETADIV	= 2* BOCA_90DEGREES_DIVISIONS ,  // divisions  of the Theta angle range [0, 2*PI radians);
//	LEGIANDRE_NRADIUSDIV	= 100,
//	LEGIANDRE_NTHETADIV	= 360,
	MAXMCTRACKS		= 100,
	MAXMVDPIXELHITS		= 500,
	MAXMVDPIXELHITSINTRACK	= 10,
	MAXMVDSTRIPHITS		= 500,
	MAXMVDSTRIPHITSINTRACK	= 10,
	MAXMVDTRACKSPEREVENT	= 400,
	MAXSCITILHITS		= 200, // max SciTil hits total.
	MAXSCITILHITSINTRACK	= 2,   // max SciTil hits in one track.
	MAXSTTHITS		= 900,
	MAXSTTHITSINTRACK	= 40,
	MAXTRACKSPEREVENT	= 200,
//	NFIDIVCONFORMAL		= (Short_t) (3.141592654 * 45./0.5),
	NFIDIVCONFORMAL		= 282,
	NRDIVCONFORMAL		= 10,
	NUMBER_STRAWS		= 4542; // the straw numbers
				// start at 1 and goes up to 4542 included;

  bool
	doMcComparison,
	fSingleHitListStt[MAXSTTHITS],
	iplotta,
	fMvdAloneTracking,
	fYesCleanMvd,
	fYesCleanStt,
	fYesSciTil,
	fInclusionListStt[MAXSTTHITS],
	fInclusionListSciTil[MAXSCITILHITS],
	finMvdTrackCandPixel[MAXMVDPIXELHITS],
	finMvdTrackCandStrip[MAXMVDSTRIPHITS],
	fTypeConf[MAXTRACKSPEREVENT];

  /** object persistence **/
  Bool_t
		fPersistence;



  /**  Branch names to be used to fetch the hits of the backgound mixed events  **/

  char
	fSttBranch[200],
	fMvdPixelBranch[200],
	fMvdStripBranch[200];




  Short_t
	fnAxialOuterRight,  //  number of axial Stt, outer, on the right (looking into the beam);
	fnAxialInnerRight,  //  number of axial Stt, inner, on the right (looking into the beam);
	fnAxialOuterLeft,  //  number of axial Stt, outer, on the left (looking into the beam);
	fnAxialInnerLeft,  //  number of axial Stt, inner, on the left (looking into the beam);

	fListAxialOuterRight[NUMBER_STRAWS],  //  list of axial Stt, outer, on the right (looking into the beam);
	fListAxialInnerRight[NUMBER_STRAWS],  //  list of axial Stt, inner, on the lright (looking into the beam);
	fListAxialOuterLeft[NUMBER_STRAWS],  //  list of axial Stt, outer, on the left (looking into the beam);
	fListAxialInnerLeft[NUMBER_STRAWS],  //  list of axial Stt, inner, on the left (looking into the beam);

	fnSkewRight,  //  number of skew Stt, on the right (looking into the beam);
	fnSkewLeft,  //  number of skew Stt, on the right (looking into the beam);
	fListSkewRight[NUMBER_STRAWS],  //  list of axial Stt, inner, on the lright (looking into the beam);
	fListSkewLeft[NUMBER_STRAWS],  //  list of axial Stt, outer, on the left (looking into the beam);

	fListHitMvdTrackCand[MAXMVDTRACKSPEREVENT][MAXMVDPIXELHITSINTRACK+MAXMVDSTRIPHITSINTRACK],
	fListHitTypeMvdTrackCand[MAXMVDTRACKSPEREVENT][MAXMVDPIXELHITSINTRACK+MAXMVDSTRIPHITSINTRACK],
	fListMvdDSPixelHitNotTrackCand[MAXMVDPIXELHITS],
	fListMvdDSStripHitNotTrackCand[MAXMVDSTRIPHITS],
	fListMvdPixelHitsinTrack[MAXTRACKSPEREVENT][MAXMVDPIXELHITSINTRACK],
	fListMvdStripHitsinTrack[MAXTRACKSPEREVENT][MAXMVDSTRIPHITSINTRACK],
	fListMvdUSPixelHitNotTrackCand[MAXMVDPIXELHITS],
	fListMvdUSStripHitNotTrackCand[MAXMVDSTRIPHITS],
	fListParContiguous[NUMBER_STRAWS][6],
	fListSciTilHitsinTrack[MAXTRACKSPEREVENT][MAXSCITILHITSINTRACK],
	fListSttParHitsinTrack[MAXTRACKSPEREVENT][MAXSTTHITSINTRACK],
	fListSttParHits[MAXSTTHITS],
	fListSttSkewHitsinTrack[MAXTRACKSPEREVENT][MAXSTTHITSINTRACK],
	fListSttSkewHitsinTrackSolution[MAXTRACKSPEREVENT][MAXSTTHITSINTRACK],
	fListSttSkewHits[MAXSTTHITS],
	fListTrackCandHit[MAXTRACKSPEREVENT]
				[MAXSTTHITSINTRACK+
				MAXMVDPIXELHITSINTRACK+
				MAXMVDSTRIPHITSINTRACK+
				MAXSCITILHITSINTRACK],
  //  type = 0 --> Mvd Pixel;  type = 1 --> Mvd Strip; type = 1 --> Mvd Strip; type = 2 --> Stt Parallel
  //  type = 3 --> Stt Straw; type 1001 --> SciTil;  type -1 -->  noise.
	fListTrackCandHitType[MAXTRACKSPEREVENT]
				[MAXSTTHITSINTRACK+
				MAXMVDPIXELHITSINTRACK+
				MAXMVDSTRIPHITSINTRACK+
				MAXSCITILHITSINTRACK],
	fnHitMvdTrackCand[MAXMVDTRACKSPEREVENT],
	fnMCTracks,
	fMCtrack_of_Pixel[MAXMVDPIXELHITS],
	fMCtrack_of_Strip[MAXMVDSTRIPHITS],
	fnMvdDSPixelHitNotTrackCand,
	fnMvdDSStripHitNotTrackCand,
	fnMvdPixelHit,
	fnMvdPixelHitsinTrack[MAXTRACKSPEREVENT],
	fnMvdStripHit,
	fnMvdStripHitsinTrack[MAXTRACKSPEREVENT],
	fnMvdTrackCand,
	fnMvdUSPixelHitNotTrackCand,
	fnMvdUSStripHitNotTrackCand,
	fnParContiguous[NUMBER_STRAWS],
	fnSciTilHits,
	fnSciTilHitsinTrack[MAXTRACKSPEREVENT],
	fnSttParHitsinTrack[MAXTRACKSPEREVENT],
	fnSttSkewHitsinTrack[MAXTRACKSPEREVENT],
	fnTrackCandHit[MAXTRACKSPEREVENT],
	fStrawCode[NUMBER_STRAWS],
	fStrawCode2[NUMBER_STRAWS],
	fTubeID[MAXSTTHITS];


  int
	fNevents_to_plot,
	istampa,
	IVOLTE ;

  Double_t
	fALFA[MAXTRACKSPEREVENT],
	fBETA[MAXTRACKSPEREVENT],
	fBFIELD,

	fCandidatePixelDriftRadius[MAXMVDPIXELHITS],
	fCandidatePixelErrorDriftRadius[MAXMVDPIXELHITS],
	fCandidatePixelS[MAXMVDPIXELHITS],
	fCandidatePixelZ[MAXMVDPIXELHITS],

	fCandidateStripDriftRadius[MAXMVDSTRIPHITS],
	fCandidateStripErrorDriftRadius[MAXMVDSTRIPHITS],
	fCandidateStripS[MAXMVDSTRIPHITS],
	fCandidateStripZ[MAXMVDSTRIPHITS],


	fCandidateSciTilDriftRadius,
	fCandidateSciTilErrorDriftRadius,
	fCandidateSciTilS,
	fCandidateSciTilZ,

	fCandidateSkewS[2*MAXSTTHITS],	// here 2*MAXSTTHITS because in principle the skew straw may have
	fCandidateSkewZ[2*MAXSTTHITS],  // 2 intersections wit the helix cylinder;
	fCandidateSkewZDrift[2*MAXSTTHITS],  //            "
	fCandidateSkewZError[2*MAXSTTHITS],  //            "

	fCosine[LEGIANDRE_NTHETADIV],
	fCxMC[MAXMCTRACKS],
	fCyMC[MAXMCTRACKS],
	fDELTATHETA,
	fFimin,
	fGAMMA[MAXTRACKSPEREVENT],
	fMCSkewAloneX[MAXSTTHITS],
	fMCSkewAloneY[MAXSTTHITS],
	fMCtruthTrkInfo[15][MAXMCTRACKS],
	fOx[MAXTRACKSPEREVENT],
	fOy[MAXTRACKSPEREVENT],
	fposizSciTil[MAXSCITILHITS][3],
	fpSciTilx[MAXSCITILHITS],
	fpSciTily[MAXSCITILHITS],
	fpSciTilz[MAXSCITILHITS],
	fR[MAXTRACKSPEREVENT],
	frefindexMvdPixel[MAXMVDPIXELHITS],
	frefindexMvdStrip[MAXMVDSTRIPHITS],
	fradiaConf[NRDIVCONFORMAL],
	fR_MC[MAXMCTRACKS],
	fsigmaXMvdPixel[MAXMVDPIXELHITS],
	fsigmaYMvdPixel[MAXMVDPIXELHITS],
	fsigmaZMvdPixel[MAXMVDPIXELHITS],
	fsigmaXMvdStrip[MAXMVDSTRIPHITS],
	fsigmaYMvdStrip[MAXMVDSTRIPHITS],
	fsigmaZMvdStrip[MAXMVDSTRIPHITS],
	fSinus[LEGIANDRE_NTHETADIV],
	fS_SciTilHitsinTrack[MAXTRACKSPEREVENT][MAXSCITILHITS],
	fxTube[NUMBER_STRAWS],
	fxxyyTube[NUMBER_STRAWS],
	fyTube[NUMBER_STRAWS],
	fzTube[NUMBER_STRAWS],
	fXMvdPixel[MAXMVDPIXELHITS],
	fXMvdStrip[MAXMVDSTRIPHITS],
	fYMvdPixel[MAXMVDPIXELHITS],
	fYMvdStrip[MAXMVDSTRIPHITS],
	fZMvdPixel[MAXMVDPIXELHITS],
	fZMvdStrip[MAXMVDSTRIPHITS];

//----------------------- real, cputime stuff;

  Double_t
	fctime,
	fctime2,
	frtime,
	frtime2;
  TStopwatch
	ftimer,
	ftimer2;
//--------------------------------


//  FairRootManager *ioman;

  FILE
	* HANDLE,
	* HANDLE2;

  TH1F
	*hdeltaRPixel,
	*hdeltaRStrip,
	*hdeltaRPixel2,
	*hdeltaRStrip2;


  TClonesArray
  /** MC Track Array  **/
	*fMCTrackArray,
 /** Input array of MC points  of Mvd**/
	*fMvdMCPointArray,
 /** Input array of MvdPixelHitArray **/
	*fMvdPixelHitArray,
 /** Input array of MvdStripHitArray **/
	*fMvdStripHitArray,
 /** Input array of PndTracksCand of Mvd**/
	*fMvdTrackCandArray,
 /** SciTil Hit Array **/
	*fSciTHitArray,
 /** SciTil MC Point Array **/
	*fSciTPointArray,
 /** Input array of PndSttTube (map of STT tubes) **/
	*fSttTubeArray,
 /** Input array of PndSttPoints **/
	*fSttPointArray,
 /** Input array of PndSttHit **/
	*fSttHitArray,
 /** Input array of PndSttTracks **/
	*fSttTrackArray,
 /** Input array of PndTracksCand of Stt **/
	*fSttTrackCandArray,
 /** Output array of PndSttMvd  PndTrackCand **/
	*fSttMvdPndTrackCandArray,
 /** Output array of PndSttMvd   PndTrack **/
	*fSttMvdPndTrackArray;

  PndGeoSttPar *fSttParameters;  //  CHECK added



  void Initialization_ClassVariables();


  bool  AcceptHitsConformal(
	Double_t  distance,
	Double_t  DriftConfR, //drift radius in conformal space
	Double_t  StrawConfR  // straw radius in conformal space
	);



  Short_t AssociateBetterAfterFitSkewHitsToXYTrack(
	Short_t TemporarynSttSkewhitinTrack,  //  input
	Short_t SkewList[][2], // input,  list of selected skew hits (in skew numbering)
	Double_t *S,       //  input,  S coordinate of selected Skew hit
	Double_t *Z,       //  input,  Z coordinate of selected Skew hit
	Double_t *ZDrift,  //  input,  drift distance IN Z DIRECTION only, of selected Skew hit
	Double_t *ZError,   //  input,  error (in SZ space) IN Z DIRECTION only, of selected Skew hit
	Double_t KAPPA,    // input, KAPPA result of fit
	Double_t FI0,    // input, FI0 result of fit
	Short_t *tempore,  //  output result, associated skew hits
	Double_t *temporeS,  //  output, associated skew hit  S
	Double_t *temporeZ,  //  output, associated skew hits Z
	Double_t *temporeZDrift,  //  output, associated skew hit Z drift
	Double_t *temporeZError,  //  output, associated skew hits Z error;
	Short_t  *STATUS   // output
	);




  Short_t AssociateSciTilHit(
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Short_t *List, // output, list of SciTil hits associated (max. 2);
	Double_t *esse // output, list of  S of the SciTil hits associated.
	);


  Short_t AssociateSkewHitsToXYTrack(
	bool *InclusionListSkew,
	Short_t NSkewhits,
	Short_t *infoskew,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Double_t info[][7],
	Double_t *WDX,
	Double_t *WDY,
	Double_t *WDZ,
	Double_t Fi_low_limit,
	Double_t Fi_up_limit,
	Short_t  Charge,
	Short_t SkewList[][2], // output,list of selected skew hits (skew numbering)
	Double_t *S,       //  output,  S coordinate of selected Skew hit
	Double_t *Z,       //  output,  Z coordinate of selected Skew hit
	Double_t *ZDrift,   //  output,  drift distance IN Z DIRECTION only,
			   // of selected Skew hit
	Double_t *ZError   //  output,  error (in SZ space) IN Z DIRECTION only, of selected Skew hit.
	);



  Short_t AssociateSkewHitsToXYTrack2(
	bool *InclusionListSkew,
	Short_t NSkewhits,
	Short_t *infoskew,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Double_t info[][7],
	Double_t *WDX,
	Double_t *WDY,
	Double_t *WDZ,
	Double_t Fi_low_limit,
	Double_t Fi_up_limit,
	Short_t  Charge,
	Short_t SkewList[][2], // output,list of selected skew hits (skew numbering)
	Double_t *S,       //  output,  S coordinate of selected Skew hit
	Double_t *Z,       //  output,  Z coordinate of selected Skew hit
	Double_t *ZDrift,   //  output,  drift distance IN Z DIRECTION only,
			   // of selected Skew hit
	Double_t *ZError   //  output,  error (in SZ space) IN Z DIRECTION only, of selected Skew hit.
	);


  bool BadTrack_ParStt(
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Short_t Charge,
	Double_t Xcross[2],  // Xcross[0]=point of entrance;
	//  Xcross[1]=point of exit.
	Double_t Ycross[2],
	Short_t nHits,
	Short_t* ListHits,
	Double_t info[][7],
	Double_t cut,
	Short_t maxnum,
	Short_t islack // uncertainty allowed as far as
	// the n. of hits that should be present.
	);

 void CalculateSinandCosin(
	);


  void CollectParSttHitsagain(
	Vec <bool>& keepit,
	Vec <bool>& Mvdhits,
	Double_t info[][7],
	Short_t nSttParHit,
	Short_t StartTrackCand,
	Short_t EndTrackCand,
	Double_t *KAPPA,
	Double_t *FI0,
	Double_t *Fi_low_limit,
	Double_t *Fi_up_limit,
	Short_t *fnSttParHitsinTrack, // input/output
	Short_t fListSttParHitsinTrack[][MAXSTTHITSINTRACK] // input/output
	);

  Short_t CompareTracks(
	Short_t first_track,
	Short_t second_track
			);


  bool EliminateClones(
	Short_t nTotalCandidates, // input;
	Double_t fraction,// input; raction of common hits to declare the two tracks clones;
	bool * keepit  // input and output;
			);


  void EliminateSpuriousSZ_bis(
  	Short_t ncand,
	Short_t MaxTurnofTracks,
	Double_t signPz,
	Double_t *SchosenPixel,
	Double_t *SchosenStrip,
	Double_t *SchosenSkew,
	Double_t *ZchosenPixel,
	Double_t *ZchosenStrip,
	Double_t *ZchosenSkew,
	Double_t *ErrorchosenPixel,
	Double_t *ErrorchosenStrip,
	Double_t *ErrorchosenSkew,
	Double_t KAPPA,
	Double_t FI0,
	Double_t Rr
	);



  void EliminateSpuriousSZ_ter(
  	Short_t ncand,
	Short_t MaxTurnofTracks,
	Double_t signPz,
	Double_t *SchosenPixel,
	Double_t *SchosenStrip,
	Double_t *SchosenSkew,
	Double_t *ZchosenPixel,
	Double_t *ZchosenStrip,
	Double_t *ZchosenSkew,
	Double_t *ErrorchosenPixel,
	Double_t *ErrorchosenStrip,
	Double_t *ErrorchosenSkew,
	Double_t KAPPA,
	Double_t FI0,
	Double_t Rr
	);
 void ExtractInfoFromMvdTrackCand();




  void FindCharge(
	Double_t oX,
	Double_t oY,
	Short_t nHits,
	Double_t *X,
	Double_t *Y,
	Short_t  * Charge
	);




  Short_t FindTrackStrictCollection(
	Short_t NFiCELLDISTANCE,
	//  seed track (original notation) as far as the Fi angle is concerned
	Short_t iSeed,
	//  n. of hits to search in ListHitsinTrackinWhichToSearch
	Short_t NParallelToSearch,
	Short_t *ListHitsinTrackinWhichToSearch,
	bool *InclusionList,
	Short_t *FiConformalIndex,
	Short_t  *OutputListHitsinTrack
	);


  void FixDiscontinuitiesFiangleinSZplane(
	Short_t TemporarynSkewHitsinTrack,
	Vec <Double_t> & S,
	Double_t *Fi_initial_helix_referenceframe,
	Short_t Charge
	);


  void InfoXYZParal(
	Double_t info[][7],
	Short_t infopar,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Double_t KAPPA,
	Double_t FI0,
	Short_t Charge,
	Double_t *Posiz	//  output
	);


  void  Initial_SttParHits_DecreasingR_Ordering(
	Double_t info[][7],
	Short_t *ListSttParHi,
	Int_t nSttParHit
	);


  void LoadPndTrack_TrackCand(
	bool* keepit,
	bool* SttSZfit,
	Short_t nTotalCandidates,
	Short_t* Charge,
	Int_t nSttTrackCand,
	Double_t *FI0,
	Double_t *KAPPA,
	Double_t info[][7],
	Double_t SchosenSkew[][MAXSTTHITS],
	Double_t ZchosenSkew[][MAXSTTHITS],
	Short_t *daTrackFoundaTrackMC
	);


  void LoadSZetc_forSZfit(
	Short_t ncand,	// input
	Short_t nhitsinfit,
//	Double_t * TemporaryS,		// input
//	Double_t * TemporaryZ,		// input
//	Double_t * TemporaryZDrift,	// input
//	Double_t * TemporaryZError,	// input

	Vec <Double_t>& ErrorDriftRadius,	 // output
	Double_t * ErrorDriftRadiusbis,	 // output
	Vec <Double_t>& DriftRadius,		 // output
	Double_t * DriftRadiusbis,	 // output
	Vec <Double_t> & S,			 // output
	Double_t * Sbis,		 // output
	Vec <Double_t>& ZED,			 // output
	Double_t * ZEDbis		 // output
	);



  void LoadSZetc_forSZfit2(
	Short_t ncand,	// input
	Short_t nhitsinfit,
	Double_t * S_Skew,		// input
	Double_t * TemporaryZ,		// input
	Double_t * TemporaryZDrift,	// input
	Double_t * TemporaryZError,	// input

	Vec <Double_t>& ErrorDriftRadius,	 // output
	Double_t * ErrorDriftRadiusbis,	 // output
	Vec <Double_t>& DriftRadius,		 // output
	Double_t * DriftRadiusbis,	 // output
	Vec <Double_t> & S,			 // output
	Double_t * Sbis,		 // output
	Vec <Double_t>& ZED,			 // output
	Double_t * ZEDbis		 // output
	);

  void MakeInclusionListStt(
	Int_t nSttHit,
	Short_t * TubeID,
	Double_t info[][7]
	);



  void MatchMvdHitsToSttTracks(
	Vec<bool> &keepit,
	Double_t delta,
	Double_t highqualitycut,
	Short_t nSttTrackCand,
	Double_t *FI0,
	Double_t *Fifirst,
	Vec <Short_t>& CHARGE,
	Short_t *nPixelHitsinTrack, // output
	Short_t ListPixelHitsinTrack[][MAXMVDPIXELHITSINTRACK], // output
	Short_t *nStripHitsinTrack, // output
	Short_t ListStripHitsinTrack[][MAXMVDSTRIPHITSINTRACK] // output
	);




  void MatchMvdHitsToSttTracksagain(
	Vec <bool>& keepit,
	Vec <bool>& Mvdhits,
	Double_t delta,
	Double_t highqualitycut,
	Short_t nSttTrackCand,
	Double_t *FI0,
	Double_t *Fifirst,
	Vec <Short_t>& CHARGE,

	Short_t *nPixelHitsinTrack, // output
	Short_t ListPixelHitsinTrack[][MAXMVDPIXELHITSINTRACK], // output
	Short_t *nStripHitsinTrack, // output
	Short_t ListStripHitsinTrack[][MAXMVDSTRIPHITSINTRACK] // output
	);




  void MatchMvdHitsToSttTracks2(
	bool *keepit,
	Double_t delta,
	Double_t highqualitycut,
	Short_t nSttTrackCand,
	Double_t *FI0,
	Double_t *Fifirst,
	Short_t *CHARGE,
	Short_t *nPixelHitsinTrack, // output
	Short_t ListPixelHitsinTrack[][MAXMVDPIXELHITSINTRACK], // output
	Short_t *nStripHitsinTrack, // output
	Short_t ListStripHitsinTrack[][MAXMVDSTRIPHITSINTRACK] // output
	);



  void OrderingConformal_Loading_ListTrackCandHit(
	bool* keepit,
	Short_t ncand,
	Double_t info[][7],
	Double_t Trajectory_Start[][2],
	Short_t* CHARGE,
	Double_t SchosenSkew[][MAXSTTHITS]
	);

  void OrderingR_Loading_ListTrackCandHit(
	bool* keepit,
	Short_t ncand,
	Double_t info[][7]
	);

  void OrderingSttSkewandSttParallel(
	Double_t oX,
	Double_t oY,
	Double_t Rr,
	Short_t nSkewhit,
	Short_t *ListSkewHits,
	Double_t *SList, // it is related to the skew hits. IMPORTANT :
		// the index must be the ORIGINAL skew hit number,
		// therefore SList[ListSkewHits[*]].
	Short_t  Charge,
	Short_t nParHits,
	Short_t *ListParHits,
	Double_t *U,
	Double_t *V,
	Short_t *BigList // final ordered Parallel+Skew list;
		// already in NATIVE hit number.
	);

  void OrderingUsingConformal(
	Double_t oX,
	Double_t oY,
	Double_t Traj_Sta[2],
	Int_t nHits,
	Double_t XY[][2], // XY[*][0] = X position, XY[*][0] = Y position.
	Short_t  Charge,  // input
	Int_t *ListHits
	);




  void Ordering_Loading_ListTrackCandHit(
//	Vec <bool>& keepit,
	bool * keepit,
	Short_t FirstCandidate,
	Short_t LastCandidate,
	Double_t info[][7],
	Double_t Trajectory_Start[][2],
//	Vec <Short_t>& CHARGE,
	Short_t * CHARGE,
	Double_t SchosenSkew[][MAXSTTHITS]
	);
  void RefitMvdStt(
	Short_t nCandHit,
	Short_t *fListTrackCandHit,
	Short_t *fListTrackCandHitType,
	Double_t info[][7],
	Double_t rotationangle,
	Double_t trajectory_vertex[2],
	Short_t iexcl,
	Double_t *pAlfa, // output of the fit
	Double_t *pBeta, // output of the fit
	Double_t *pGamma,// set at zero always for now
	bool *status    // fit status; true = successful
	  );

  void StartFromSciTil(
	Short_t * Charge,
	Short_t * FiConformalIndex,
	Double_t *Fi_final_helix_referenceframe,
	Double_t *Fi_initial_helix_referenceframe,
	Double_t * Fi_low_limit,
	Double_t * Fi_up_limit,
	Short_t HitsinBoxConformal[][NRDIVCONFORMAL][NFIDIVCONFORMAL],
	Double_t info[][7],
	Double_t infoparalConformal[][5],
	Short_t nBoxConformal[][NFIDIVCONFORMAL],
	Int_t nSttParHit,
	Int_t &nSttTrackCand,
	Short_t * RConformalIndex,
	Double_t *trajectory_vertex,
	Double_t *UU,
	Double_t *VV
	);

  void SeparateInnerOuterParallel(

	// input
	Short_t nHits,
	Short_t *ListHits,
	Double_t info[][7],
	Double_t RStrawDetInnerParMax,

	// output
	Short_t *nInnerHits,
	Short_t *ListInnerHits,
	Short_t *nOuterHits,
	Short_t *ListOuterHits,

	Short_t *nInnerHitsLeft,
	Short_t *ListInnerHitsLeft,
	Short_t *nInnerHitsRight,
	Short_t *ListInnerHitsRight,

	Short_t *nOuterHitsLeft,
	Short_t *ListOuterHitsLeft,
	Short_t *nOuterHitsRight,
	Short_t *ListOuterHitsRight
	);

  bool SttParalCleanup(
	Double_t GAP,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Short_t Charge,
	Double_t Start[3],
	Double_t FI0,
	Double_t FiLimitAdmissible,
	Short_t nHits,
	Short_t *Listofhits,
	Double_t info[][7],
	Double_t RStrawDetMin,
	Double_t RStrawDetInnerParMax,
	Double_t RStrawDetOuterParMin,
	Double_t RStrawDetMax
	);


  bool SttSkewCleanup(
	Double_t GAP,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Short_t  Charge,
	Double_t Start[3],
	Double_t FI0,
	Double_t FiLimitAdmissible,
	Short_t nHits,
	Short_t *Listofhits,
	Double_t *S,
	Double_t info[][7],
	Double_t RminStrawSkew,
	Double_t RmaxStrawSkew,
	Double_t cut,
	Short_t maxnum
	);


  bool TrackCleanup(
	Double_t GAP,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Double_t KAPPA,
	Double_t FI0,
	Short_t  Charge,
	Double_t Start[3],
	Short_t &nHitsPar,
	Short_t *ListHitsPar,
	Short_t &nHitsSkew,
	Short_t *ListHitsSkew,
	Double_t *auxS,
	Double_t info[][7],
	Double_t RStrawDetMin,
	Double_t ApotemaMaxInnerPar,
	Double_t ApotemaMinSkew,
	Double_t ApotemaMaxSkew,
	Double_t ApotemaMinOuterPar,
	Double_t RStrawDetMax
	);


  Short_t TrkAssociatedParallelHitsToHelixQuater(
	bool *ExclusionList,
	Double_t m,
	Double_t q,
	Short_t Status,
	Short_t nHitsinTrack,
	Short_t *ListHitsinTrack,
	Int_t NhitsParallel,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Double_t info[][7],
	Double_t infoparalConformal[][5],
	Short_t *RConformalIndex,
	Short_t *FiConformalIndex,
	Short_t nBoxConformal[][NFIDIVCONFORMAL],
	Short_t HitsinBoxConformal[][NRDIVCONFORMAL][NFIDIVCONFORMAL],
	Short_t *auxListHitsinTrack
	);



  Short_t TrkAssociatedParallelHitsToHelix5(
	bool *ExclusionList,
	Int_t NhitsParallel,
	Double_t Oxx,
	Double_t Oyy,
	Double_t Rr,
	Double_t info[][7],
	Double_t Fi_low,
	Double_t Fi_up,
	Short_t *auxListHitsinTrack
	);


  void StoreSZ_MvdScitil(Short_t ncand);


  ClassDef(PndTrkTracking2,1);

};

#endif

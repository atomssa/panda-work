/**
 @class PndFtsHoughTrackFinder

 @author Martin J. Galuska <martin [dot] j [dot] galuska (at) physik [dot] uni (minus) giessen [dot] de>

 @brief Implementation of the Hough transform based FTS PR. Creates Hough spaces, finds peaks (=tracklets) and combines them to track candidates.

 This is a class version of the HoughTest.C macro PR test implementation
 minus all the plotting stuff.
 Take a look at the notes of the macro version for further details.

 Recent Changes
 Major code cleanup and deletion of test code / unneeded code
 use PndFtsHoughTrackCand to store information about track candidates and Hough transforms
 Find all peaks with a minimum height, analysing peak shapes.
 Moved all Hough space related code to PndFtsHoughSpace -> Major code simplification, better maintainability
 Fill PndTrackCands and PndTrack for output

 This class is loosely modeled after the
 mvd/MvdTracking/PndRiemannTrackFinder
 sttmvdtracking/PndMvdSttGemRiemannTrackFinder
 classes


 TODO
 Match straight line for stations 5+6 to parabola
 Add skewed hits
 Add drift circles
 Adaptive Hough

 Created: 18.06.2013
 */




#ifndef PndFtsHoughTrackFinder_H
#define PndFtsHoughTrackFinder_H

//#include "TClonesArray.h"
#include "Rtypes.h" // for Double_t, Int_t, etc
#include "FairLogger.h" // for FairLogger, MESSAGE_ORIGIN

#include <cmath>
#include "TMath.h"
#include <math.h>
#include <algorithm>
#include <set>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

#include "PndFtsHit.h"
#include "PndFtsHoughTrackerTask.h"
#include "PndFtsHoughSpace.h"
#include "PndFtsHoughTracklet.h"
#include "PndFtsHoughTrackCand.h"
#include "PndTrackCand.h"
#include "PndTrack.h"

// For error throwing
#include "TString.h"
#include <stdexcept>

static const Double_t meinpi = 3.14159265359;
///< sets where the apex of the parabola is supposed to be
static const Double_t fZLineParabola = 368.; // the value should coincide with the start of the dipole field // 368. was ok
static const Double_t fZParabolaLine = 605.; // the value should coincide with the end of the dipole field // TODO determine this value

class PndFtsHoughTrackFinder
{
public:
	PndFtsHoughTrackFinder(PndFtsHoughTrackerTask *trackerTask); ///< @brief Set pointer to tracker task (super important as it provides functionality such as the array of all FTS hits and the branchId of FTS, magnetic field, etc.)
	virtual ~PndFtsHoughTrackFinder(); ///< @brief Destructor

	void FindTracks();	///< @brief Performs the track finding.


	// Output
	Int_t NTracks() const { return fHoughTrackCands.size(); };		///< @brief Returns the number of found tracks
	/**@brief Returns the track cand. with index i.
	 *
	 * Note: Method calculates first and last parameters of the PndTrack object, but uses an empty PndTrackCand which has to be set lateron using SetTrackCandRef!
	 * @param i Index of requested track cand.
	 * @return Track cand. to index i.
	 */
	PndTrack GetPndTrack(int i){ return fHoughTrackCands[i].getPndTrack(); };
	/**@brief Returns the track cand. with index i.
	 *
	 * Note: Use this to add a PndTrackCand to the corresponding PndTrack object with SetTrackCandRef!
	 * @param i Index of requested track cand.
	 * @return Track cand. to index i.
	 */
	PndTrackCand GetPndTrackCand(int i) { return fHoughTrackCands[i].getPndTrackCand(); };
	/**@brief Returns the track cand. with index i.
	 *
	 * Note: For debugging only.
	 * @param i Index of requested track cand.
	 * @return Track cand. to index i.
	 */
	PndFtsHoughTrackCand GetHoughTrack(int i) const { return fHoughTrackCands[i]; };


	// Parameters
	//	void SetMinPeakHeightZxLineParabola(UInt_t val){ fMinPeakHeightZxLineParabola = val; };
	//	void SetMinPeakHeightZxParabola(UInt_t val){ fMinPeakHeightZxParabola = val; };
	//	void SetMinPeakHeightZxParabolaLine(UInt_t val){ fMinPeakHeightZxParabolaLine = val; };
	//	void SetMinPeakHeightZyLine(UInt_t val){ fMinPeakHeightZyLine= val; };




private:
	/// @brief Task which handles PandaRoot input/output and provides settings.
	/// Has to be set using the constructor.
	PndFtsHoughTrackerTask *fTrackerTask;

	/** @brief For error reporting */
	void throwError(const TString s){ throw std::runtime_error(s.Data()); };

	inline void PrintFoundTracklets(const std::vector<PndFtsHoughTracklet>& tracklets, const TString& option) const{
		std::cout << tracklets.size() << " peaks found for " << option << '\n';
		if (10<fVerbose){
			for (UInt_t i=0; i< tracklets.size(); ++i)
			{
				tracklets[i].Print();
			}
		}
	};

	//	Int_t   fFtsBranchId; // needed for saving and accessing hits
	//	TClonesArray *fFtsHitArray; ///< @brief Input array of all FTS hits.



	///< @brief Minimum required height for peaks in Hough spaces.
	UInt_t    fMinPeakHeightZxLineParabola;					///< zx line before dipole field
	///< @brief Minimum required height for peaks in Hough spaces.
	UInt_t    fMinPeakHeightZxParabola;					///< zx parabola within dipole field
	///< @brief Minimum required height for peaks in Hough spaces.
	UInt_t    fMinPeakHeightZxParabolaLine;					///< zx line after dipole field
	///< @brief Minimum required height for peaks in Hough spaces.
	UInt_t    fMinPeakHeightZyLine;					///< zy line


	Int_t fVerbose;
	Bool_t fSaveDebugInfo;



	//
	//	// for B field access
	//	FairField* fField;





	// for Hough
	//-----------
	PndFtsHoughSpace* fHoughSpaceZxLineBeforeDipole;
	PndFtsHoughSpace* fHoughspaceZxParabola;
	PndFtsHoughSpace* fHoughSpaceZxLineBehindDipole;
	PndFtsHoughSpace* fHoughspaceZyLine;
	std::vector<PndFtsHoughTrackCand> fHoughTrackCandsNew;	///< For temporary internal storing of track cands.
	std::vector<PndFtsHoughTrackCand> fHoughTrackCands;		///< For internal storing of track cands.
	//	std::vector<PndTrackCand> fTrackCand; // resulting track candidates, also used for returning PndTracks


	//static const Double_t meinpi = 3.14159265359;
	/////< sets where the apex of the parabola is supposed to be
	//static const Double_t fZLineParabola = 368.; // the value should coincide with the start of the dipole field // 368. was ok
	//static const Double_t fZParabolaLine = 605.; // the value should coincide with the end of the dipole field // TODO determine this value







	// takes the heighest peak (according to peak finder)
	// of all peaks that share > maxSameHits
	/**@brief Filters a vector of tracklets based on the number of shared hits.
	 *
	 * It will only keep the heighest peaks if 2 or more peaks share more than maxAcceptableSharedHits hits. If two peaks have the same height, both are kept.
	 * @param maxAcceptableSharedHits Defines how many hits two tracklets / peaks are allowed to share.
	 * @param[in,out] tracklets Vector containing tracklets which are supposed to be filtered. The vector will be modified.
	 * @return
	 */
	Bool_t FilterTrackletsBasedOnSharedHits(
			UInt_t maxAcceptableSharedHits,
			std::vector<PndFtsHoughTracklet> &tracklets
	);



	ClassDef(PndFtsHoughTrackFinder,1);
};

#endif /*PndFtsHoughTrackFinder_H*/

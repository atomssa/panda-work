#ifndef SIMCOMPLETE_H
#define SIMCOMPLETE_H


#include <TString.h>

void SimComplete(Int_t nEvents, TString const &simEngine, Double_t momentum,
		 Bool_t useEvtGen, Bool_t useDpm, Bool_t useBoxGenerator,
		 Double_t beamMomentum, TString const &outFile, TString const &outParamsFile,
		 TString const &inDigiParamsFile, TString const &trackDetector);


#endif

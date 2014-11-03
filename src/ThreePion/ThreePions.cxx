// The header file
#include "ThreePions.h"

// C++ headers
#include <vector>
#include <string>
#include <iostream>

// FAIR headers
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairRadLenPoint.h"

// ROOT headers
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"

// other headers
//#include "PndEmcBump.h"
//#include "PndEmcCluster.h"
//#include "PndEmcBump.h"
//#include "PndEmcRecoHit.h"
//#include "PndEmcMapper.h"
#include "PndMCTrack.h"

using std::cout;
using std::endl;
using std::vector;

ThreePions::ThreePions() :
  FairTask("Radiation Length Profiler"),nevt(0) {
}

ThreePions::~ThreePions() { }

InitStatus ThreePions::Init() {

  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- ThreePions::Init: " << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  track_array = dynamic_cast<TClonesArray *> (ioman->GetObject("MCTrack"));
  if ( ! track_array ) {
    cout << "-W- ThreePions::Init: " << "No PndMCTrack array!" << endl;
    return kERROR;
  }

  return kSUCCESS;

}

void ThreePions::Exec(Option_t* opt)
{
  if (nevt%1000==0)
    cout << "===== ThreePions::Exec -- Event " << nevt << " ====="<< endl;
}

void ThreePions::Finish()
{

}

ClassImp(ThreePions)

// The header file
#include "RadLenData.h"

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

RadLenData::RadLenData() :
  FairTask("Radiation Length Profiler"),nevt(0) { 
}

RadLenData::~RadLenData() { }

InitStatus RadLenData::Init() 
{		

  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ){
    cout << "-E- RadLenData::Init: " << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  track_array = dynamic_cast<TClonesArray *> (ioman->GetObject("MCTrack"));
  if ( ! track_array ) {
    cout << "-W- RadLenData::Init: " << "No PndMCTrack array!" << endl;
    return kERROR;
  }

  rlp_array = dynamic_cast<TClonesArray *> (ioman->GetObject("RadLen"));
  if ( ! rlp_array ) {
    cout << "-W- RadLenData::Init: " << "No FairRadLenPoint array!" << endl;
    return kERROR;
  }

  for (int idet=0; idet<ndet; ++idet) {
    vProf.push_back(new TProfile(Form("RadLenProf_Det%d",idet),Form("RadLenProf_Det%d",idet),200,0.0,180. ) ) ;
    vProf2d.push_back(new TProfile2D(Form("RadLenProf_ThePhi_Det%d",idet),Form("RadLenProf_ThePhi_Det%d",idet),200,0.0,180.,200,0.0,360.) ) ;
  }
      
  return kSUCCESS;

}


void RadLenData::Exec(Option_t* opt)
{

  if (nevt%1000==0)
    cout << "===== RadLenData::Exec -- Event " << nevt << " ====="<< endl;

  //cout << "=========================================================================" << endl;
  //cout << "================       new event                      ===================" << endl;  
  //cout << "=========================================================================" << endl;
  
  if ( ! rlp_array ) Fatal("Exec", "No Track Array");

  const int nrlp = rlp_array->GetEntriesFast();

  // Assumes only one Geantino per event
  PndMCTrack *track = (PndMCTrack*) track_array->At(0);
  const double theta_deg = 180 * track->GetMomentum().Theta() / TMath::Pi();
  const double phi_deg = 180 * track->GetMomentum().Phi() / TMath::Pi();
  //cout << "Track: phi= " << phi_deg << " theta= " << theta_deg << endl;

  
  vector<double> eff_len_tot(ndet,0.0);

  for (int irlp=0; irlp<nrlp; ++irlp) {

    FairRadLenPoint *rlp = (FairRadLenPoint*) rlp_array->At(irlp);
    
    const int det_id = rlp->GetDetectorID();
    const double rad_len = rlp->GetRadLength();

    TVector3 pos_in = rlp->GetPosition();
    TVector3 pos_out = rlp->GetPositionOut();
    const TVector3 dist_vect = pos_out-pos_in;
    const double eff_len = dist_vect.Mag()/rad_len;
    eff_len_tot[det_id] += eff_len;
    
    //if ( det_id == 6 || det_id == 7 ) { 
    //  //cout << "irlp= " << irlp
    //	// << " ievt= " << rlp->GetEventID()
    //	//   << " itrk= " << rlp->GetTrackID()
    //  cout   << " idet= " << rlp->GetDetectorID()
    //	//<< " dist= " << dist_vect.Mag()
    //	     << " rin= " << pos_in.Mag()
    //	     << " rout= " << pos_out.Mag()
    //	     << " rout-rin= " << pos_out.Mag()-pos_in.Mag()
    //	     << " eff_len= " << eff_len
    //	     << " eff_len_tot= " << eff_len_tot[det_id]
    //	//<< " eloss= " << rlp->GetEnergyLoss()
    //	     << " A= " <<	rlp->GetA()
    //	     << " Zm= " << rlp->GetZm()
    //	     << " RadLen= " << rlp->GetRadLength()
    //	     << " den= " << rlp->GetDensity() << endl;
    //
    //
    //  //cout << "     pos_in= ";      pos_in.Print();
    //  //cout << "     pos_out= ";      pos_out.Print();
    //  //cout << "     dist_vect= ";     dist_vect.Print();
    //
    //}

  }

  // for each detector, add the effective length if it's different from zero?
  for (int idet=0; idet<ndet; ++idet) {
    if (eff_len_tot[idet]>0) {
      vProf[idet]->Fill(theta_deg, eff_len_tot[idet]);
      vProf2d[idet]->Fill(theta_deg, phi_deg, eff_len_tot[idet]);
    }
  }

  ++nevt;

}

void RadLenData::Finish()
{	
  for (int idet=0; idet<ndet; ++idet) {
    vHistCumul.push_back(vProf[idet]->ProjectionX(Form("RadLenProf_CumulUpToDet%d",idet)));
    vHistCumul2d.push_back(vProf2d[idet]->ProjectionXY(Form("RadLenProf_ThePhi_CumulUpToDet%d",idet)));    
    //new TProfile(,Form("RadLenProf_CumulDet%d",idet),200,0.0, 180. ) ) ;
    //vHistCumul[idet]->Add(vProf[idet]);
    if (idet>0) {
      vHistCumul[idet]->Add(vHistCumul[idet-1]);
      vHistCumul2d[idet]->Add(vHistCumul2d[idet-1]);
    }
  }

  for (int idet=0; idet<ndet; ++idet) {
    vProf[idet]->Write();
    vHistCumul[idet]->Write();
    vProf2d[idet]->Write();
    vHistCumul2d[idet]->Write();    
  }
}

ClassImp(RadLenData)

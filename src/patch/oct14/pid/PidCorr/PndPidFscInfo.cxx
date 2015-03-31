#include "PndDetectorList.h"
#include "PndPidCorrelator.h"
#include "PndPidCandidate.h"
#include "PndMCTrack.h"
#include "PndTrack.h"
#include "PndTrackID.h"

#include "PndEmcBump.h"
#include "PndEmcDigi.h"
#include "PndEmcStructure.h"
#include "PndEmcXtal.h"
#include "PndEmcErrorMatrix.h"
#include "PndEmcClusterCalibrator.h"
#include "PndEmcClusterEnergySums.h"

#include "FairTrackParH.h"
#include "FairMCApplication.h"
#include "FairRunAna.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"

Bool_t PndPidCorrelator::GetFscInfo(FairTrackParH* helix, PndPidCandidate* pidCand)
{
  if(!helix)
  {
	  std::cerr << "<Error> PndPidCorrelator FSCINFO: FairTrackParH NULL pointer parameter."<<std::endl;
	  return kFALSE;
  }
  if(!pidCand)
  {
	  std::cerr << "<Error> PndPidCorrelator FSCINFO: pidCand NULL pointer parameter."<<std::endl;
	  return kFALSE;
  }
  //if(helix->GetZ() < fCorrPar->GetZLastPlane())
  if(helix->GetZ() < 747.5)
    {
      if (fVerbose>0) std::cout << "-W- PndPidCorrelator::GetFscInfo :: Skipping forward tracks propagation in forward direction" << std::endl;
      return kFALSE;
    }
  FairGeanePro *fProFsc = new FairGeanePro();
  if (!fCorrErrorProp) fProFsc->PropagateOnlyParameters();
  //---
  Float_t trackTheta = helix->GetMomentum().Theta()*TMath::RadToDeg();

  Int_t fscEntries = fEmcCluster->GetEntriesFast();
  Int_t emcIndex = -1, emcModuleCorr = -1, emcNCrystals = -1;
  Float_t emcEloss = 0., emcElossCorr = 0., emcGLength = -1000;
  Float_t emcQuality = 1000000;
  Float_t chi2 = 0;
  TVector3 vertex(0., 0., 0.); TVector3 emcPos(0., 0., 0.);// TVector3 momentum(0., 0., 0.);

  // Cluster zenike moment
  Double_t Z20 = 0.0, Z53 = 0.0, secLatM = 0.00, E1 = 0., E9 = 0., E25 = 0.;

  for (Int_t ee = 0; ee<fscEntries; ee++)
  {
      PndEmcCluster *fscHit = (PndEmcCluster*)fEmcCluster->At(ee);
      Int_t emcModule = fscHit->GetModule();
      if (emcModule!=5) continue;

      if (fIdeal)
      {
    	  std::vector<Int_t> mclist = fscHit->GetMcList();
    	  if (mclist.size()==0) continue;
    	  if (mclist[0]!=pidCand->GetMcIndex()) continue;
      }

      if (fscHit->energy() < fCorrPar->GetEmc12Thr()) continue;

      emcPos = fscHit->where();
      if (fGeanePro)
      { // Overwrites vertex if Geane is used
		  fProFsc->SetPoint(emcPos);
		  fProFsc->PropagateToPCA(1, 1);
		  vertex.SetXYZ(-10000, -10000, -10000); // reset vertex
		  FairTrackParH *fRes= new FairTrackParH();
		  Bool_t rc =  fProFsc->Propagate(helix, fRes, fPidHyp*pidCand->GetCharge()); // First propagation at module
		  if (!rc) continue;

		  emcGLength = fProFsc->GetLengthAtPCA();
		  vertex.SetXYZ(fRes->GetX(), fRes->GetY(), fRes->GetZ());
		  //std::map<PndEmcTwoCoordIndex*, PndEmcXtal*> tciXtalMap=PndEmcStructure::Instance()->GetTciXtalMap();
		  //PndEmcDigi *lDigi= (PndEmcDigi*)emcHit->Maxima();
		  //PndEmcXtal* xtal = tciXtalMap[lDigi->GetTCI()];
		  //emcPos = xtal->frontCentre();
      }

      Float_t dist = (emcPos-vertex).Mag2();
      if ( emcQuality > dist )
      {
		  emcIndex = ee;
		  emcQuality = dist;
		  emcEloss = fscHit->energy();
		  emcElossCorr = fEmcCalibrator->Energy(fscHit);
		  emcModuleCorr = emcModule;
		  emcNCrystals = fscHit->NumberOfDigis();
		  Z20 = fscHit->Z20();// Z_{n = 2}^{m = 0}
		  Z53 = fscHit->Z53();// Z_{n = 5}^{m = 3}
		  secLatM = fscHit->LatMom();
		  if (fEmcDigi)
		  {
			  PndEmcClusterEnergySums esum(*fscHit, fEmcDigi);
			  E1  = esum.E1();
			  E9  = esum.E9();
			  E25 = esum.E25();
		  }
      }

      if ((fClusterQ[ee]<0) || (dist < fClusterQ[ee]))
	// If the track-fsc distance is less than the previous stored value (or still not initialized)
      {
    	  fClusterQ[ee] = dist; // update the param
      }

      if (fDebugMode)
      {
	Float_t ntuple[] = {static_cast<Float_t>(vertex.X()), static_cast<Float_t>(vertex.Y()), static_cast<Float_t>(vertex.Z()), static_cast<Float_t>(vertex.Phi()),
			    static_cast<Float_t>(helix->GetMomentum().Mag()), static_cast<Float_t>(helix->GetQ()), static_cast<Float_t>(helix->GetMomentum().Theta()), static_cast<Float_t>(helix->GetZ()),
			    static_cast<Float_t>(emcPos.X()), static_cast<Float_t>(emcPos.Y()), static_cast<Float_t>(emcPos.Z()), static_cast<Float_t>(emcPos.Phi()),
			    dist, static_cast<Float_t>(vertex.DeltaPhi(emcPos)), static_cast<Float_t>(fscHit->energy()), emcGLength, static_cast<Float_t>(emcModule)};
    	  // Float_t ntuple[] = {vertex.X(), vertex.Y(), vertex.Z(), vertex.Phi(),
	  // 		    helix->GetMomentum().Mag(), helix->GetQ(), helix->GetMomentum().Theta(), helix->GetZ(),
	  // 		    emcPos.X(), emcPos.Y(), emcPos.Z(), emcPos.Phi(),
	  // 		    dist, vertex.DeltaPhi(emcPos), fscHit->energy(), emcGLength, emcModule};
    	  fscCorr->Fill(ntuple);
      }
  }// End for(ee = 0;)

  if ((emcQuality < fCorrPar->GetEmc12Cut()) || (fIdeal && emcIndex!=-1))
  {
	  fClusterList[emcIndex] = kTRUE;
	  pidCand->SetEmcQuality(emcQuality);
	  pidCand->SetEmcRawEnergy(emcEloss);
	  pidCand->SetEmcCalEnergy(emcElossCorr);
	  pidCand->SetEmcIndex(emcIndex);
	  pidCand->SetEmcModule(emcModuleCorr);
	  pidCand->SetEmcNumberOfCrystals(emcNCrystals);
	//=======
	  pidCand->SetEmcClusterZ20(Z20);
	  pidCand->SetEmcClusterZ53(Z53);
	  pidCand->SetEmcClusterLat(secLatM);
	  pidCand->SetEmcClusterE1(E1);
	  pidCand->SetEmcClusterE9(E9);
	  pidCand->SetEmcClusterE25(E25);
		//=====
  }

  return kTRUE;
}

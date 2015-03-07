#include "PndPidCorrelator.h"
#include "PndEmcClusterEnergySums.h"
#include <cmath>

//_________________________________________________________________
Bool_t PndPidCorrelator::GetEmcInfo(FairTrackParH* helix, PndPidCandidate* pidCand) {
  if(! helix){
    std::cerr << "<Error> PndPidCorrelator EMCINFO: FairTrackParH NULL pointer parameter."
	      <<std::endl;
    return kFALSE;
  }
  if(! pidCand){
    std::cerr << "<Error> PndPidCorrelator EMCINFO: pidCand NULL pointer parameter."
	      <<std::endl;
    return kFALSE;
  }
  //---
  Float_t trackTheta = helix->GetMomentum().Theta()*TMath::RadToDeg();

  Int_t emcEntries = fEmcCluster->GetEntriesFast();
  Int_t emcIndex = -1, emcModuleCorr = -1, emcNCrystals = -1, emcNBumps = -1;
  Float_t emcEloss = 0., emcElossCorr = 0., emcGLength = -1000;
  Float_t emcQuality = 1000000;
  Float_t chi2 = 0;
  TVector3 vertex(0., 0., 0.); TVector3 emcPos(0., 0., 0.);// TVector3 momentum(0., 0., 0.);

  // Cluster zenike moments
  Double_t Z20 = 0.0, Z53 = 0.0, secLatM = 0.00, E1 = 0., E9 = 0., E25 = 0.;

  for (Int_t ee = 0; ee<emcEntries; ee++)
    {
      PndEmcCluster *emcHit = (PndEmcCluster*)fEmcCluster->At(ee);

      if ( fIdeal )
	{
	  std::vector<Int_t> mclist = emcHit->GetMcList();
	  if (mclist.size()==0) continue;
	  if (mclist[0]!=pidCand->GetMcIndex()) continue;
	}

      if (emcHit->energy() < fCorrPar->GetEmc12Thr()) continue;
      Int_t emcModule = emcHit->GetModule();

      if (emcModule>4) continue; // skip FSC
      if ( (emcModule<3) && (helix->GetZ()>150.)  ) continue; // not consider tracks after emc barrel for BARREL
      //if ( (emcModule==3) && (helix->GetZ()<165.) ) continue; // consider tracks only from last gem plane for FWD
      if ( (emcModule==4) && (helix->GetZ()>-30.) ) continue; // consider tracks only ending at the back of STT for BKW

      emcPos = emcHit->where();
      if (fGeanePro)
	{ // Overwrites vertex if Geane is used
	  fGeanePropagator->SetPoint(emcPos);
	  fGeanePropagator->PropagateToPCA(1, 1);
	  vertex.SetXYZ(-10000, -10000, -10000); // reset vertex
	  FairTrackParH *fRes= new FairTrackParH();
	  Bool_t rc =  fGeanePropagator->Propagate(helix, fRes, fPidHyp*pidCand->GetCharge()); // First propagation at module
	  if (!rc) continue;

	  emcGLength = fGeanePropagator->GetLengthAtPCA();
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
	  emcEloss = emcHit->energy();
	  emcElossCorr = fEmcCalibrator->Energy(emcHit);
	  emcModuleCorr = emcModule;
	  emcNCrystals = emcHit->NumberOfDigis();
          emcNBumps  = emcHit->NBumps();
	  Z20 = emcHit->Z20();// Z_{n = 2}^{m = 0}
	  Z53 = emcHit->Z53();// Z_{n = 5}^{m = 3}
	  secLatM = emcHit->LatMom();
	  if (fEmcDigi)
	    {
	      PndEmcClusterEnergySums esum(*emcHit, fEmcDigi);
	      E1  = esum.E1();
	      E9  = esum.E9();
	      E25 = esum.E25();
	    }
	}

      if ( (fClusterQ[ee]<0) || (dist < fClusterQ[ee]))
	// If the track-emc distance is less than the previous stored value (or still not initialized)
	{
	  fClusterQ[ee] = dist; // update the param
	}

      if (fDebugMode){
	  Float_t ntuple[] = {static_cast<Float_t>(vertex.X()), static_cast<Float_t>(vertex.Y()), static_cast<Float_t>(vertex.Z()), static_cast<Float_t>(vertex.Phi()),
			      static_cast<Float_t>(helix->GetMomentum().Mag()), static_cast<Float_t>(helix->GetQ()), static_cast<Float_t>(helix->GetMomentum().Theta()), static_cast<Float_t>(helix->GetZ()),
			      static_cast<Float_t>(emcPos.X()), static_cast<Float_t>(emcPos.Y()), static_cast<Float_t>(emcPos.Z()), static_cast<Float_t>(emcPos.Phi()),
			      dist, static_cast<Float_t>(vertex.DeltaPhi(emcPos)), static_cast<Float_t>(emcHit->energy()), emcGLength, static_cast<Float_t>(emcModule)};
	// Float_t ntuple[] = {vertex.X(), vertex.Y(), vertex.Z(), vertex.Phi(),
	// 		    helix->GetMomentum().Mag(), helix->GetQ(), helix->GetMomentum().Theta(), helix->GetZ(),
	// 		    emcPos.X(), emcPos.Y(), emcPos.Z(), emcPos.Phi(),
	// 		    dist, vertex.DeltaPhi(emcPos), emcHit->energy(), emcGLength, emcModule};
	emcCorr->Fill(ntuple);
      }
    }// End for(ee = 0;)

  if ( (emcQuality < fCorrPar->GetEmc12Cut()) || ( fIdeal && emcIndex!=-1) ){
    fClusterList[emcIndex] = kTRUE;
    pidCand->SetEmcQuality(emcQuality);
    pidCand->SetEmcRawEnergy(emcEloss);
    pidCand->SetEmcCalEnergy(emcElossCorr);
    pidCand->SetEmcIndex(emcIndex);
    pidCand->SetEmcModule(emcModuleCorr);
    pidCand->SetEmcNumberOfCrystals(emcNCrystals);
    pidCand->SetEmcNumberOfBumps(emcNBumps);
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

ClassImp(PndPidCorrelator)

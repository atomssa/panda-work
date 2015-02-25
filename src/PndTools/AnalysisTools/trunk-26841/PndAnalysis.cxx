// PndAnalysis
// Needs a FairRunAna set up in the macro for file & parameter I/O

#include "PndAnalysis.h"

#include <string>
#include <iostream>

using std::cout;
using std::endl;

//Root stuff
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "RhoParticleSelectorBase.h"
#include "FairRecoCandidate.h"

//RHO stuff
#include "RhoFactory.h"
#include "RhoCandidate.h"
#include "RhoCandList.h"

#include "FairTrackParP.h"
#include "FairTrackParH.h"
#include "FairGeanePro.h"
//#include "FairRunAna.h"
#include "FairField.h"

#include "PndTrack.h"
#include "PndPidCandidate.h"
#include "PndPidProbability.h"
#include "PndAnaPidSelector.h"
#include "PndAnaPidCombiner.h"
#include "PndMCTrack.h"
#include "PndAnaCovTool.h" // using a cov matrix tool
#include "PndAnalysisCalcTools.h"
#include "PndPidBremCorrected4Mom.h"

ClassImp ( PndAnalysis );

PndAnalysis::PndAnalysis ( TString tname1, TString tname2, TString algnamec, TString algnamen ) :
    fRootManager ( FairRootManager::Instance() ),
    fPidSelector ( 0 ),
    fEvtCount ( 0 ),
    fChainEntries ( 0 ),
    fEventRead ( false ),
    fBuildMcCands ( false ),
    fVerbose(0),
    fPhotosMax(0), fPhotosThresh(0.05),
    fChargedPidName ( algnamec ),
    fNeutralPidName ( algnamen ),
    fTracksName ( tname1 ),
    fTracksName2 ( tname2 )
{
  if ( 0 == fRootManager ) {
    std::cout << "-E- PndAnalysis: RootManager not instantiated!" << std::endl;
    return;
  }

  Init();
}

PndAnalysis::~PndAnalysis()
{
  if ( 0!=fPidSelector ) {
    delete fPidSelector;
  }
}

TClonesArray* PndAnalysis::ReadTCA ( TString tcaname )
{
  TClonesArray* tca = ( TClonesArray* ) fRootManager->GetObject ( tcaname.Data() );
  if ( ! tca ) {
    std::cout << "-I- PndAnalysis::ReadTCA(): No "<<tcaname.Data() <<" array found." << std::endl;
  }

  return tca;
}

void PndAnalysis::Init()
{
  Reset();

  //read arrays
  fChargedCands = ReadTCA ( "PidChargedCand" );
  fBremCorr = ReadTCA( "BremCorrected4Mom" );
  fChargedProbability = ReadTCA ( fChargedPidName.Data() );
  fNeutralCands = ReadTCA ( "PidNeutralCand" );
  fNeutralProbability = ReadTCA ( fNeutralPidName.Data() );
  // -- Barrel Part
  std::cout << "-I- PndAnalysis::Init(): Trying "<<fTracksName.Data() <<" now." << std::endl;
  fTracks = ReadTCA ( fTracksName );

  if ( !fTracks ) {
    std::cout << "-I- PndAnalysis::Init(): Trying SttMvdGenTrack now." << std::endl;
    fTracks = ReadTCA ( "SttMvdGenTrack" );
  }

  if ( !fTracks ) {
    std::cout << "-I- PndAnalysis::Init(): Trying SttMvdGemGenTrack now." << std::endl;
    fTracks = ReadTCA ( "SttMvdGemGenTrack" );
  }

  if ( !fTracks ) {
    std::cout << "-I- PndAnalysis::Init(): Trying BarrelGenTrack now." << std::endl;
    fTracks = ReadTCA ( "BarrelGenTrack" );
  }

  if ( !fTracks ) {
    std::cout << "-E- PndAnalysis::Init(): No track inpt array. Make a Zero TclonesArray." << std::endl;
  } else {
    fTracks->Print();
  }

  // -- Forward part
  std::cout << "-I- PndAnalysis::Init(): Second: Trying "<<fTracksName2.Data() <<" now." << std::endl;

  fTracks2 = ReadTCA ( fTracksName2 );

  if ( !fTracks2 ) {
    std::cout << "-I- PndAnalysis::Init(): Second: Trying FTSGenTrack now." << std::endl;
    fTracks2 = ReadTCA ( "FTSGenTrack" );
  }

  if ( !fTracks2 ) {
    std::cout << "-I- PndAnalysis::Init(): Second: Trying FTSTrkIdeal now." << std::endl;
    fTracks2 = ReadTCA ( "FTSTrkIdeal" );
  }

  // -- MC Tracks
  fBuildMcCands = false;

  fMcCands = 0;//ReadTCA ( "PndMcTracks" ); // try already built RhoCandidates

  if ( ! fMcCands ) {
    std::cout << "-I- PndAnalysis::Init(): Trying mc stack now." << std::endl;
    fMcTracks = ( TClonesArray* ) fRootManager->GetObject ( "MCTrack" );

    if ( ! fMcTracks && fVerbose ) {
      std::cout << "-W- PndAnalysis::Init(): No \"MCTrack\" array found. No MC info available." << std::endl;
    }

    fMcCands =new TClonesArray ( "RhoCandidate" );

    // next line commented by KG, 07/2012
    fRootManager->Register ( "PndMcTracks","PndMcTracksFolder", fMcCands, kFALSE );
    fBuildMcCands = true;
  }

  fChainEntries = ( fRootManager->GetInChain() )->GetEntries();

  //TODO default constructor here?
  fPidCombiner = new PndAnaPidCombiner();
  fPidSelector = new PndAnaPidSelector();

  fPdg = TDatabasePDG::Instance();

}

void PndAnalysis::Rewind()
{
  fEvtCount=0;
}

void PndAnalysis::Cleanup()
{
    // do a safe cleanup
  fAllCandList.Cleanup();
  fChargedCandList.Cleanup();
  fNeutralCandList.Cleanup();
  fMcCandList.Cleanup();
  RhoFactory::Instance()->Reset();

}

void PndAnalysis::ReadCandidates()
{
  ReadRecoCandidates();
  BuildMcCands();
  // now fill carged and neutral lists.
  // MC association to reconstructed particles done at this point and copying is ok
  for(int ik=0;ik<fAllCandList.GetLength();ik++)
  {
    if(fAllCandList[ik]->GetCharge()==0){
      fNeutralCandList.Add ( fAllCandList[ik] );
    }else{
      fChargedCandList.Add ( fAllCandList[ik] );
    }
  }
  return;
}

void PndAnalysis::GetEventInTask()
{
  Cleanup();
  ReadCandidates();
  //std::cout<<"Marke A"<<std::endl;
  return;
}

Int_t PndAnalysis::GetEvent ( Int_t n )
{
  Cleanup();

  if ( n>=0 ) {
    fEvtCount=n+1;
  } else {
    fEvtCount++;
  }

  if ( fEvtCount>fChainEntries ) {
    fEvtCount=fChainEntries;
    if(fVerbose) Info("PndAnalysis::GetEvent()","Maximum number of entried in the file chain reached: %i.",fEvtCount);
    return 0;
  }
  fRootManager->ReadEvent ( fEvtCount-1 );

  ReadCandidates();

  if(fVerbose) Info("PndAnalysis::GetEvent()","Finished loading event fEvtCount=%i.",fEvtCount);
  //std::cout<<"PndAnalysis::fAllCandList:   "<<fAllCandList<<std::endl;
  return fEvtCount;
}

FairMCEventHeader* PndAnalysis::GetEventHeader()
{
  if ( !fEventRead ) {
    fRootManager->ReadEvent ( fEvtCount-1 );
    fEventRead=kTRUE;
  }

  FairMCEventHeader*  evthead = ( FairMCEventHeader* ) FairRootManager::Instance()->GetObject ( "MCEventHeader." );

  return evthead;
}

Bool_t PndAnalysis::FillList ( RhoCandList& l, TString listkey, TString pidTcaNames )
{
  // Reads the specified List for the current event
  l.Cleanup();

  // Set which PID information should be used.
  if ( pidTcaNames!="" ) {
    fPidCombiner->SetTcaNames ( pidTcaNames );
  } else {
    fPidCombiner->SetDefaults();
  }

  // Get or build Monte-Carlo truth list
  if ( listkey=="McTruth" ) {
    return GetMcCandList(l);
  }

  //Info("PndAnalysis::FillList","key=%s",listkey.Data());
  // acceleration: just give the large lists directly
  if ( listkey=="All" ) {
    fPidCombiner->Apply ( fAllCandList );
    l=fAllCandList;
    return kTRUE;
  }

  if ( listkey=="Neutral" ) {
    fPidCombiner->Apply(fNeutralCandList);
    l=fNeutralCandList;
    return kTRUE;
  }

  if ( listkey=="Charged" ) {
    fPidCombiner->Apply ( fChargedCandList );
    l=fChargedCandList;
    return kTRUE;
  }

  const bool doBremCorr = listkey.Contains("Brem");
  if (doBremCorr) listkey.ReplaceAll("Brem","");

  // Real selection requested:
  // set the base list for the PID list maker
  Bool_t checkcrit = fPidSelector->SetCriterion ( listkey );
  if (!checkcrit) return kFALSE;

  if ( listkey.Contains ( "Electron" ) ||listkey.Contains ( "Muon" ) ||listkey.Contains ( "Pion" )
       || listkey.Contains ( "Kaon" ) ||listkey.Contains ( "Proton" )
       || listkey.Contains ( "Plus" ) ||listkey.Contains ( "Minus" ) ||listkey.Contains ( "Charged" ) ) {

    if ( doBremCorr ) {
      if (fBremCorr==0) {
	if(fVerbose)
	  Warning("PndAnalysis::FillList","Brem requested but no PndPidBremCorrected4Mom found on input file. Brem Correction can't be done.");
      } else {
    	for (int j=0; j<fChargedCandList.GetLength(); ++j) {
    	  int trk_id = fChargedCandList[j]->GetTrackNumber();
    	  int nBremCorr = fBremCorr->GetEntriesFast();
    	  if (nBremCorr!=fChargedCandList.GetLength())
	    if(fVerbose) Warning("PndAnalysis::FillList","Warning: BermCorr list size diff. from chargeCandList");
    	  PndPidBremCorrected4Mom *bremCorr = (PndPidBremCorrected4Mom*) fBremCorr->At(trk_id);
    	  fChargedCandList[j]->SetP3(bremCorr->GetMomentum());
    	}
      }
    }

    l=fChargedCandList;
    fPidCombiner->Apply ( l );
    l.Select(fPidSelector);
    //fPidSelector->Select ( fChargedCandList,l );
    return kTRUE;
  }

  if ( listkey.Contains ( "Neutral" ) ) {
    fPidCombiner->Apply ( fNeutralCandList );
    fPidSelector->Select ( fNeutralCandList,l );
    return kTRUE;
  }


  Error ( "FillList", "Unknown list key: %s",listkey.Data() );
  return kFALSE;
}

Bool_t PndAnalysis::GetMcCandList(RhoCandList& l)
{
  l.Clear();
  // Put all candidates from the mctruth list to the candlist and set the mother-daughter relations
  if ( !fMcCands ) return kFALSE;

  RhoCandidate* truth=0;
  for (int i=0;i<fMcCands->GetEntriesFast();i++)
  {
    // copy candidates via put
    truth = (RhoCandidate*) fMcCands->At(i);
    l.Put(truth);
  }

  // now set genealogy inside the list
  RhoCandidate* truthmother=0;
  for (int k=0;k<l.GetLength();k++)
  {
    // get mother track
    PndMCTrack* part = (PndMCTrack*) fMcTracks->At(k);
    Int_t mcMotherID = part->GetMotherID();
    if(mcMotherID<0) mcMotherID=part->GetSecondMotherID();

    // SetMotherLink does the deep mother-daughter relation
    if (mcMotherID<0) continue; // no mother there, go on...
    if (mcMotherID>=l.GetLength()) continue; // something bad hapened to the indices

    // do the linking
    truthmother = (RhoCandidate*) l[mcMotherID];
    l[k]->SetMotherLink(truthmother, false);
  }

  return kTRUE;
}

Int_t PndAnalysis::GetEntries()
{
  if ( fRootManager ) {
    return ( fRootManager->GetInChain() )->GetEntries();
  } else {
    return 0;
  }
}

void PndAnalysis::ReadRecoCandidates()
{
  UInt_t _uid=0;
  fAllCandList.Cleanup();
  fChargedCandList.Cleanup();
  fNeutralCandList.Cleanup();
  if ( fNeutralCands ) {
    for ( Int_t i1=0; i1<fNeutralCands->GetEntriesFast(); i1++ ) {
        FairRecoCandidate* mic = ( FairRecoCandidate* ) fNeutralCands->At ( i1 );
        _uid++; // uid will start from 1
        RhoCandidate tc ( *mic,_uid );
        tc.SetTrackNumber ( -1 );//(i1);
	tc.SetType( 22 );     // default PDG code for neutrals is gamma = 22
        // TODO: Do we want to set something here? It is neutrals anyway.

        if ( 0!=fNeutralProbability && i1<fNeutralProbability->GetEntriesFast() ) {
          PndPidProbability* neuProb = ( PndPidProbability* ) fNeutralProbability->At ( i1 );

          if ( neuProb == 0 ) {
            Error ( "FillList", "Neutral PID Probability object not found, skip setting pid for candidate %i.",i1 );
            continue;
          }

          // numbering see PndPidListMaker
          tc.SetPidInfo ( 0,neuProb->GetElectronPidProb() );
          tc.SetPidInfo ( 1,neuProb->GetMuonPidProb() );
          tc.SetPidInfo ( 2,neuProb->GetPionPidProb() );
          tc.SetPidInfo ( 3,neuProb->GetKaonPidProb() );
          tc.SetPidInfo ( 4,neuProb->GetProtonPidProb() );
        }
        fAllCandList.Add ( &tc );
      }
    } else {
      if(fVerbose) Warning("PndAnalysis::ReadRecoCandidates()","No neutral reco array found.");
    }

    if ( fChargedCands) {
      for ( Int_t i2=0; i2<fChargedCands->GetEntriesFast(); i2++ ) {
        _uid++; // uid will start from (n_neutrals + 1)
        FairRecoCandidate* mic = ( FairRecoCandidate* ) fChargedCands->At ( i2 );
        RhoCandidate tc ( *mic,_uid );
        tc.SetTrackNumber ( i2 );
	tc.SetType( tc.Charge()*211 );  // default PDG code for charged is pi = +-211
        // TODO: Check that no i+1 is requested anymore elsewhere!!!
        fAllCandList.Add ( &tc );
      }
    } else {
      if(fVerbose) Warning("PndAnalysis::ReadRecoCandidates()","No charged reco array found.");
    }
  return;
}



void PndAnalysis::BuildMcCands()
{
  int i;
  // Make Monte-carlo truth candidates by the reconstructed particles up to the initial state (if available)
  if ( !fMcCands ){
    Warning("PndAnalysis::BuildMcCands","No array to store candidates...");
    return;
  }
  if ( !fBuildMcCands ) {
    if(fVerbose) Info("PndAnalysis::BuildMcCands","No mc to build...");
    return;
  }
  if ( fMcCands->GetEntriesFast() != 0 ) {
    fMcCands->Delete();
  }

  if ( fMcTracks == 0 ) {
    Error ( "BuildMcCands","MC track Array does not exist." );
    return;
  }

  //loop all MCTracks
  for (i=0;i<fMcTracks->GetEntriesFast();i++)
  {
    // fetch particle properties
    PndMCTrack* part = (PndMCTrack*) fMcTracks->At(i);
    TLorentzVector p4 = part->Get4Momentum();
    TVector3    stvtx = part->GetStartVertex();
    TParticlePDG* ppdg = fPdg->GetParticle(part->GetPdgCode());
    double charge=0.0;
    if ( ppdg ) {
       charge=ppdg->Charge();
    } else if (fVerbose) {
       cout <<"-W- CreateMcCandidate: strange PDG code:"<<part->GetPdgCode() <<endl;
    }
       if ( fabs(charge) >2 ) {
       charge/=3.;
    }
    // create mc candidate
    RhoCandidate* pmc=new ( (*fMcCands)[i] ) RhoCandidate(p4,charge);
    pmc->SetPos(stvtx);
    pmc->SetType(part->GetPdgCode()); //this overwrites our generator's mass information
    pmc->SetP4(p4);
    pmc->SetTrackNumber(i);
  }

  //write correctly assigned copy of mc truth candidates
  GetMcCandList(fMcCandList);

  // Assign MC truth to reconstructed allCandnds
  RhoCandidate* truth=0;
  for(int icand=0;icand<fAllCandList.GetLength();icand++){
    RhoCandidate* currentcand=fAllCandList.Get(icand);
    //   get reco candidate
    FairRecoCandidate* reco = currentcand->GetRecoCandidate();
    if(!reco) {
      if (fVerbose) Info("BuildMcCands","reco object to candidate %i (%p) missing.",icand,currentcand);
      continue;
    }
    // get the mctruth
    Int_t mcidx = reco->GetMcIndex();
    if (mcidx>fMcCandList.GetLength() || mcidx<0) continue;
    truth = fMcCandList[mcidx];
    currentcand->SetMcTruth(truth);
    if(fVerbose)Info("PndAnalysis::BuildMcCands()","Now setting truth index %i (%p) to candidate (uid=%i)", mcidx,truth,currentcand->Uid());
  }

  //std::cout<<"BuildMcCands():  "<<fAllCandList<<std::endl;
}

Bool_t PndAnalysis::PropagateToIp ( RhoCandidate* cand )
{
  //Propagate from the tracks first parameter set to the POCA from (0,0,0)
  TVector3 ip( 0.,0.,0. );
  return PropagateToPoint ( cand, ip );
}

Bool_t PndAnalysis::PropagateToZAxis ( RhoCandidate* cand )
{
  if ( !cand ) {
    Error ( "PropagateToZAxis","Candidate not found: %p",cand );
    return kFALSE;
  }

  PndPidCandidate* pidCand = (PndPidCandidate*)cand->GetRecoCandidate();

  PndTrack* track = ( PndTrack* ) fTracks->At ( pidCand->GetTrackIndex() );

  if ( !track ) {
    Warning ( "PropagateToZAxis","Could not find track object of index %d",pidCand->GetTrackIndex() );
    return kFALSE;
  }

  FairTrackParP tStart = track->GetParamFirst();

  return Propagator ( 2,tStart,cand );
}

Bool_t PndAnalysis::PropagateToPoint ( RhoCandidate* cand, TVector3 mypoint )
{
  //Propagate from the tracks first parameter set to the POCA from mypoint
  //The candidate is updated but the track not touched
  //Only the uncorrelated errors are propagated,
  //TODO: implement a real cov matrix
  if ( !cand ) {
    Error ( "PropagateToPoint","Candidate not found: %p",cand );
    return kFALSE;
  }

  PndPidCandidate* pidCand = (PndPidCandidate*)cand->GetRecoCandidate();

  PndTrack* track = ( PndTrack* ) fTracks->At ( pidCand->GetTrackIndex() );

  if ( !track ) {
    Warning ( "PropagateToPoint","Could not find track object of index %d",pidCand->GetTrackIndex() );
    return kFALSE;
  }

  FairTrackParP tStart = track->GetParamFirst();

  return Propagator ( 1,tStart,cand,mypoint );
}

Bool_t PndAnalysis::PropagateToPlane(RhoCandidate* cand, TVector3 origin, TVector3 dj, TVector3 dk)
{
  //Propagate from the tracks first parameter set to a (detector) plane
  //The candidate is updated but the track not touched
  //Only the uncorrelated errors are propagated,
  //TODO: implement a real cov matrix
  if ( !cand ) {
    Error ( "PropagateToPlane","Candidate not found: %p",cand );
    return kFALSE;
  }

  PndPidCandidate* pidCand = (PndPidCandidate*)cand->GetRecoCandidate();

  PndTrack* track = ( PndTrack* ) fTracks->At ( pidCand->GetTrackIndex() );

  if ( !track ) {
    Warning ( "PropagateToPlane","Could not find track object of index %d",pidCand->GetTrackIndex() );
    return kFALSE;
  }

  FairTrackParP tStart = track->GetParamFirst();
  return Propagator ( 3,tStart,cand,origin,kFALSE,kFALSE,dj,dk );
}


FairTrackParP PndAnalysis::GetFirstPar ( RhoCandidate* cand )
{
  if ( !cand ) {
    Error ( "GetFirstPar","Candidate not found: %p",cand );
    FairTrackParP dummy;
    return dummy;
  }

  PndPidCandidate* pidCand = (PndPidCandidate*)cand->GetRecoCandidate();

  PndTrack* track = ( PndTrack* ) fTracks->At ( pidCand->GetTrackIndex() );

  if ( !track ) {
    Warning ( "GetFirstPar","Could not find track object of index %d",pidCand->GetTrackIndex() );
    FairTrackParP dummy;
    return dummy;
  }

  FairTrackParP tStart = track->GetParamFirst();

  return tStart;
}

Bool_t PndAnalysis::ResetDaughters ( RhoCandidate* cand )
{
  Bool_t success=kTRUE;

  for ( Int_t daug =0; daug<cand->NDaughters(); daug++ ) {
    RhoCandidate* a=cand->Daughter ( daug );
    success = success && ResetCandidate ( a );
  }

  return success;
}

Bool_t PndAnalysis::ResetCandidate ( RhoCandidate* cand )
{
  FairTrackParP firstpar = GetFirstPar ( cand );
  Double_t globalCov[6][6];
  firstpar.GetMARSCov ( globalCov );
  TMatrixD err ( 6,6 );

  for ( Int_t ii=0; ii<6; ii++ ) for ( Int_t jj=0; jj<6; jj++ ) {
      err[ii][jj]=globalCov[ii][jj];
    }

  //if(fVerbose>2){ std::cout<<"MARS cov (px,py,pz,E,x,y,z): ";err.Print();}
  TLorentzVector lv = cand->P4();

  static PndAnaCovTool covTool; // external tool to convert from a 6x6 (p3,v) cov matrix to the 7x7(p4,v) cov matrix

  TMatrixD covPosMom = covTool.GetConverted7 ( covTool.GetFitError ( lv, err ) );

  //if(fVerbose>2){ std::cout<<"covPosMom (x,y,z,px,py,pz,E): ";covPosMom.Print();}

  cand->SetPosition ( firstpar.GetPosition() );

  cand->SetP3 ( firstpar.GetMomentum() ); // implicitly uses the candidates mass to set P4

  cand->SetCov7 ( covPosMom );

  return kTRUE;
}


Bool_t PndAnalysis::Propagator ( int mode, FairTrackParP& tStart, RhoCandidate* cand, TVector3 mypoint, Bool_t skipcov, Bool_t overwrite, TVector3 planej, TVector3 planek )
{
  //Propagate from the tracks first parameter set to the POCA from mypoint
  //The candidate is updated but the track not touched
  //Only the uncorrelated errors are propagated,
  //TODO: implement a real cov matrix

  Bool_t rc = kFALSE;
  static PndAnaCovTool covTool; // external tool to convert from a 6x6 (p3,v) cov matrix to the 7x7(p4,v) cov matrix
  FairGeanePro* geaneProp = new FairGeanePro();
  Int_t pdgcode = cand->PdgCode();

  if ( fVerbose>0 ) {
    cout<<"Try mode "<<mode<<" with pdgCode "<<pdgcode<<endl;
  }

  if ( fVerbose>2 ) {
    std::cout<<"Start Params are:"<<std::endl; tStart.Print();
  }

  Double_t startCov[6][6];

  tStart.GetMARSCov ( startCov );
  TMatrixD errst ( 6,6 );

  for ( Int_t ii=0; ii<6; ii++ ) for ( Int_t jj=0; jj<6; jj++ ) {
      errst[ii][jj]=startCov[ii][jj];
    }

  if ( fVerbose>2 ) {
    std::cout<<"Start MARS cov: "; errst.Print();
  }

  if ( 1==mode ) { // to point
    geaneProp->BackTrackToVertex(); //set where to propagate
    geaneProp->SetPoint ( mypoint );
  } else if ( 2==mode ) { // to line
    geaneProp->PropagateToPCA ( 2, -1 );// track back to z axis
    TVector3 ex1 ( 0.,0.,-50. ); // virtual wire of arbitrarily chosen size
    TVector3 ex2 ( 0.,0.,100. );
    geaneProp->SetWire ( ex1,ex2 );
  } else if ( 3==mode ) { // to plane
    geaneProp->PropagateToPlane(mypoint,planej,planek);
  } else {
    Error ( "Propagator()","Use mode 1 (to a TVector3) or mode 2 (to z axis) or mode 3 (to plane). (Mode=%i)",mode );
    return kFALSE;
  }

  if(skipcov) geaneProp->PropagateOnlyParameters();

  FairTrackParH* myResult=0;
  // now we propagate
  if(mode==3){
    FairTrackParP* tResult = new FairTrackParP();
    rc = geaneProp->Propagate ( &tStart, tResult,pdgcode );
    myResult = new FairTrackParH(*tResult);
  }else{
    myResult = new FairTrackParH();
    FairTrackParH* myStart = new FairTrackParH ( tStart );
    rc = geaneProp->Propagate ( myStart, myResult,pdgcode );
  }

  if ( !rc ) {
    if ( fVerbose>0 ) {
      Warning ( "Propagator()","Geane propagation failed" );
    }
    return kFALSE;
  }

  int ierr=0;
  TVector3 di = myResult->GetMomentum();
  di.SetMag ( 1. );
  TVector3 dj = di.Orthogonal();
  TVector3 dk = di.Cross ( dj );
  FairTrackParP* myParab = new FairTrackParP ( myResult, dj, dk, ierr );

  TVector3 pos( myResult->GetX(),myResult->GetY(),myResult->GetZ() ); // I want to be sure...
  cand->SetPosition( pos );
  cand->SetP3( myResult->GetMomentum() ); // implicitly uses the candidates mass to set P4

  //printout for checks
  if ( fVerbose>1 ) {
    TVector3 vecdiff=tStart.GetPosition() - myResult->GetPosition();
    std::cout<<"position start     :";
    tStart.GetPosition().Print();
    std::cout<<"position ip        :";
    myResult->GetPosition().Print();
    std::cout<<"position difference:";
    vecdiff.Print();
    vecdiff=tStart.GetMomentum()-myResult->GetMomentum();
    std::cout<<"momentum start     :";
    tStart.GetMomentum().Print();
    std::cout<<"momentum ip        :";
    myResult->GetMomentum().Print();
    std::cout<<"momentum difference:";
    vecdiff.Print();
  }

  if(kFALSE==skipcov){
    Double_t globalCov[6][6];
    myParab->GetMARSCov ( globalCov );
    TMatrixD err ( 6,6 );

    for ( Int_t ii=0; ii<6; ii++ ) for ( Int_t jj=0; jj<6; jj++ ) {
        err[ii][jj]=globalCov[ii][jj];
      }

    if ( fVerbose>2 ) {
      std::cout<<"MARS cov (px,py,pz,E,x,y,z): "; err.Print();
    }

    TLorentzVector lv = cand->P4();

    TMatrixD covPosMom = covTool.GetConverted7 ( covTool.GetFitError ( lv, err ) );

    if ( fVerbose>2 ) {
      std::cout<<"covPosMom (x,y,z,px,py,pz,E): "; covPosMom.Print();
    }

    cand->SetCov7 ( covPosMom );
  }
//  rc = PndAnalysisCalcTools::FillHelixParams(cand,skipcov);
//  if (!rc) {Warning("Propagator()","P7toHelix failed"); return kFALSE;}

  if ( fVerbose>2 ) {
    std::cout<<" :::::::::::  Printout in PndAnalysis::Propagator() :::::::::::  "<<std::endl;

    //std::cout<<"Start Params:"<<std::endl;
    //myStart->Print();

    std::cout<<"SC system params:"
             <<"\nq/p    = "<<myResult->GetQp()
             <<"\nLambda = "<<myResult->GetLambda()
             <<"\nPhi    = "<<myResult->GetPhi()
             <<"\nX_sc   = "<<myResult->GetX_sc()
             <<"\nY_sc   = "<<myResult->GetY_sc()
             <<"\nZ_sc   = "<<myResult->GetZ_sc()
             <<std::endl;

    std::cout<<"some values:"
             <<"\n Z0 ?= z_sc / cos(lambda)  = "<< myResult->GetZ_sc() / cos ( myResult->GetLambda() )
             <<"\n Z0 ?= sqrt(x_sc^2+z_sc^2) = "<<sqrt ( myResult->GetX_sc() *myResult->GetX_sc() +myResult->GetZ_sc() *myResult->GetZ_sc() )
             <<std::endl;
  }

  // COMMENT:
  // When taking the Trackparams in the POCA to the z-axis, the SC system from GEANE matches the common helix params easier, i.e:
  // D0 = y_sc and Z0 = sqrt(x_sc^2 + z_sc^2) = z_sc*tan(Lambda)

  if(overwrite){
    if ( fVerbose>1 ) {
      Info ( "Propagator  ","overwriting start parameter state with result");
    }
    tStart.SetTrackPar(myParab->GetV(), myParab->GetW(),
                       myParab->GetTV(), myParab->GetTW(),
                       myParab->GetQp(), myParab->GetCov(),
                       myParab->GetOrigin(),
                       myParab->GetIVer(),
                       myParab->GetJVer(),
                       myParab->GetKVer(),
                       myParab->GetSPU()
                       );
  }

  if ( fVerbose>1 ) {
    Info ( "Propagator  ","Succsess=%i",rc );
  }

  return kTRUE;
}



//////////////////////////////////////////////////////
//  MC Truth matching                               //
//////////////////////////////////////////////////////



Bool_t PndAnalysis::McTruthMatch(RhoCandidate* cand, Int_t level, bool verbose)
{
    return MctMatch(cand,fMcCandList,level,verbose);
}

Int_t PndAnalysis::McTruthMatch(RhoCandList& list, Int_t level, bool verbose)
{
  Int_t ifound = 0;
  for(int icand=0;icand<list.GetLength();icand++){
    if( true == MctMatch(list[icand],fMcCandList,level,verbose) ) {
      ifound++;
    }
  }
  return ifound;
}

Bool_t PndAnalysis::MctMatch ( RhoCandidate* c, RhoCandList& mct, Int_t level, bool verbose )
{
  Int_t nd  = c->NDaughters();
  Int_t pdg = c->PdgCode();
  Int_t nmct = mct.GetLength();

  if ( 0==nd ) { // final state particle
    RhoCandidate* mccnd = c->GetMcTruth();
    if ( !mccnd ) {
      if(verbose) Info("PndMcTruthMatch::MctMatch","rejected final state by nonexistent mc truth pointer");
      return false;
    }
    if ( mccnd->PdgCode() == pdg ) {
      if(verbose) Info("PndMcTruthMatch::MctMatch","accepted final state by PDG code (pdg=%i)",pdg);
      //if(verbose) std::cout<<*c<<std::endl;
      return true;
    } else {
      if(verbose) Info("PndMcTruthMatch::MctMatch","rejected final state by PDG Code (pdg=%i|mcpdg=%i)",pdg, mccnd->PdgCode());
      return false;
    }
  }

    // check recursively whether all daughter trees match
  for ( Int_t i=0; i<nd; i++ ) {
    if ( !MctMatch (  c->Daughter ( i ) , mct, level, verbose ) ) {
      if(verbose) Info("PndMcTruthMatch::MctMatch","rejected composite (pdg=%i) by non-matching daughter: idau=%i",pdg,i);
      return false;
    }
  }

  // ***
  // *** MATCH LEVEL 0 reached: only PID of final state particles are matched
  // ***
  if ( 0==level ) {
    return true;
  }
  //if(verbose)std::cout<<"going Level 1"<<std::endl;
  // reset mc truth pointer
  RhoCandidate* mccnd = c->GetMcTruth();
  if(mccnd) {
    if(verbose) Warning("PndMcTruthMatch::MctMatch","Existing MC truth found. Will reset it now.");
    c->SetMcTruth(0);
    //return false;
  }

  // find this particle's truth in the mc decay tree
  RhoCandidate* dauzero = c->Daughter(0);
  if (!dauzero) {
    if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by not existing daughter zero");
    return false;
  }
  RhoCandidate* mcdauzero = dauzero->GetMcTruth();
  if (!mcdauzero) {
    if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by not existing MC truth of daughter zero");
    return false;
  }
  RhoCandidate* mcdauzeromother=mcdauzero->TheMother();
  if (!mcdauzeromother) {
    if(verbose) {Info("PndMcTruthMatch::MctMatch","rejected by not existing mother of MC truth of daughter zero");
    cout <<*mcdauzero<<endl;}

    return false;
  }

  //now check the tree structure:

  int  ndaudiff  = mcdauzeromother->NDaughters() - nd;

  // we still might accept if only photos photons are missing in the reco tree
  // if overall difference is small enough, we check further below what particles are missing
  if( ndaudiff<0 || ndaudiff>fPhotosMax )
  {
    if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by differing daughter count: cand:%i mc:%i",c->NDaughters(),mcdauzeromother->NDaughters());
    return false;
  }
  // count daughter photons (pdg==22) from MC mother with E < photos_thresh
  int nphall = 0;
  if (ndaudiff>0)
    for (int idau=0;idau<mcdauzeromother->NDaughters();++idau)
      if (mcdauzeromother->Daughter(idau)->PdgCode()==22 && mcdauzeromother->Daughter(idau)->E()<fPhotosThresh) nphall++;

  // reco'd photons with E<photos_thresh
  int nphreco = 0;
  //  now if all daughters MC-Mother is the same
  for(int idau=1;idau<nd;idau++){
    // look if all daughters mc mothers are the same
    RhoCandidate* dau = c->Daughter(idau);
    if (!dau) {
      if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by not existing daughter %i",idau);
      return false;
    }
    RhoCandidate* mcdau = dau->GetMcTruth();
    if (!mcdau) {
      if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by not existing MC truth of daughter %i",idau);
      //if(verbose) std::cout<<*dau<<std::endl;
      return false;
    }
    RhoCandidate* mcdaumother=mcdau->TheMother();
    if (!mcdaumother) {
      if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by not existing mother of MC truth of daughter %i",idau);
      return false;
    }
    if(mcdaumother!=mcdauzeromother){
      if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by mc mother of daughter %i(%p) is different to mc mother of daughter zero(%p) -> Tree does not match"
       			,idau,mcdaumother,mcdauzeromother);
      return false;
    }

    // daughter was correctly matched; now check for low energetic photon for photos accept
    if (ndaudiff>0 && mcdau->PdgCode()==22 && mcdau->E()<fPhotosThresh) nphreco++;
  }
  // difference in #phot<max allowed and has to be exactly the total number of particle difference
  if ( (nphall-nphreco)>fPhotosMax || (nphall-nphreco)!=ndaudiff )
  {
	if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by differing daughter count not being photos photons: cand:%i mc:%i",c->NDaughters(),mcdauzeromother->NDaughters());
	return false;
  }

  // ***
  // *** MATCH LEVEL 1: PID of the leaves and tree topology are matched
  // ***
  if ( 1==level ) {
    if (verbose) Info("PndMcTruthMatch::MctMatch","accepted composite(pdg=%i) on match-level 1",pdg);
    c->SetMcTruth(mcdauzeromother);
    if (verbose)  cout <<*c->GetMcTruth()<<endl;
    return true;
  }
  //if(verbose)std::cout<<"going Level 2"<<std::endl;
  // check whether all daughter's mother has correct PDG code
  if ( pdg != mcdauzeromother->PdgCode() ) {
    if(verbose) Info("PndMcTruthMatch::MctMatch","rejected by nonmatching pdg code in tree (pdgcode|mcpdgcode) (%i|%i)",pdg,mcdauzeromother->PdgCode());
    return false;
  }
  // ***
  // *** MATCH LEVEL 2: PID of leaves, tree topology and intermediate particle types are matched
  // ***
  if (verbose) Info("PndMcTruthMatch::MctMatch","accepted composite(pdg=%i) on match-level 2", pdg);
  c->SetMcTruth(mcdauzeromother);
  if (verbose) cout <<*c->GetMcTruth()<<endl;
  return true;  // c's tree matches!
}

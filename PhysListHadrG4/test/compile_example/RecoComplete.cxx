#include "RecoComplete.h"

#include <PndEmcMakeCluster.h>
#include <PndEmcHdrFiller.h>
#include <PndEmcMakeBump.h>
#include <PndEmcMakeRecoHit.h>
#include <PndDchFindTracks.h>
#include <PndDchTrackFinderIdealCylHit.h>
#include <PndDchMatchTracks.h>
#include <PndMvdClusterTask.h>
#include <PndGemFindTracks.h>
#include <PndGemTrackFinderOnHits.h>
#include <PndGemTrackFinderQA.h>
#include <PndSttTrackFinderIdeal.h>
#include <PndSttFindTracks.h>
#include <PndSttMatchTracks.h>
#include <PndSttTrackFitter.h>
#include <PndSttHelixTrackFitter.h>
#include <PndSttFitTracks.h>
#include <PndSttHelixHitProducer.h>
#include <PndTpcClusterFinderTask.h>
#include <PndTpcIdealTrackingTask.h>

#include <KalmanTask.h>
#include <TrackFitStatTask.h>

#include <FairGeane.h>
#include <FairBoxGenerator.h>
#include <FairEvtGenGenerator.h>
#include <FairModule.h>
#include <FairParAsciiFileIo.h>
#include <FairParRootFileIo.h>
#include <FairPrimaryGenerator.h>
#include <FairRunAna.h>
#include <FairRunSim.h>
#include <FairRuntimeDb.h>

#include <TROOT.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include <iostream>


using namespace std;


void RecoComplete(TString const &mcFile, TString const &dgFile,
		  TString const &parFile, TString const &digiFile,
		  TString const &outFile)
{
  gDebug = 0;

  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0; // just forget about it, for the moment
  
  // Number of events to process
  Int_t nEvents = 0;  // if 0 all the vents will be processed
  
  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  // ------------------------------------------------------------------------
  
  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(mcFile);
  fRun->AddFriend(dgFile);
  fRun->SetOutputFile(outFile);


  FairGeane *Geane = new FairGeane();
  fRun->AddTask(Geane);


  // -----  Parameter database   --------------------------------------------
  TString emcDigiFile = gSystem->Getenv("VMCWORKDIR");
  emcDigiFile += "/macro/params/";
  emcDigiFile += digiFile;
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(parFile.Data());
  
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(emcDigiFile.Data(),"in");
        
  rtdb->setFirstInput(parInput1);
  rtdb->setSecondInput(parIo1);
  

  // -----   EMC hit producers   ---------------------------------
  // The file name should be the same of the geometry file which was used for the simulation
  
  PndEmcMakeCluster* emcMakeCluster= new PndEmcMakeCluster(iVerbose);
  emcMakeCluster->SetStorageOfData(kFALSE);
  fRun->AddTask(emcMakeCluster);
  PndEmcHdrFiller* emcHdrFiller = new PndEmcHdrFiller();
  fRun->AddTask(emcHdrFiller); // ECM header
  PndEmcMakeBump* emcMakeBump= new PndEmcMakeBump();
  fRun->AddTask(emcMakeBump);
  PndEmcMakeRecoHit* emcMakeRecoHit= new PndEmcMakeRecoHit();
  fRun->AddTask(emcMakeRecoHit);

  // trackfinding ....
  PndSttTrackFinderIdeal* sttTrackFinder = new PndSttTrackFinderIdeal(iVerbose);
  PndSttFindTracks* sttFindTracks = new PndSttFindTracks("Track Finder", "FairTask", sttTrackFinder, iVerbose);
  sttFindTracks->AddHitCollectionName("STTHit", "STTPoint");
  fRun->AddTask(sttFindTracks);
  // trackmatching ....
  PndSttMatchTracks* sttTrackMatcher = new PndSttMatchTracks("Match tracks", "STT", iVerbose);
  sttTrackMatcher->AddHitCollectionName("STTHit", "STTPoint");
  fRun->AddTask(sttTrackMatcher);  
  // trackfitting ....
  PndSttHelixTrackFitter* sttTrackFitter = new PndSttHelixTrackFitter(0);
  PndSttFitTracks* sttFitTracks = new PndSttFitTracks("STT Track Fitter", "FairTask", sttTrackFitter); 
  sttFitTracks->AddHitCollectionName("STTHit");
  fRun->AddTask(sttFitTracks);
  // helix hit production ....
  PndSttHelixHitProducer* sttHHProducer = new PndSttHelixHitProducer();
  fRun->AddTask(sttHHProducer);

  // -----   TPC Reco Sequence  --------------------------------------------
  PndTpcClusterFinderTask* tpcCF = new PndTpcClusterFinderTask();
  tpcCF->SetMode(1); // individual timeslice
  tpcCF->SetPersistence();
  tpcCF->timeslice(20); // = 4 sample times = 100ns @ 40MHz
  //tpcCF->SetTrivialClustering();
  fRun->AddTask(tpcCF);

  
  PndTpcIdealTrackingTask* tpcIPR = new PndTpcIdealTrackingTask();
  tpcIPR->useGeane(true);
  tpcIPR->useDistSorting(true);
  fRun->AddTask(tpcIPR);
  tpcIPR->SetPersistence();

  //------ Ideal DCH track finder --------------------
  PndDchFindTracks* finderTask = new PndDchFindTracks("dchFindTracks");
  finderTask->SetUseHitOrDigi("chit");
  fRun->AddTask(finderTask);
  // ------------------------------------------------- 
  PndDchTrackFinderIdealCylHit* mcTrackFinder = new  PndDchTrackFinderIdealCylHit();
  mcTrackFinder->SetPrimary(1);  // 1 = Only primary tracks are processed, 0 = all (default)
  finderTask->UseFinder(mcTrackFinder);
  //--------------------------------------------------
  PndDchMatchTracks *matchTask = new PndDchMatchTracks();//match PndDchTracks and MCTracks
  matchTask->SetUseHitOrDigi("chit");
  fRun->AddTask(matchTask);

  //----- MVD Hit Reco -----
  PndMvdClusterTask* mvdmccls = new PndMvdClusterTask();
  fRun->AddTask(mvdmccls);
        
  //------ GEM Realistic Track finder --------------------
  //Create and add finder task
  PndGemFindTracks* gemFinderTask = new  
      PndGemFindTracks("PndGemFindTracks");
  gemFinderTask->SetUseHitOrDigi("hit"); // hit = (default), digi
  fRun->AddTask(gemFinderTask);
        
  PndGemTrackFinderOnHits* gemTrackFinder = new   PndGemTrackFinderOnHits();
  gemTrackFinder->SetVerbose(0);  // verbosity level
  gemTrackFinder->SetPrimary(0);  // 1 = Only primary tracks are  processed, 0 = all (default)
  gemFinderTask->UseFinder(gemTrackFinder);
        
  PndGemTrackFinderQA* gemTrackFinderQA = new PndGemTrackFinderQA();
  gemTrackFinderQA->SetVerbose(0);
  fRun->AddTask(gemTrackFinderQA);

  KalmanTask* kalman =new KalmanTask();
  kalman->SetPersistence();
  kalman->SetNumIterations(3); // number of fitting iterations (back and forth)
  fRun->AddTask(kalman);

  TrackFitStatTask* fitstat=new TrackFitStatTask();
  fitstat->SetPersistence();
  fitstat->SetMCPCut(10); // in sigma dp/p
  fitstat->SetMCCuts(0.05, // pmin
                     10., // pmax
                     -TMath::Pi(),   // thetamin 5deg
                     TMath::Pi(),  // thetamax
                     5); // nPndTpcPoints
  //fitstat->SetPdgSelection(321);
  //fitstat->DoResiduals();
  fRun->AddTask(fitstat);
          
  // -----   Intialise and run   --------------------------------------------
  cout << "fRun->Init()" << endl;
  fRun->Init();

  timer.Start();
  fRun->Run(0,nEvents);
  // ------------------------------------------------------------------------


  // -----   Finish   -------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished successfully." << endl;
  cout << "Output file is "    << outFile << endl;
  cout << "Parameter file is " << parFile << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------
  
}

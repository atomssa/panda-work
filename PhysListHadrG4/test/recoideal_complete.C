void recoideal_complete()
{
  // Macro created 20/09/2006 by S.Spataro
  // It loads a simulation file and digitize hits for EMC

  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");  
  rootlogon();
  
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0; // just forget about it, for the moment
  
	// Number of events to process
  Int_t nEvents = 0;  // if 0 all the vents will be processed
  
  // Parameter file
  TString parFile = "simparams.root"; // at the moment you do not need it
  
  // Digitisation file (ascii)
  TString digiFile = "all.par";
  
  // Output file
  TString outFile = "reco_complete.root";
  
  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
    // ------------------------------------------------------------------------
  
  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile("sim_complete.root");
  fRun->AddFriend("digi_complete.root");
  fRun->SetOutputFile(outFile);
  fRun->SetWriteRunInfoFile(kFALSE);
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

  // ------------------------------------------------------------------------
  PndSttMvdGemTrackingIdeal* trackStt = new PndSttMvdGemTrackingIdeal();
  trackStt->SetRelativeMomentumSmearing(0.05);
  trackStt->SetVertexSmearing(0.05, 0.05, 0.05);
  trackStt->SetTrackingEfficiency(1.);
  trackStt->SetTrackOutput("SttMvdGemIdealTrack");
  trackStt->SetPersistence(kFALSE);
  fRun->AddTask(trackStt);
 
  /* 
  PndMCTrackAssociator* trackMC = new PndMCTrackAssociator();
  trackMC->SetTrackInBranchName("SttMvdGemIdealTrack");
  trackMC->SetTrackOutBranchName("SttMvdGemIdealTrackID");
  fRun->AddTask(trackMC);
  */

  PndRecoKalmanTask* recoKalman = new PndRecoKalmanTask();
  recoKalman->SetTrackInBranchName("SttMvdGemIdealTrack");
  //recoKalman->SetTrackInIDBranchName("SttMvdGemIdealTrackID");
  recoKalman->SetTrackOutBranchName("SttMvdGemGenTrack");
  recoKalman->SetBusyCut(50); // CHECK to be tuned
  //recoKalman->SetIdealHyp(kTRUE);
  //recoKalman->SetNumIterations(3);
  fRun->AddTask(recoKalman);

  PndMCTrackAssociator* trackMC2 = new PndMCTrackAssociator();
  trackMC2->SetTrackInBranchName("SttMvdGemGenTrack"); 
  trackMC2->SetTrackOutBranchName("SttMvdGemGenTrackID");
  fRun->AddTask(trackMC2);
 
  PndFtsTrackerIdeal* trackFts = new PndFtsTrackerIdeal();
  trackFts->SetRelativeMomentumSmearing(0.05);
  trackFts->SetVertexSmearing(0.05, 0.05, 0.05);
  trackFts->SetTrackingEfficiency(1.);
  trackFts->SetTrackOutput("FtsIdealTrack");
  trackFts->SetPersistence(kFALSE);
  fRun->AddTask(trackFts);

  PndMCTrackAssociator* trackMCfwd = new PndMCTrackAssociator();
  trackMCfwd->SetTrackInBranchName("FtsIdealTrack");
  trackMCfwd->SetTrackOutBranchName("FtsIdealTrackID");
  fRun->AddTask(trackMCfwd);

  PndRecoKalmanTask* recoKalmanFwd = new PndRecoKalmanTask();
  recoKalmanFwd->SetTrackInBranchName("FtsIdealTrack");
  //recoKalmanFwd->SetTrackInIDBranchName("FtsIdealTrackID");
  recoKalmanFwd->SetTrackOutBranchName("FtsIdealGenTrack");
  recoKalmanFwd->SetBusyCut(50); // CHECK to be tuned
  //recoKalmanFwd->SetIdealHyp(kTRUE);
  //recoKalmanFwd->SetNumIterations(3);
  fRun->AddTask(recoKalmanFwd);

  PndMCTrackAssociator* trackMC3 = new PndMCTrackAssociator();
  trackMC3->SetTrackInBranchName("FtsIdealGenTrack");
  trackMC3->SetTrackOutBranchName("FtsIdealGenTrackID");
  fRun->AddTask(trackMC3);

  // -----   Intialise and run   --------------------------------------------
  PndEmcMapper::Init(1);
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
  cout << " Test passed" << endl;
  cout << " All ok " << endl;
  exit(0);
}

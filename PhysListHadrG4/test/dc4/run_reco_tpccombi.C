{
  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0;

  // Input file
  TString inDigiFile = "digi_tpccombi.root";
  TString inSimFile = "points_tpccombi.root";

  // Parameter file
  TString parFile = "params_tpccombi.root";

  // Output file
  TString outFile = "reco_tpccombi.root";

  // Number of events to process
  Int_t nEvents = 0;
 
  // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  TString sysFile = gSystem->Getenv("VMCWORKDIR");
  // ------------------------------------------------------------------------
  // In general, the following parts need not be touched
  // ========================================================================

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------

  // -----   Digitization run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(inDigiFile);
  fRun->AddFriend(inSimFile);
  fRun->SetOutputFile(outFile);
  FairGeane *Geane = new FairGeane();
  fRun->AddTask(Geane);
  // ------------------------------------------------------------------------

  // -----  Parameter database   --------------------------------------------
   TString allDigiFile = sysFile+"/macro/params/all.par";

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(parFile.Data());
	
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(allDigiFile.Data(),"in");
        
  rtdb->setFirstInput(parInput1);
  rtdb->setSecondInput(parIo1);
  PndGeoHandling* geoH = PndGeoHandling::Instance();
  // ------------------------------------------------------------------------
  // -----   LHETRACK  ---------------------------------
  
  PndLheHitsMaker* trackMS = new PndLheHitsMaker("Tracking routine");
  trackMS->SetTpcMode(2);  // 0 OFF, 1 TpcPoint, 2 TpcCluster // TpcPoint smearing [cm], if negative no smearing
  trackMS->SetMvdMode(2);  // 0 OFF, 1 MVDPoint, 2 MVDHit     // MVDPoint smearing [cm], if negative no smearing
  trackMS->SetGemMode(2);  // 0 OFF, 1 GEMPoint, 2 GEMHit     // GEMPoint smearing [cm], if negative no smearing
  fRun->AddTask(trackMS);
  
  PndLheTrackFinder* trackFinder    = new PndLheTrackFinder();
  //PndLheTrackFinderIdeal* trackFinder    = new PndLheTrackFinderIdeal();
  fRun->AddTask(trackFinder);
  
  PndLheTrackFitter* trackFitter    = new PndLheTrackFitter("fitting");
  fRun->AddTask(trackFitter);
  
  PndRecoKalmanTask* recoKalman = new PndRecoKalmanTask();
  recoKalman->SetTrackInBranchName("LheTrack");
  recoKalman->SetTrackOutBranchName("LheGenTrack");
  //recoKalman->SetNumIterations(3);
  fRun->AddTask(recoKalman);
 
  // -----   Intialise and run   --------------------------------------------
  PndEmcMapper::Init(6);
  fRun->Init();
  fRun->Run(0, nEvents);

  rtdb->saveOutput();
  rtdb->print();

  // ------------------------------------------------------------------------

  // -----   Finish   -------------------------------------------------------

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file is "    << outFile << endl;
  cout << "Parameter file is " << parFile << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------


}

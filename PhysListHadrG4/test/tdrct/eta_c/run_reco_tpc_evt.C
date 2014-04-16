{
  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0;

  // Input file
  TString inDigiFile = "evt_digi_tpc.root";
  TString inSimFile = "evt_points_tpc.root";

  // Parameter file
  TString parFile = "evt_params_tpc.root";

  // Output file
  TString outFile = "evt_reco_tpc.root";

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
 
  // ------- RECO procedure ------------------------------------------------
 
  //correct for unfortunate shift in TPC digi
  //PndTpcRoughAlignmentTask* align = new PndTpcRoughAlignmentTask();
  //align->SetShift(TVector3(0.,0.,-3.71357e-01));   //old PSA
  //align->SetShift(TVector3(0.,0.,5.6E-2));      //new PSA
  //fRun->AddTask(align);

  //find PndTpcRiemannTracks in the TPC alone
  PndTpcRiemannTrackingTask* tpcSPR = new PndTpcRiemannTrackingTask();
  //tpcSPR->SetPersistence();
  //tpcSPR->SetVerbose(1);
  fRun->AddTask(tpcSPR);

  //build GFTracks from PndTpcRiemannTracks
  PndTpcTrackInitTask* trackInit=new PndTpcTrackInitTask();
  trackInit->SetPersistence();
  //trackInit->SetVerbose(1);
  trackInit->SetMCPid(); // use ideal particle identification
  //trackInit->SetPDG(211);
  //trackInit->useGeane(); // uses RKTrackrep and GeaneTrackrep
  trackInit->SetSmoothing(true);
  fRun->AddTask(trackInit);
  
  KalmanTask* kalman =new KalmanTask();
  kalman->SetPersistence();
  kalman->SetNumIterations(3); // number of fitting iterations (back and forth)
  fRun->AddTask(kalman);       // creates TrackPostFit branch

  //correlate fitted track with MVD pixels and strips
  PndTpcMVDCorrelatorTask* corr = new PndTpcMVDCorrelatorTask();
  corr->SetMatchDistance(0.18);   //cm
  corr->SetMinMVDHits(1);
  corr->SetOutTrackBranchName("TrackPreFitMVD");
  corr->SetPersistence(true);
  fRun->AddTask(corr);

  
  //fit after MVD corr
  KalmanTask* kalman2 =new KalmanTask();
  kalman2->SetPersistence();
  kalman2->SetNumIterations(3); // number of fitting iterations (back and forth)
  kalman2->SetTrackBranchName("TrackPreFitMVD");
  kalman2->SetOutBranchName("TrackPostFitMVD");
  fRun->AddTask(kalman2);

  PndTpcGEMCorrelatorTask* corrG = new PndTpcGEMCorrelatorTask();
  corrG->SetMatchDistance(0.5);   //cm
  corrG->SetMinGEMHits(2);
  corrG->SetTrackBranchName("TrackPostFitMVD");
  corrG->SetOutTrackBranchName("TrackPreFitGEM");
  corrG->SetPersistence(true);
  fRun->AddTask(corrG);

  //final fit
  KalmanTask* kalman3 =new KalmanTask();
  kalman3->SetPersistence();
  kalman3->SetNumIterations(3); // number of fitting iterations (back and forth)
  kalman3->SetTrackBranchName("TrackPreFitGEM");
  kalman3->SetOutBranchName("TrackPostFitComplete");
  fRun->AddTask(kalman3);

  PndGFTrackToPndTrackConvertorTask* converter =new PndGFTrackToPndTrackConvertorTask();
  converter->SetTrackInBranchName("TrackPostFitComplete");
  converter->SetTrackOutBranchName("PndTrackPostFitComplete");
  fRun->AddTask(converter);

  PndMCTrackAssociator* trackMC = new PndMCTrackAssociator();
  trackMC->SetTrackInBranchName("PndTrackPostFitComplete"); 
  trackMC->SetTrackOutBranchName("TrackPostFitCompleteID");
  fRun->AddTask(trackMC);

  // -----   Intialise and run   --------------------------------------------
  PndEmcMapper::Init(6);
  fRun->Init();
  fRun->Run(0, nEvents);

  rtdb->saveOutput();
  rtdb->print();

  corr->WriteHistograms("MVDRes.root");
  corrG->WriteHistograms("GEMRes.root");

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

{
  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0;

  // Input file (MC events)
  TString inFile = "points_tpccombi.root";

  // Parameter file
  TString parFile = "params_tpccombi.root";

  // Output file
  TString outFile = "digi_tpccombi.root";

  // Number of events to process
  Int_t nEvents = 0;
 
  TString mcMode = "TGeant3";
  // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  TString sysFile = gSystem->Getenv("VMCWORKDIR");
  // ------------------------------------------------------------------------

  // ---  Now choose concrete engines for the different tasks   -------------
  // ------------------------------------------------------------------------

  // In general, the following parts need not be touched
  // ========================================================================

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------

  // -----   Digitization run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(inFile);
  fRun->SetOutputFile(outFile);
  
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
  // ------------------------------------------------------------------------

  // -----   TPC digi producers   ---------------------------------
  PndTpcClusterizerTask* tpcClusterizer = new PndTpcClusterizerTask();
   if(mcMode=="TGeant3") tpcClusterizer->SetMereChargeConversion();
  //tpcClusterizer->SetPersistence();
  fRun->AddTask(tpcClusterizer);
 
  PndTpcDriftTask* tpcDrifter = new PndTpcDriftTask();
  // tpcDrifter->SetPersistence();
  tpcDrifter->SetDistort(false);
  fRun->AddTask(tpcDrifter);

  PndTpcGemTask* tpcGem = new PndTpcGemTask();
  //tpcGem->SetPersistence();
  fRun->AddTask(tpcGem);

  PndTpcPadResponseTask* tpcPadResponse = new PndTpcPadResponseTask();
  //tpcPadResponse->SetPersistence();
  fRun->AddTask(tpcPadResponse);

  PndTpcElectronicsTask* tpcElec = new PndTpcElectronicsTask();
  tpcElec->SetPersistence();
  fRun->AddTask(tpcElec);

  PndTpcClusterFinderTask* tpcCF = new PndTpcClusterFinderTask();
  tpcCF->SetPersistence();
  tpcCF->timeslice(20); // = 4 sample times = 100ns @ 40MHz
  //
  // Following three lines added on request by Felix, 20/02/2011
  //
  tpcCF->SetMode(1);
  tpcCF->SetDiffFactor(1.3);
  tpcCF->SetErrorPars(600,300);
  //
  fRun->AddTask(tpcCF);

  // -----   MDV digi producers   --------------------------------- 
  PndMvdDigiTask* mvddigi = new PndMvdDigiTask();
  mvddigi->SetVerbose(iVerbose);
  fRun->AddTask(mvddigi);

  PndMvdClusterTask* mvdmccls = new PndMvdClusterTask();
  mvdmccls->SetVerbose(iVerbose);
  fRun->AddTask(mvdmccls); 
  // -----   EMC hit producers   ---------------------------------
  //PndEmcHitProducer* emcHitProd = new PndEmcHitProducer();
  //fRun->AddTask(emcHitProd); // hit production 

  //PndEmcMakeDigi* emcMakeDigi=new PndEmcMakeDigi();
  //fRun->AddTask(emcMakeDigi); // fast digitization

  PndEmcHitsToWaveform* emcHitsToWaveform= new PndEmcHitsToWaveform(iVerbose);
  PndEmcWaveformToDigi* emcWaveformToDigi=new PndEmcWaveformToDigi(iVerbose);
  emcHitsToWaveform->SetStorageOfData(kFALSE);
  emcWaveformToDigi->SetStorageOfData(kFALSE);
  fRun->AddTask(emcHitsToWaveform);  // full digitization
  fRun->AddTask(emcWaveformToDigi);  // full digitization

  PndEmcMakeCluster* emcMakeCluster= new PndEmcMakeCluster(iVerbose);
  fRun->AddTask(emcMakeCluster);

  PndEmcMakeBump* emcMakeBump= new PndEmcMakeBump();
  fRun->AddTask(emcMakeBump);

  PndEmcHdrFiller* emcHdrFiller = new PndEmcHdrFiller();
  fRun->AddTask(emcHdrFiller); // ECM header
  
  // -----   MDT hit producers   ---------------------------------
  PndMdtHitProducerIdeal* mdtHitProd = new PndMdtHitProducerIdeal();
  mdtHitProd->SetPositionSmearing(.3); // position smearing [cm]
  fRun->AddTask(mdtHitProd);
  
  PndMdtTrkProducer* mdtTrkProd = new PndMdtTrkProducer();
  fRun->AddTask(mdtTrkProd);

  // -----   DRC hit producers   ---------------------------------
  PndDrcHitProducerIdeal* drchit = new PndDrcHitProducerIdeal();
  drchit->SetVerbose(iVerbose);
  fRun->AddTask(drchit);
  
  // -----   GEM hit producers   ---------------------------------
  Int_t verboseLevel = 0;
  PndGemDigitize* gemDigitize = new PndGemDigitize("GEM Digitizer", verboseLevel);
  fRun->AddTask(gemDigitize);

  PndGemFindHits* gemFindHits = new PndGemFindHits("GEM Hit Finder", verboseLevel);
  fRun->AddTask(gemFindHits);

  // -----   Intialise and run   --------------------------------------------
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

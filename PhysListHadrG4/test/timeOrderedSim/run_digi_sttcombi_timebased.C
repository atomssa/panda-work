{
  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0;

  TString inFile = "Mvd_Sim.root";



  // Number of events to process
  Int_t nEvents = 100;
 
  // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  TString sysFile = gSystem->Getenv("VMCWORKDIR");
  // ------------------------------------------------------------------------

  // Output file
  PndFileNameCreator creator(inFile.Data());
  TString parFile = creator.GetParFileName().c_str();
  TString outFile = creator.GetDigiFileName("timebased").c_str();
  std::cout << "DigiFileName: " << outFile.Data() << std::endl;

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
  fRun->SetEventMeanTime(50);
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

   // -----   STT digi producers   --------------------------------- 
  PndSttHitProducerRealFull* sttHitProducer = new PndSttHitProducerRealFull();
  sttHitProducer->RunTimeBased();
  fRun->AddTask(sttHitProducer);

  PndSttHitSorterTask* sttSorter = new PndSttHitSorterTask(5000, 50, "STTHit", "STTSortedHits", "PndSTT");
  fRun->AddTask(sttSorter);
 
  // -----   MDV digi producers   --------------------------------- 
  PndMvdDigiTask* mvddigi = new PndMvdDigiTask();
  mvddigi->RunTimeBased();
  mvddigi->SetVerbose(iVerbose);
  fRun->AddTask(mvddigi);


  // -----   EMC hit producers   ---------------------------------
  PndEmcHitsToWaveform* emcHitsToWaveform= new PndEmcHitsToWaveform(iVerbose);
  PndEmcWaveformToDigi* emcWaveformToDigi=new PndEmcWaveformToDigi(iVerbose);
  emcHitsToWaveform->SetStorageOfData(kTRUE);
  //emcWaveformToDigi->SetStorageOfData(kFALSE);
  fRun->AddTask(emcHitsToWaveform);  // full digitization
  fRun->AddTask(emcWaveformToDigi);  // full digitization

  // -----   SciT hit producers   ---------------------------------
  PndSciTHitProducerIdeal* tofhit = new PndSciTHitProducerIdeal();
  tofhit->SetVerbose(iVerbose);
  fRun->AddTask(tofhit);

  // -----   MDT hit producers   ---------------------------------
  PndMdtHitProducerIdeal* mdtHitProd = new PndMdtHitProducerIdeal();
  mdtHitProd->SetPositionSmearing(.3); // position smearing [cm]
  fRun->AddTask(mdtHitProd);
  
  // -----   DRC hit producers   ---------------------------------
  PndDrcHitProducerIdeal* drchit = new PndDrcHitProducerIdeal();
  drchit->SetVerbose(iVerbose);
  fRun->AddTask(drchit);
  // -----   GEM hit producers   ---------------------------------
  Int_t verboseLevel = 0;
  PndGemDigitize* gemDigitize = new PndGemDigitize("GEM Digitizer", verboseLevel);
  fRun->AddTask(gemDigitize);

  // -----   FTS hit producers   ---------------------------------
  PndFtsHitProducerRealFull* ftsHitProducer = new PndFtsHitProducerRealFull();
  ftsHitProducer->RunTimeBased();
  //PndFtsHitProducerIdeal* ftsHitProducer = new PndFtsHitProducerIdeal();
  //PndFtsHitProducerRealFull* ftsHitProducer = new PndFtsHitProducerRealFull();
  fRun->AddTask(ftsHitProducer);

  PndFtsHitSorterTask* ftsSorter = new PndFtsHitSorterTask(5000, 50, "FTSHit", "FTSSortedHits", "PndFTS");
//  fRun->AddTask(ftsSorter);
  // -----   Ftof hit producers   ---------------------------------
  PndFtofHitProducerIdeal* ftofhit = new PndFtofHitProducerIdeal();
  ftofhit->SetVerbose(iVerbose);
  fRun->AddTask(ftofhit);
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

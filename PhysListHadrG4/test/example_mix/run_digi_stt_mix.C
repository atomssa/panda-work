{
  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0;

  // Parameter file
  TString parFile = "simparams_stt.root";

  // Output file
  TString outFile = "mixed_digi_stt.root";

  // Number of events to process
  Int_t nEvents = 120;
 
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
  fRun->SetOutputFile(outFile);
  
    //----- Set log level   ----------------------------------------
  
    //FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
  
  
    //----- Set Input files ----------------------------------------
  
    //** Set BG file */
  fRun->SetBackgroundFile("sim_stt_bg.root");
    //** Set first signal file */
  fRun->SetSignalFile("sim_stt_s1.root",1);
    //** Set second signal file */
  fRun->SetSignalFile("sim_stt_s2.root",2);
  
  /** Chained files can be added using 
   fRun->AddBackgroundFile(TString name ) or for signal files 
   fRun->AddSignalFile(TString name, UInt_t identifier )  identifier=1 will add to chain 1
   
   */
  
  
    //----- Mix using entries  ----------------------------------------
  
  /** for each ~20 entries background 1 entry from signal chain  1 will be read  */
    //fRun->BGWindowWidthNo(20,1);
  /** for each ~30 entries background 1 entry from signal chain  2 will be read  */
    //fRun->BGWindowWidthNo(30,2);
  
  
    //----- Mix using time       ----------------------------------------
  
  /**Set the event mean time, event time will be a random number generated from (1/T)exp(-x/T) */
  fRun->SetEventMeanTime(10);
  
  /** each ~100 ns background 1 entry from signal chain  1 will be read  */
  fRun->BGWindowWidthTime(100,1);
  /** each ~60 ns background 1 entry from signal chain  2 will be read  */
  fRun->BGWindowWidthTime(60,2);
  
  
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
  PndSttHitProducerRealFast* sttHitProducer = new PndSttHitProducerRealFast();
  fRun->AddTask(sttHitProducer);
 
  // -----   MDV digi producers   --------------------------------- 
  PndMvdDigiTask* mvddigi = new PndMvdDigiTask();
  mvddigi->SetVerbose(iVerbose);
  fRun->AddTask(mvddigi);

  PndMvdClusterTask* mvdmccls = new PndMvdClusterTask();
  mvdmccls->SetVerbose(iVerbose);
  fRun->AddTask(mvdmccls); 

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

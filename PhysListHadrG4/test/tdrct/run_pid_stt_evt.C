{
  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0;
  Int_t nEvents = 0;
  // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  TString sysFile = gSystem->Getenv("VMCWORKDIR");
  // ------------------------------------------------------------------------
  // Output file
  TString parFile = "evt_params_stt.root";
  TString inSimuFile = "evt_points_stt.root";
  TString inDigiFile = "evt_digi_stt.root";
  TString inRecoFile = "evt_reco_stt.root";

  TString outFile = "evt_pid_stt.root";
   
  // In general, the following parts need not be touched
  // ========================================================================

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------
  
  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(inSimuFile);
  fRun->AddFriend(inDigiFile);
  fRun->AddFriend(inRecoFile);
  fRun->SetOutputFile(outFile.Data());
  FairGeane *Geane = new FairGeane();
  fRun->AddTask(Geane);
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
  
  PndPidCorrelator* corr = new PndPidCorrelator();
  //corr->SetVerbose();
  corr->SetInputBranch("SttMvdGemGenTrack");
  corr->SetInputIDBranch("SttMvdGemGenTrackID");
  //corr->SetDebugMode(kTRUE);
  fRun->AddTask(corr);
 
  PndPidIdealAssociatorTask *assMC= new PndPidIdealAssociatorTask();
  fRun->AddTask(assMC);

  PndPidMvdAssociatorTask *assMvd= new PndPidMvdAssociatorTask();
  fRun->AddTask(assMvd);

  PndPidMdtHCAssociatorTask *assMdt= new PndPidMdtHCAssociatorTask();
  fRun->AddTask(assMdt);

  PndPidDrcAssociatorTask *assDrc= new PndPidDrcAssociatorTask();
  fRun->AddTask(assDrc);

  PndPidDiscAssociatorTask *assDisc= new PndPidDiscAssociatorTask();
  fRun->AddTask(assDisc);

  PndPidSttAssociatorTask *assStt= new PndPidSttAssociatorTask();
  fRun->AddTask(assStt);
 
  // -----   Intialise and run   --------------------------------------------
  PndEmcMapper::Init(6);
  fRun->Init();
  fRun->Run(0,nEvents);
  // ------------------------------------------------------------------------
  rtdb->print();
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

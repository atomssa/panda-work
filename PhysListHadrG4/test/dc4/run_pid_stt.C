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
  TString parFile = "params_sttcombi.root";
  TString inSimuFile = "points_sttcombi.root";
  TString inDigiFile = "digi_sttcombi.root";
  TString inRecoFile = "reco_sttcombi.root";

  TString outFile = "pid_sttcombi.root";
   
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
  corr->SetInputBranch("SttMvdGenTrack");
  //corr->SetInputIDBranch("LheTrackID");
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

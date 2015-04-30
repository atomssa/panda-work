void brem(int ibinsize=0 ) {

  TString dir = gSystem->Getenv("ODIR");
  // Macro created 02/10/2012 by S.Spataro
  // It loads a reconstruction file and compute PID informations

  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0; // just forget about it, for the moment

  // Number of events to process
  Int_t nEvents = 0;  // if 0 all the vents will be processed

  // Parameter file
  TString parFile = dir+"/simparams.root"; // at the moment you do not need it

  // Digitisation file (ascii)
  TString digiFile = "all.par";

  // Output file
  TString outFile = dir+"/brem_complete.root";

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  // ------------------------------------------------------------------------

  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(dir+"/pid_complete.root");
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

  //PndPidBremCorrectorNT *bremCorr = new PndPidBremCorrectorNT();
  //fRun->AddTask(bremCorr);
  BremPidReader *bpr = new BremPidReader(ibinsize);
  bpr->set_output_name(dir+Form("/bremcorr.ibs.%d.root",ibinsize));
  fRun->AddTask(bpr);

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

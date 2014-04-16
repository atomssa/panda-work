void anaPatternRecoQuality(Int_t nEvents = 0)
{
  // Macro created 20/09/2006 by S.Spataro
  // It loads a digi file and performs tracking

  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0; // just forget about it, for the moment

	// Number of events to process
//  Int_t nEvents = 0;  // if 0 all the vents will be processed

  // Parameter file
  TString parFile = "simparams.root"; // at the moment you do not need it

  // Digitisation file (ascii)
  TString digiFile = "all.par";

  // Output file
  TString outFile = "ana_complete.root";

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
    // ------------------------------------------------------------------------

  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile("sim_complete.root");
  fRun->AddFriend("digi_complete.root");
  fRun->AddFriend("reco_complete.root");
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

  PndMCMatchCreatorTask* mcMatch = new PndMCMatchCreatorTask();
  mcMatch->SetVerbose(0);
  fRun->AddTask(mcMatch);

//  PndMCTestPatternRecoQuality* mcCheck = new PndMCTestPatternRecoQuality("SttMvdGemGenTrack");
//  fRun->AddTask(mcCheck);

  PndTrackingQualityTask* mcCheck = new PndTrackingQualityTask("SttMvdGemGenTrack", kTRUE);
  fRun->AddTask(mcCheck);

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
//  exit(0);
}







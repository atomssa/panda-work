///-- E.A. Run with the suggested temporary fix of putting MVD module addition at the end after all detectors

// ج

// Macro for running Panda simulation  with Geant3  or Geant4 (M. Al-Turany)
// This macro is supposed to run the full simulation of the panda detector
// to run the macro:
// root  sim_complete.C  or in root session root>.x  sim_complete_stt.C
// to run with different options:(e.g more events, different momentum, Geant4)
// root  sim_complete.C"(100, "TGeant4",2)"

sim_complete_prod(Int_t batch = 0, TString SimEngine="TGeant4" /* or TGeant3 */ )
{

  Int_t nEvents = 200000;
  //-----User Settings:-----------------------------------------------
  TString  OutputFile     =Form("output/sim_complete_%s_%d.root",SimEngine.Data(), batch);
  TString  ParOutputfile  =Form("output/simparams_%s_%d.root",SimEngine.Data(),batch);
  TString  MediaFile      ="media_pnd.geo";
  gDebug                  = 0;
  TString digiFile        = "all.par"; //The emc run the hit producer directly 

  double BeamMomentum   =15.0; // ** change HERE if you run Box generator

  //TString  SimEngine ="TGeant3";
  
  //------------------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  gRandom->SetSeed(); 

  // Create the Simulation run manager--------------------------------
  FairRunSim *fRun = new FairRunSim();
  fRun->SetName(SimEngine.Data() );
  fRun->SetOutputFile(OutputFile.Data());
  fRun->SetWriteRunInfoFile(kFALSE);
  fRun->SetBeamMom(BeamMomentum);
  fRun->SetMaterials(MediaFile.Data());

  fRun->SetRadLenRegister(true);
  //fRun->SetRadMapRegister(true);
  //fRun->SetRadGridRegister(true);

  FairRuntimeDb *rtdb=fRun->GetRuntimeDb();
  
  // Set the parameters 
  //-------------------------------
  TString allDigiFile = gSystem->Getenv("VMCWORKDIR");
  allDigiFile += "/macro/params/";
  allDigiFile += digiFile;
 
  //-------Set the parameter output --------------------
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(allDigiFile.Data(),"in");
  rtdb->setFirstInput(parIo1);        

  //---------------------Set Parameter output      ---------- 
  Bool_t kParameterMerged=kTRUE;
  FairParRootFileIo* output=new FairParRootFileIo(kParameterMerged);
  output->open(ParOutputfile.Data());
  rtdb->setOutput(output);

  // Create and add detectors

  //-------------------------  CAVE      ----------------- 
  FairModule *Cave= new PndCave("CAVE");
  Cave->SetGeometryFileName("pndcave.geo");
  fRun->AddModule(Cave); 
  
  //-------------------------  Magnet   ----------------- 
  FairModule *Magnet= new PndMagnet("MAGNET");
  //Magnet->SetGeometryFileName("FullSolenoid_V842.root");
  Magnet->SetGeometryFileName("FullSuperconductingSolenoid_v831.root");
  fRun->AddModule(Magnet);
  FairModule *Dipole= new PndMagnet("MAGNET");
  Dipole->SetGeometryFileName("dipole.geo");
  fRun->AddModule(Dipole);
  
  //-------------------------  Pipe     -----------------
  FairModule *Pipe= new PndPipe("PIPE");
  Pipe->SetGeometryFileName("beampipe_201309.root");
  fRun->AddModule(Pipe);

  //-------------------------  STT       -----------------
  FairDetector *Stt= new PndStt("STT", kTRUE);
  Stt->SetGeometryFileName("straws_skewed_blocks_35cm_pipe.geo");
  fRun->AddModule(Stt);
  
  //-------------------------  GEM       -----------------
  FairDetector *Gem = new PndGemDetector("GEM", kTRUE);
  Gem->SetGeometryFileName("gem_3Stations.root");
  fRun->AddModule(Gem);
  
  //-------------------------  EMC       -----------------
  PndEmc *Emc = new PndEmc("EMC",kTRUE);
  Emc->SetGeometryVersion(1);
  Emc->SetStorageOfData(kFALSE);
  fRun->AddModule(Emc);

  //-------------------------  MVD       -----------------
  FairDetector *Mvd = new PndMvdDetector("MVD", kTRUE);
  Mvd->SetGeometryFileName("Mvd-2.1_FullVersion.root");
  fRun->AddModule(Mvd);
  
  //-------------------------  SCITIL    -----------------
  FairDetector *SciT = new PndSciT("SCIT",kTRUE);
  SciT->SetGeometryFileName("barrel-SciTil_07022013.root");
  fRun->AddModule(SciT);
  
  //-------------------------  DRC       -----------------
  PndDrc *Drc = new PndDrc("DIRC", kTRUE);
  Drc->SetGeometryFileName("dirc_l0_p0_updated.root"); 
  Drc->SetRunCherenkov(kFALSE);
  fRun->AddModule(Drc); 
  
  //-------------------------  DISC      -----------------
  PndDsk* Dsk = new PndDsk("DSK", kTRUE);
  Dsk->SetStoreCerenkovs(kFALSE);
  Dsk->SetStoreTrackPoints(kFALSE);
  fRun->AddModule(Dsk);
  
  //-------------------------  MDT       -----------------
  PndMdt *Muo = new PndMdt("MDT",kTRUE);
  Muo->SetBarrel("fast");
  Muo->SetEndcap("fast");
  Muo->SetMuonFilter("fast");
  Muo->SetForward("fast");
  Muo->SetMdtMagnet(kTRUE);
  Muo->SetMdtMFIron(kTRUE);
  fRun->AddModule(Muo);
  
  //-------------------------  FTS       -----------------
  FairDetector *Fts= new PndFts("FTS", kTRUE);
  Fts->SetGeometryFileName("fts.geo");
  fRun->AddModule(Fts); 
  
  //-------------------------  FTOF      -----------------
  FairDetector *FTof = new PndFtof("FTOF",kTRUE);
  FTof->SetGeometryFileName("ftofwall.root");
  fRun->AddModule(FTof);
  
  //-------------------------  RICH       ----------------
  FairDetector *Rich= new PndRich("RICH",kFALSE);
  Rich->SetGeometryFileName("rich_v2_shift.geo");
  fRun->AddModule(Rich);

  
  // Create and Set Event Generator
  //-------------------------------
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  fRun->SetGenerator(primGen);

  FairBoxGenerator* boxGen = new FairBoxGenerator(0, 1);
  boxGen->SetPRange(2.0,2.0); // GeV/c
  boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
  boxGen->SetThetaRange(0., 180.); // Polar angle in lab system range [degree]
  boxGen->SetXYZ(0., 0., 0.); // cm
  primGen->AddGenerator(boxGen);

  //---------------------Create and Set the Field(s)---------- 
  PndMultiField *fField= new PndMultiField("FULL");
  fRun->SetField(fField);

  // EMC Hit producer
  //-------------------------------
  PndEmcHitProducer* emcHitProd = new PndEmcHitProducer();
  fRun->AddTask(emcHitProd);
  
  //-------------------------  Initialize the RUN  -----------------  
  fRun->Init();
  //-------------------------  Run the Simulation  -----------------   
  fRun->Run(nEvents);
  //-------------------------  Save the parameters ----------------- 
  rtdb->saveOutput();
  //------------------------Print some info and exit----------------     
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
  
  cout << " Test passed" << endl;
  cout << " All ok " << endl;
  
  //exit(0);

}  
  

// Macro for running Panda simulation  with Geant3  or Geant4 (M. Al-Turany)
// This macro is supposed to run the full simulation of the panda detector and
// to store data for the calculation of radiation lenght
// to run the macro:
// root  rad_complete.C  or in root session root>.x  rad_complete.C
// to run with different options:(e.g more events, different momentum, Geant4)
// root  rad_complete.C"(100, "TGeant4",2)"

rad_complete(Int_t nEvents = 10, TString  SimEngine ="TGeant3", Float_t mom = 7.24)
{
  //-----User Settings:-----------------------------------------------
  TString  OutputFile     ="sim_complete.root";
  TString  ParOutputfile  ="simparams.root";
  Double_t BeamMomentum   =15.0; // beam momentum ONLY for the scaling of the dipole field. For the generator use "mom"
  TString  MediaFile      ="media_pnd.geo";
  gDebug                  = 0;
  
  //------------------------------------------------------------------
  
  TStopwatch timer;
  timer.Start();
 
  // Load basic libraries---------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  gRandom->SetSeed(); 
  // Create the Simulation run manager--------------------------------
  FairRunSim *fRun = new FairRunSim();
  fRun->SetName(SimEngine.Data() );
  fRun->SetOutputFile(OutputFile.Data());
  fRun->SetBeamMom(BeamMomentum);
  fRun->SetMaterials(MediaFile.Data());
  FairRuntimeDb *rtdb=fRun->GetRuntimeDb();
  
  //---------------------Set Parameter output      ---------- 
  Bool_t kParameterMerged=kTRUE;
  FairParRootFileIo* output=new FairParRootFileIo(kParameterMerged);
  output->open(ParOutputfile.Data());
  rtdb->setOutput(output);

 fRun->SetRadLenRegister(kTRUE);

   // Create and add detectors

 //-------------------------  CAVE      -----------------

  FairModule *Cave= new PndCave("CAVE");
  Cave->SetGeometryFileName("pndcave.geo");
  fRun->AddModule(Cave); 
  //-------------------------  Magnet   ----------------- 
  /*
  FairModule *Magnet= new PndMagnet("MAGNET");
  //Magnet->SetGeometryFileName("FullSolenoid_V842.root");
  Magnet->SetGeometryFileName("FullSuperconductingSolenoid_v831.root");
  fRun->AddModule(Magnet);
  FairModule *Dipole= new PndMagnet("MAGNET");
  Dipole->SetGeometryFileName("dipole.geo");
  fRun->AddModule(Dipole);
  */
  //-------------------------  Pipe     -----------------
  FairModule *Pipe= new PndPipe("PIPE");
  Pipe->SetGeometryFileName("beampipe_201112.root");
  fRun->AddModule(Pipe);
  //-------------------------  STT       -----------------
  FairDetector *Stt= new PndStt("STT", kFALSE);
  Stt->SetGeometryFileName("straws_skewed_blocks_35cm_pipe.geo");
  fRun->AddModule(Stt);
  //-------------------------  MVD       -----------------
  FairDetector *Mvd = new PndMvdDetector("MVD", kFALSE);
  Mvd->SetGeometryFileName("Mvd-2.1_FullVersion.root");
  fRun->AddModule(Mvd);
  //-------------------------  GEM       -----------------
  FairDetector *Gem = new PndGemDetector("GEM", kFALSE);
  Gem->SetGeometryFileName("gem_3Stations.root");
  fRun->AddModule(Gem);
  //-------------------------  EMC       -----------------
  PndEmc *Emc = new PndEmc("EMC",kFALSE);
  Emc->SetGeometryVersion(1);
  Emc->SetStorageOfData(kFALSE);
  fRun->AddModule(Emc);
  //-------------------------  SCITIL    -----------------
  FairDetector *SciT = new PndSciT("SCIT",kFALSE);
  SciT->SetGeometryFileName("barrel-SciTil_18122012.root");
  fRun->AddModule(SciT);
  //-------------------------  DRC       -----------------
  PndDrc *Drc = new PndDrc("DIRC", kFALSE);
  Drc->SetGeometryFileName("dirc_l0_p0_updated.root"); 
  Drc->SetRunCherenkov(kFALSE);
  fRun->AddModule(Drc); 
  //-------------------------  DISC      -----------------
  PndDsk* Dsk = new PndDsk("DSK", kFALSE);
  Dsk->SetStoreCerenkovs(kFALSE);
  Dsk->SetStoreTrackPoints(kFALSE);
  fRun->AddModule(Dsk);
  //-------------------------  MDT       -----------------
  PndMdt *Muo = new PndMdt("MDT",kFALSE);
  Muo->SetBarrel("fast");
  Muo->SetEndcap("fast");
  Muo->SetMuonFilter("fast");
  Muo->SetForward("fast");
  Muo->SetMdtMagnet(kTRUE);
  Muo->SetMdtMFIron(kTRUE);
  fRun->AddModule(Muo);
  //-------------------------  FTS       -----------------
  FairDetector *Fts= new PndFts("FTS", kFALSE);
  Fts->SetGeometryFileName("fts.geo");
  fRun->AddModule(Fts); 
  //-------------------------  FTOF      -----------------
  FairDetector *FTof = new PndFtof("FTOF",kFALSE);
  FTof->SetGeometryFileName("ftofwall.root");
  fRun->AddModule(FTof);
  //-------------------------  RICH       ----------------
  FairDetector *Rich= new PndRich("RICH",kFALSE);
  Rich->SetGeometryFileName("rich_v2.geo");
  fRun->AddModule(Rich);

  // Create and Set Event Generator
  //-------------------------------
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  fRun->SetGenerator(primGen);
	 
  FairBoxGenerator* boxGen = new FairBoxGenerator(0, 10); // 0 = geantino; 10 = multipl.
  boxGen->SetPtRange(mom,mom); // GeV/c
  boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
  boxGen->SetThetaRange(0., 90.); // Polar angle in lab system range [degree]
  boxGen->SetXYZ(0., 0., 0.); // mm o cm ??
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
  
  exit(0);

}  
  

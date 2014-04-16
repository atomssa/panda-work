// Macro created 03/05/2011 by S.Spataro
// It creates a DPM background simulation for the tracking TDR
run_sim_stt_dpm(Int_t nEvents=10, Float_t mom = 4.0, Int_t mode =1, UInt_t seed=0)
{
  gRandom->SetSeed(seed);

  TStopwatch timer;
  timer.Start();
  gDebug=0;
  // Load basic libraries
  // If it does not work,  please check the path of the libs and put it by hands
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  
  TString digiFile = "all.par";
  TString parFile = "dpm_params_stt.root";
  
  FairRunSim *fRun = new FairRunSim();

  // set the MC version used
  // ------------------------

  fRun->SetName("TGeant3");
  //fRun->SetName("TGeant4");

  fRun->SetOutputFile("dpm_points_stt.root");

  // Set the parameters
  //-------------------------------
  TString allDigiFile = gSystem->Getenv("VMCWORKDIR");
  allDigiFile += "/macro/params/";
  allDigiFile += digiFile;
 
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(allDigiFile.Data(),"in");
  rtdb->setFirstInput(parIo1);        
  Bool_t kParameterMerged=kTRUE;
	
  FairParRootFileIo* output=new FairParRootFileIo(kParameterMerged);
  output->open(parFile);
  rtdb->setOutput(output);
  
  // Set Material file Name
  //-----------------------
  fRun->SetMaterials("media_pnd.geo");

  // Create and add detectors
  //-------------------------
  FairModule *Cave= new PndCave("CAVE");
  Cave->SetGeometryFileName("pndcave.geo");
  fRun->AddModule(Cave); 

  FairModule *Magnet= new PndMagnet("MAGNET");
  //Magnet->SetGeometryFileName("FullSolenoid_V842.root");
  Magnet->SetGeometryFileName("FullSuperconductingSolenoid_v831.root");
  fRun->AddModule(Magnet);

  FairModule *Pipe= new PndPipe("PIPE");
  fRun->AddModule(Pipe);

  FairDetector *Stt= new PndStt("STT", kTRUE);
  Stt->SetGeometryFileName("straws_skewed_blocks_35cm_pipe.geo");
  fRun->AddModule(Stt);

  FairDetector *Mvd = new PndMvdDetector("MVD", kTRUE);
  Mvd->SetGeometryFileName("Mvd-2.1_FullVersion.root");
  fRun->AddModule(Mvd);

  PndEmc *Emc = new PndEmc("EMC",kFALSE);
  Emc->SetGeometryVersion(19); 
  Emc->SetStorageOfData(kFALSE);
  fRun->AddModule(Emc);
/*
  PndMdt *Muo = new PndMdt("MDT",kFALSE);
  Muo->SetBarrel("fast");
  Muo->SetEndcap("fast");
  Muo->SetMuonFilter("fast");
  Muo->SetMdtMagnet(kTRUE);
  Muo->SetMdtMFIron(kTRUE);
  fRun->AddModule(Muo);
*/
  FairDetector *Gem = new PndGemDetector("GEM", kFALSE);
  Gem->SetGeometryFileName("gem_3Stations.root");
  fRun->AddModule(Gem);

  PndDsk* Dsk = new PndDsk("DSK", kFALSE);
  Dsk->SetGeometryFileName("dsk.root");
  Dsk->SetStoreCerenkovs(kFALSE);
  Dsk->SetStoreTrackPoints(kFALSE);
  fRun->AddModule(Dsk);

  PndDrc *Drc = new PndDrc("DIRC", kFALSE);
  Drc->SetGeometryFileName("dirc_l0_p0.root");
  Drc->SetRunCherenkov(kFALSE); // for fast sim Cherenkov -> kFALSE
  fRun->AddModule(Drc);
  
  // Create and Set Event Generator
  //-------------------------------

  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  primGen->SetTarget(0., 0.5/2.355);
  primGen->SmearVertexZ(kTRUE);
  primGen->SmearGausVertexZ(kTRUE);
  primGen->SetBeam(0., 0., 0.1, 0.1);
  primGen->SmearVertexXY(kTRUE);
  fRun->SetGenerator(primGen);

  PndDpmDirect *dpmGen = new PndDpmDirect(mom,mode);
  primGen->AddGenerator(dpmGen);

  // Create and Set Magnetic Field
  //-------------------------------
  fRun->SetBeamMom(mom);
  PndMultiField *fField= new PndMultiField("FULL");
  fRun->SetField(fField);
/*
  // -----   STT digi producers   --------------------------------- 
  PndSttHitProducerRealFast* sttHitProducer = new PndSttHitProducerRealFast();
  fRun->AddTask(sttHitProducer);
 
  // -----   MDV digi producers   --------------------------------- 
  PndMvdDigiTask* mvddigi = new PndMvdDigiTask();
  mvddigi->SetVerbose(0);
  fRun->AddTask(mvddigi);

  PndMvdClusterTask* mvdmccls = new PndMvdClusterTask();
  mvdmccls->SetVerbose(0);
  fRun->AddTask(mvdmccls); 
  */  
  /**Initialize the session*/
  fRun->Init();
  
  rtdb->setOutput(output);
  rtdb->saveOutput();
  rtdb->print();

  // Transport nEvents
  // -----------------
  fRun->Run(nEvents);

  timer.Stop();

  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

}

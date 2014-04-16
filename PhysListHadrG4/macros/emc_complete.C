// Macro created 20/09/2006 by S.Spataro
// It creates a geant simulation file for emc
/**
 *  Modified 17.12.08 by M. Al-Turany
 *  The hit producing is added to the simulation session, so that the EmcHits are produced on the fly 
 *  during the simulation
 */

emc_complete(Int_t nEvents = 10, Float_t mom = 1., Int_t charge = 1,
	     TString phys_list, Bool_t full_panda,
	     TString out_dat, TString out_par){

  TStopwatch timer;
  timer.Start();
  gDebug=0;

  // Load basic libraries
  // If it does not work,  please check the path of the libs and put it by hands
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
  rootlogon();
  basiclibs();
  //gSystem->ListLibraries();
  
  FairRunSim *fRun = new FairRunSim();
	
  // set the MC version used
  // ------------------------
  Bool_t G3 = strncmp(phys_list.Data(),"G3_",3)==0;
  cout << "Setting up MC engine to " << (G3?"TGeant3":"TGeant4") << " with " << (full_panda?"full PANDA":"EMCal only")<< endl;
  fRun->SetName(G3?"TGeant3":"TGeant4");

  fRun->SetOutputFile(out_dat);
	
  /**Get the run time data base for this session and set the needed input*/
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();

  /**Set the digitization parameters */
  TString emcDigiFile = gSystem->Getenv("VMCWORKDIR");
  emcDigiFile += "/macro/params/emc.par";
  
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(emcDigiFile.Data(),"in");
  rtdb->setFirstInput(parIo1);        

  /**Parameters created for this simulation goes to the out put*/
  Bool_t kParameterMerged=kTRUE;
  FairParRootFileIo* output=new FairParRootFileIo(kParameterMerged);
  output->open(out_par);
  rtdb->setOutput(output);
	
  // Set Material file Name
  //-----------------------
  fRun->SetMaterials("media_pnd.geo");

  // Create and add detectors
  //-------------------------
  FairModule *Cave= new PndCave("CAVE");
  Cave->SetGeometryFileName("pndcave.geo");
  fRun->AddModule(Cave);
  if (full_panda) {
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
    Pipe->SetGeometryFileName("beampipe_201112.root");
    fRun->AddModule(Pipe);
    //-------------------------  STT       -----------------
    FairDetector *Stt= new PndStt("STT", kTRUE);
    Stt->SetGeometryFileName("straws_skewed_blocks_35cm_pipe.geo");
    fRun->AddModule(Stt);
    //-------------------------  MVD       -----------------
    FairDetector *Mvd = new PndMvdDetector("MVD", kTRUE);
    Mvd->SetGeometryFileName("Mvd-2.1_FullVersion.root");
    fRun->AddModule(Mvd);
    //-------------------------  GEM       -----------------
    FairDetector *Gem = new PndGemDetector("GEM", kTRUE);
    Gem->SetGeometryFileName("gem_3Stations.root");
    fRun->AddModule(Gem);
  }

  //-------------------------  EMC       -----------------
  PndEmc *Emc = new PndEmc("EMC",kTRUE);
  Emc->SetGeometryVersion(1);
  Emc->SetStorageOfData(kFALSE);
  fRun->AddModule(Emc);

  if (full_panda) {
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
    Rich->SetGeometryFileName("rich_v2.geo");
    fRun->AddModule(Rich);
  }
  
  // Create and Set Event Generator
  //-------------------------------
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  fRun->SetGenerator(primGen);
	
  // Box Generator. first number: PDG particle code: 2nd number: particle multiplicity per event
  FairBoxGenerator* boxGen = new FairBoxGenerator(charge*211, 1); // 13 = muon // 1 = multipl. // 211 = pi+ // -211 = pi-
  boxGen->SetPRange(mom,mom); // GeV/c
  // boxGen->SetPtRange(1.,1.); // GeV/c
  boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
  boxGen->SetThetaRange(85., 95.); // Polar angle in lab system range [degree] - restrict to small rapidity
  boxGen->SetXYZ(0., 0., 0.); // vertex coordinates [mm]
  primGen->AddGenerator(boxGen);  
	
  fRun->SetStoreTraj(kTRUE); // to store particle trajectories 
  fRun->SetBeamMom(15);	

  //---------------------Create and Set the Field(s)---------- 
  PndMultiField *fField= new PndMultiField("FULL");
  fRun->SetField(fField);
	
  //----------- Add Hit producer task to the simulation ------
  PndEmcHitProducer* emcHitProd = new PndEmcHitProducer();
  emcHitProd->SetStorageOfData(kFALSE);
  fRun->AddTask(emcHitProd);

  PndEmcHitsToWaveform* emcHitsToWaveform= new PndEmcHitsToWaveform(0);
  PndEmcWaveformToDigi* emcWaveformToDigi=new PndEmcWaveformToDigi(0);
  //emcHitsToWaveform->SetStorageOfData(kFALSE);
  //emcWaveformToDigi->SetStorageOfData(kFALSE);
  fRun->AddTask(emcHitsToWaveform);  // full digitization
  fRun->AddTask(emcWaveformToDigi);  // full digitization
 
  PndEmcMakeCluster* emcMakeCluster= new PndEmcMakeCluster(0);
  //emcMakeCluster->SetStorageOfData(kFALSE);
  fRun->AddTask(emcMakeCluster);

  PndEmcHdrFiller* emcHdrFiller = new PndEmcHdrFiller();
  fRun->AddTask(emcHdrFiller); // ECM header

  PndEmcMakeBump* emcMakeBump= new PndEmcMakeBump();
  //emcMakeBump->SetStorageOfData(kFALSE);
  fRun->AddTask(emcMakeBump);

  PndEmcMakeRecoHit* emcMakeRecoHit= new PndEmcMakeRecoHit();
  fRun->AddTask(emcMakeRecoHit);
	
  /**Initialize the session*/
  fRun->Init();
  PndEmcMapper *emcMap = PndEmcMapper::Init(1);

  /**After initialization now we can save the field parameters */
  PndMultiFieldPar* Par = (PndMultiFieldPar*) rtdb->getContainer("PndMultiFieldPar");
  if (fField) {  Par->SetParameters(fField); }
  Par->setInputVersion(fRun->GetRunId(),1);
  Par->setChanged();

  /**All parameters are initialized and ready to be saved*/
  rtdb->saveOutput();
  rtdb->print();
		
  // Transport nEvents
  // -----------------
  fRun->Run(nEvents);
	
  timer.Stop();
	
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",timer.RealTime(),timer.CpuTime());
	
}  


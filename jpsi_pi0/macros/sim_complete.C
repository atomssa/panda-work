// Macro for running Panda simulation  with Geant3  or Geant4 (M. Al-Turany)
// This macro is supposed to run the full simulation of the panda detector
// to run the macro:
// root  sim_complete.C  or in root session root>.x  sim_complete_stt.C
// to run with different options:(e.g more events, different momentum, Geant4)
// root  sim_complete.C"(100, "TGeant4",2)"

sim_complete(int type /*1:Signal 0:Background */)
{

  //string file_in,
  Float_t mom = 5.513;

  int nEvents = 10000;

  if (type == 0 ) {
    // Figure out the number of events from input
    // assume the input file name is filt_complete.root
    // This is true for background simulation
    TFile *f = TFile::Open("filt_complete.root");
    TTree *t = (TTree*) f->Get("data");
    nEvents = t->GetEntries();
    f->Close();
  }

  //-----User Settings:-----------------------------------------------
  TString  SimEngine ="TGeant4";
  TString  OutputFile     ="sim_complete.root";
  TString  ParOutputfile  ="simparams.root";
  TString  MediaFile      ="media_pnd.geo";
  gDebug                  = 0;
  TString digiFile        = "all.par"; //The emc run the hit producer directly
  // choose your event generator
  Bool_t UseEvtGen	      =kFALSE;
  Bool_t UseEvtGenDirect      = (type==1); // Signal simulation
  Bool_t UseDpm 	      = (type==0); // Background simulation
  Bool_t UseBoxGenerator      =kFALSE;

  Double_t BeamMomentum = 0.; // beam momentum ONLY for the scaling of the dipole field.
  if (UseBoxGenerator)
    {
      BeamMomentum   =15.0; // ** change HERE if you run Box generator
    }
  else
    {
      BeamMomentum = mom;  // for DPM/EvtGen BeamMomentum is always = mom
    }

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
  //-------------------------  MVD       -----------------
  FairDetector *Mvd = new PndMvdDetector("MVD", kTRUE);
  Mvd->SetGeometryFileName("Mvd-2.1_FullVersion.root");
  fRun->AddModule(Mvd);
  //-------------------------  GEM       -----------------
  FairDetector *Gem = new PndGemDetector("GEM", kTRUE);
  Gem->SetGeometryFileName("gem_3Stations.root");
  fRun->AddModule(Gem);
  //-------------------------  EMC       -----------------
  PndEmc *Emc = new PndEmc("EMC",kTRUE);
  Emc->SetGeometryVersion(1);
  Emc->SetStorageOfData(kFALSE);
  fRun->AddModule(Emc);
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

  if(UseBoxGenerator){	// Box Generator
    FairBoxGenerator* boxGen = new FairBoxGenerator(22, 5); // 13 = muon; 1 = multipl.
    boxGen->SetPRange(mom,mom); // GeV/c
    boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
    boxGen->SetThetaRange(0., 90.); // Polar angle in lab system range [degree]
    boxGen->SetXYZ(0., 0., 0.); // cm
    primGen->AddGenerator(boxGen);
  }
  if(UseDpm){
    PndDpmGenerator* dpmGen = new PndDpmGenerator("filt_complete.root");
    primGen->AddGenerator(dpmGen);
  }
  if(UseEvtGen){
    TString  EvtInput =gSystem->Getenv("VMCWORKDIR");
    EvtInput+="/input/psi2s_jpsi2pi_1k.evt";
    FairEvtGenGenerator* evtGen = new FairEvtGenGenerator(EvtInput.Data());
    primGen->AddGenerator(evtGen);
  }
  if(UseEvtGenDirect){
    //TString  EvtInput =gSystem->Getenv("VMCWORKDIR");
    //EvtInput+="/macro/run/psi2s_Jpsi2pi_Jpsi_mumu.dec";
    TString EvtInput = "/vol0/panda/work/jpsi_pi0/macros/jpsi_pi0.dec";
    PndEvtGenDirect *EvtGen = new PndEvtGenDirect("pbarpSystem", EvtInput.Data(), mom);
    EvtGen->SetStoreTree(kTRUE);
    primGen->AddGenerator(EvtGen);
  }

  //---------------------Create and Set the Field(s)----------
  PndMultiField *fField= new PndMultiField("AUTO");
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

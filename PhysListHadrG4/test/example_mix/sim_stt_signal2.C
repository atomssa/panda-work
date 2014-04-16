// Macro for running Panda simulation  with Geant3  or Geant4 (M. Al-Turany)
// This macro is supposed to run the full simulation of the panda detector with the STT option
// to run the macro:
// root  sim_complete_stt.C  or in root session root>.x  sim_complete_stt.C
// to run with different options:(e.g more events, different momentum, Geant4)
// root  sim_complete_stt.C"(100, "TGeant4",2)"

sim_stt_signal2(Int_t nEvents = 20, TString  SimEngine ="TGeant3", Float_t mom=10.0, TString OutputFile="sim_stt_s2.root" )
{
  //-----User Settings:-----------------------------------------------
  TString  ParOutputfile  ="simparams_stt.root";
  Double_t BeamMomentum   =15.0;
  TString  MediaFile      ="media_pnd.geo";
  gDebug                  = 0;
  TString digiFile        = "emc.par"; //The emc run the hit producer directly 
  // choose your event generator 
  Bool_t UseEvtGen	      =kTRUE;     
  Bool_t UseDpm 	      =kFALSE;
  Bool_t UseBoxGenerator  =kFALSE;
  
  //------------------------------------------------------------------

  TStopwatch timer;
  timer.Start();
 
  // Load basic libraries---------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  
  // Create the Simulation run manager--------------------------------
  FairRunSim *fRun = new FairRunSim();
  fRun->SetName(SimEngine.Data() );
  fRun->SetOutputFile(OutputFile.Data());
  fRun->SetBeamMom(BeamMomentum);
  fRun->SetMaterials(MediaFile.Data());
  FairRuntimeDb *rtdb=fRun->GetRuntimeDb();
  
  // Set the parameters 
  //-------------------------------
  TString emcDigiFile = gSystem->Getenv("VMCWORKDIR");
  emcDigiFile += "/macro/params/";
  emcDigiFile += digiFile;
 
 
  //-------Set the parameter output --------------------
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(emcDigiFile.Data(),"in");
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
  fRun->SetGenerator(primGen);
	 
  if(UseBoxGenerator){	// Box Generator
     FairBoxGenerator* boxGen = new FairBoxGenerator(22, 5); // 13 = muon; 1 = multipl.
     boxGen->SetPtRange(mom,mom); // GeV/c
     boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
     boxGen->SetThetaRange(0., 90.); // Polar angle in lab system range [degree]
     boxGen->SetXYZ(0., 0., 0.); // mm o cm ??
     primGen->AddGenerator(boxGen);
  }
  if(UseDpm){
  	  PndDpmDirect *Dpm= new PndDpmDirect(mom,1);
	  primGen->AddGenerator(Dpm);
  }
  if(UseEvtGen){	
	  TString  EvtInput =gSystem->Getenv("VMCWORKDIR");
	  EvtInput+="/input/psi2s_jpsi2pi_1k.evt";	
	  FairEvtGenGenerator* evtGen = new FairEvtGenGenerator(EvtInput.Data());
	  primGen->AddGenerator(evtGen);
  }	
	

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
 
}  
  

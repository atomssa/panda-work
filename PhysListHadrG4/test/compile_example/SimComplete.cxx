#include "SimComplete.h"

#include <PndCave.h>
#include <PndDchDetector.h>
#include <PndDpmDirect.h>
#include <PndDrc.h>
#include <PndDsk.h>
#include <PndEmc.h>
#include <PndEmcHitProducer.h>
#include <PndGemDetector.h>
#include <PndMagnet.h>
#include <PndMdt.h>
#include <PndMultiField.h>
#include <PndMvdDetector.h>
#include <PndPipe.h>
#include <PndTof.h>
#include <PndStt.h>
#include <PndTpcDetector.h>

#include <FairBoxGenerator.h>
#include <FairEvtGenGenerator.h>
#include <FairModule.h>
#include <FairParAsciiFileIo.h>
#include <FairParRootFileIo.h>
#include <FairPrimaryGenerator.h>
#include <FairRunAna.h>
#include <FairRunSim.h>
#include <FairRuntimeDb.h>

#include <TROOT.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include <iostream>


using namespace std;


void SimComplete(Int_t nEvents, TString const &simEngine, Double_t momentum,
		 Bool_t useEvtGen, Bool_t useDpm, Bool_t useBoxGenerator,
		 Double_t beamMomentum, TString const &outFile, TString const &outParamsFile,
		 TString const &inDigiParamsFile, TString const &trackDetector)
{
  gDebug = 0;

  //------------------------------------------------------------------

  TStopwatch timer;
  timer.Start();
 
  // Create the Simulation run manager--------------------------------
  FairRunSim *fRun = new FairRunSim();
  fRun->SetName(simEngine.Data() );
  fRun->SetOutputFile(outFile.Data());
  fRun->SetBeamMom(beamMomentum);
  fRun->SetMaterials("media_pnd.geo");
  FairRuntimeDb *rtdb=fRun->GetRuntimeDb();

   // Create and add detectors

 //-------------------------  CAVE      -----------------

  FairModule *Cave= new PndCave("CAVE");
  Cave->SetGeometryFileName("pndcave.geo");
  fRun->AddModule(Cave); 
  //-------------------------  Magnet   ----------------- 
/* FairModule *Magnet= new PndMagnet("MAGNET");
  Magnet->SetGeometryFileName("FullSolenoid.root");
  fRun->AddModule(Magnet);
*/
  FairModule *Dipole= new PndMagnet("MAGNET");
  Dipole->SetGeometryFileName("dipole.geo");
  fRun->AddModule(Dipole);

  if (0==trackDetector.CompareTo("stt"))
  {
      //-------------------------  STT       -----------------
      FairDetector *Stt= new PndStt("STT", kTRUE);
      Stt->SetGeometryFileName("straws_skewed_blocks_pipe_120cm.geo");
      fRun->AddModule(Stt);
  }
  else if (0==trackDetector.CompareTo("tpc"))
  {
      //-------------------------  TPC       -----------------
      PndTpcDetector *PndTpc = new PndTpcDetector("TPC", kTRUE);
      PndTpc->SetGeometryFileName("tpc.geo");
      if(simEngine=="TGeant3")  PndTpc->SetAliMC();
      fRun->AddModule(PndTpc);
  }
  //-------------------------  MVD       -----------------
  FairDetector *Mvd = new PndMvdDetector("MVD", kTRUE);
  Mvd->SetGeometryFileName("MVD_v1.0_woPassiveTraps.root");
  fRun->AddModule(Mvd);
 //-------------------------  EMC       -----------------
  PndEmc *Emc = new PndEmc("EMC",kTRUE);
  Emc->SetGeometryFileNameTriple("emc_module125.dat","emc_module3new.root","emc_module4_StraightGeo24.4.root"); //MapperVersion: 6
  Emc->SetStorageOfData(kFALSE);
  fRun->AddModule(Emc);
 //-------------------------  TOF       -----------------  
  FairDetector *Tof = new PndTof("TOF",kTRUE);
  Tof->SetGeometryFileName("tofbarrel.geo");
  fRun->AddModule(Tof);
 //-------------------------  DRC       -----------------
  PndDrc *Drc = new PndDrc("DIRC", kTRUE);
  Drc->SetGeometryFileName("dirc.geo"); 
  Drc->SetRunCherenkov(kFALSE);
  fRun->AddModule(Drc); 
  //-------------------------  MDT       -----------------
  PndMdt *Muo = new PndMdt("MDT",kTRUE);
  Muo->SetBarrel("fast");
  Muo->SetEndcap("fast");
  Muo->SetMuonFilter("fast");
  Muo->SetMdtMagnet(kTRUE);
  Muo->SetMdtMFIron(kTRUE);
  fRun->AddModule(Muo);
   //-------------------------  DCH       -----------------
  FairDetector *Dch = new PndDchDetector("DCH", kTRUE);
  Dch->SetGeometryFileName("dch.root"); 
  fRun->AddModule(Dch);
 
   //-------------------------  GEM      -----------------
  FairDetector *Gem = new PndGemDetector("GEM", kTRUE);
  Gem->SetGeometryFileName("gem_3Stations.root");
  fRun->AddModule(Gem);
 
  //-------------------------  DSK      -----------------
  PndDsk* Dsk = new PndDsk("DSK", kTRUE);
  Dsk->SetGeometryFileName("dsk.geo");
  fRun->AddModule(Dsk);


	
	// Create and Set Event Generator
	//-------------------------------
	FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
	fRun->SetGenerator(primGen);
	
	if(useBoxGenerator){	// Box Generator
		FairBoxGenerator* boxGen = new FairBoxGenerator(22, 5); // 13 = muon; 1 = multipl.
		boxGen->SetPtRange(momentum,momentum); // GeV/c
		boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
		boxGen->SetThetaRange(0., 90.); // Polar angle in lab system range [degree]
		boxGen->SetXYZ(0., 0., 0.); // mm o cm ??
		primGen->AddGenerator(boxGen);
	}
	if(useDpm){
		PndDpmDirect *Dpm= new PndDpmDirect(momentum,1);
		primGen->AddGenerator(Dpm);
	}
	if(useEvtGen){	
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
  
  // Set the parameters 
  //-------------------------------
  TString emcDigiFile = gSystem->Getenv("VMCWORKDIR");
  emcDigiFile += "/macro/params/";
  emcDigiFile += inDigiParamsFile;
 
 
  //-------Set the parameter output --------------------
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(emcDigiFile.Data(),"in");
  rtdb->setFirstInput(parIo1);        

 //---------------------Set Parameter output      ---------- 
  Bool_t kParameterMerged=kTRUE;
  FairParRootFileIo* output=new FairParRootFileIo(kParameterMerged);
  output->open(outParamsFile.Data());
  rtdb->setOutput(output);

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
  
  delete fRun;
}

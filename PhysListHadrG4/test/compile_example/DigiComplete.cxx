#include "DigiComplete.h"

#include <PndEmcHitsToWaveform.h>
#include <PndEmcWaveformToDigi.h>
#include <PndDchDigiProducer.h>
#include <PndDchCylinderHitProducer.h>
#include <PndMvdDigiTask.h>
#include <PndTofHitProducerIdeal.h>
#include <PndMdtHitProducerIdeal.h>
#include <PndGemDigitize.h>
#include <PndGemFindHits.h>
#include <PndSttHitProducerRealFast.h>
#include <PndTpcClusterizerTask.h>
#include <PndTpcDriftTask.h>
#include <PndTpcGemTask.h>
#include <PndTpcPadResponseTask.h>
#include <PndTpcElectronicsTask.h>

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


void DigiComplete(TString const &inFile, TString const &parFile,
		  TString const &digiFile, TString const &outFile)
{
  gDebug = 0;

  //------------------------------------------------------------------

  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 0; // just forget about it, for the moment
      
  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  
  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(inFile);
  fRun->SetOutputFile(outFile);
  
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
  

  // -----   EMC hit producers   ---------------------------------
  // The file name should be the same of the geometry file which was used for the simulation
  
  PndEmcHitsToWaveform* emcHitsToWaveform= new PndEmcHitsToWaveform(iVerbose);
  emcHitsToWaveform->SetStorageOfData(kFALSE);
  PndEmcWaveformToDigi* emcWaveformToDigi=new PndEmcWaveformToDigi(iVerbose);
  fRun->AddTask(emcHitsToWaveform);  // full digitization
  fRun->AddTask(emcWaveformToDigi);  // full digitization

  PndSttHitProducerRealFast* sttHitProducer = new PndSttHitProducerRealFast();
  fRun->AddTask(sttHitProducer);

  QAPlotCollection* qa=new QAPlotCollection("TpcDigiQAPlots");

  PndTpcClusterizerTask* tpcClusterizer = new PndTpcClusterizerTask();
  tpcClusterizer->SetPersistence();

  //ONLY USE THIS WHEN USING ALICE SETTINGS WITH GEANT3
  tpcClusterizer->SetMereChargeConversion();  
  
  fRun->AddTask(tpcClusterizer);

  /**   use Alice Style MC    
                                make one hit per collision with atom
                                use other straggling
                WARNING:        
            1. geant3 has to be used!
            2. LOSS = 5 has to be set!
            3. DCUTE und DCUTM should be set to 10 keV. 
            4. For Digitaization: PndTpcClusterizerTask
                        tpcClusterizer->SetMereChargeConversion() has to be set!
            5. if you do not use this option make sure 2., 4. are not set!
                :-(     
            6. SetMaxNStep should be set to a high value
  */ 

  PndTpcDriftTask* tpcDrifter = new PndTpcDriftTask();
  tpcDrifter->SetPersistence();
  tpcDrifter->SetDistort(false);
  tpcDrifter->SetQAPlotCol(qa);
  fRun->AddTask(tpcDrifter);

  PndTpcPadResponseTask* tpcPadResponse = new PndTpcPadResponseTask();
  tpcPadResponse->SetPersistence();
  fRun->AddTask(tpcPadResponse);

  PndTpcElectronicsTask* tpcElec = new PndTpcElectronicsTask();
  tpcElec->SetPersistence();
  fRun->AddTask(tpcElec);

  PndDchDigiProducer* digiProducer= new PndDchDigiProducer();
  fRun->AddTask(digiProducer);

  PndDchCylinderHitProducer* cylHitProducer= new PndDchCylinderHitProducer();
  fRun->AddTask(cylHitProducer);

  PndMvdDigiTask* mvddigi = new PndMvdDigiTask();
  mvddigi->SetVerbose(iVerbose);
  fRun->AddTask(mvddigi);
  
  PndTofHitProducerIdeal* tofhit = new PndTofHitProducerIdeal();
  tofhit->SetVerbose(iVerbose);
  fRun->AddTask(tofhit);

  PndMdtHitProducerIdeal* mdtHitProd = new PndMdtHitProducerIdeal();
  mdtHitProd->SetPositionSmearing(0.2); // position smearing [cm]
  fRun->AddTask(mdtHitProd);

  PndGemDigitize* gemDigitize = new PndGemDigitize("GEM Digitizer", iVerbose);
  fRun->AddTask(gemDigitize);
        
  PndGemFindHits* gemFindHits = new PndGemFindHits("GEM Hit Finder",  iVerbose);
  fRun->AddTask(gemFindHits);
  
  // -----   Intialise and run   --------------------------------------------
  fRun->Init();

  timer.Start();
  fRun->Run(0,0);

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
}

/** Macro for running fast simulation it will only fill the 
*   PndData without any transport(M. Al-Turany)
*/
run_fast(Int_t nEvents = 1000 )
{
  TString BaseDir =  gSystem->Getenv("VMCWORKDIR");

 //-----User Settings:-----------------------------------------------
  TString  OutputFile     ="sim_fast.root";
  gDebug                  = 0;
  
  // choose your event generator 
  Bool_t UseEvtGen	      = kTRUE;     
  Bool_t UseDpm 	      = kFALSE;
  Bool_t UseBoxGenerator  = kFALSE; 
  
  TString EvtInput = BaseDir + "/input/psi2s_jpsi2pi_1k.evt"; //  Input EvtGen
                     
  Double_t MomDpm  = 7.24; // pbar momentum for DPM generator; matches psi(2S) energy
  
  Double_t MomMin  = 0.5;  // minimum momentum for box generator
  Double_t MomMax  = 2.0;  // maximum   "       "
  
  
   // Load basic libraries---------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon(); 
  // Load the rho and fast sim libraries
  gSystem->Load("libRho");
  gSystem->Load("libfsim");
  
  TStopwatch timer;
  timer.Start();
  gDebug=0;
 
  // Create the Simulation run manager--------------------------------
  FairRunSim *fRun = new FairRunSim();
  fRun->SetOutputFile(OutputFile.Data());
  
 // Create and Set Event Generator
  //-------------------------------
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  fRun->SetGenerator(primGen);
  fRun->SetName("TGeant3");
	 
  if(UseBoxGenerator){	// Box Generator
     FairBoxGenerator* boxGen = new FairBoxGenerator(211, 5); // 211 = pion; 1 = multipl.
     boxGen->SetPtRange(MomMin,MomMax); // GeV/c
     boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
     boxGen->SetThetaRange(0., 90.); // Polar angle in lab system range [degree]
     boxGen->SetXYZ(0., 0., 0.); // mm o cm ??
     primGen->AddGenerator(boxGen);
  }
  if(UseDpm){
  	  PndDpmDirect *Dpm= new PndDpmDirect(MomDpm,0);
	  primGen->AddGenerator(Dpm);
  }
  if(UseEvtGen){
	  FairEvtGenGenerator* evtGen = new FairEvtGenGenerator(EvtInput.Data());
	  primGen->AddGenerator(evtGen);
  }	
  
  // ------------- switch off the transport of particles
  primGen->DoTracking(kFALSE);
	
	
 //---------------------Create and Set the Field(s)---------- 
  PndMultiField *fField= new PndMultiField("FULL");
  fRun->SetField(fField);

  //-------- Setup the Fast Simulation Task  --------------
  //-------------------------------------------------------
  PndFastSim* fastSim = new PndFastSim();
  fastSim->SetVerbosity(0);
  fastSim->AddDetector("CmpDet");
  fastSim->EnablePropagation();
  
  fRun->AddTask(fastSim);


  //-------------------------  Initialize the RUN  -----------------  
  fRun->Init();
  //-------------------------  Run the Simulation  ----------------- 
  fRun->Run(nEvents);
  //------------------------Print some info and exit----------------  
  timer.Stop();
  
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
}  
  

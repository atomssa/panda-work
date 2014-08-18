radlen() {

  string version = "TGeant3";   // OR "TGeant3";
  
  TString outFile= Form("radlen_20k_%s.root",version.c_str());
  cout << "hist_filler outFile= " << outFile << endl;

  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");
  gSystem->Load("libradlendata");
        
  FairLogger::GetLogger()->SetLogToFile(kFALSE);
  FairRunAna* fRun = new FairRunAna();

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();

  //fRun->SetInputFile(inFile);

  //fRun->SetInputFile(Form("sim_complete_20k_upto_%s.root",upto_det.c_str()));
  
  //fRun->SetInputFile("sim_complete_upto_GEM_200k_evt.root");
  
  fRun->SetInputFile(Form("output/sim_complete_%s_0.root",version.c_str()));
  fRun->AddFile(Form("output/sim_complete_%s_1.root",version.c_str()));
  fRun->AddFile(Form("output/sim_complete_%s_2.root",version.c_str()));
  fRun->AddFile(Form("output/sim_complete_%s_3.root",version.c_str()));
  fRun->AddFile(Form("output/sim_complete_%s_4.root",version.c_str()));  

  fRun->SetOutputFile(outFile);
        
  RadLenData *rld = new RadLenData();
  fRun->AddTask(rld);

  fRun->Init(); 
  fRun->Run(0,0);

}

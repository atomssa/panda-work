radlen() {

  string upto_det = "emc_no_mvd";
  
  //TString inFile = Form("sim_complete_50k_evt.root");
  //cout << "hist_filler inFile= " << inFile << endl;
  
  TString outFile= Form("radlen_20k_upto_%s.root",upto_det.c_str());
  cout << "hist_filler outFile= " << outFile << endl;

  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");
  gSystem->Load("libradlendata");
        
  FairLogger::GetLogger()->SetLogToFile(kFALSE);
  FairRunAna* fRun = new FairRunAna();

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  //fRun->SetInputFile(inFile);

  fRun->SetInputFile(Form("sim_complete_20k_upto_%s.root",upto_det.c_str()));
  
  //fRun->SetInputFile("sim_complete_upto_GEM_200k_evt.root");
  
  //fRun->SetInputFile("output/sim_complete_0.root");
  //fRun->AddFile("output/sim_complete_1.root");
  //fRun->AddFile("output/sim_complete_2.root");
  //fRun->AddFile("output/sim_complete_3.root");
  //fRun->AddFile("output/sim_complete_4.root");  

  fRun->SetOutputFile(outFile);
        
  RadLenData *rld = new RadLenData();
  fRun->AddTask(rld);

  fRun->Init(); 
  fRun->Run(0,0);

}

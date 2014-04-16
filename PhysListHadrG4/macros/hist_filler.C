void hist_filler(const char *tag="QGSP_INCLXX_EMV_MOM_0.5_pip")
{
  
  TString inFile = Form("output/data_%s.root",tag);
  //TString outFile= Form("hists_clust_mult/hist_%s.root",tag);
  TString outFile= Form("tmp_hist_%s.root",tag);
  
  cout << "hist_filler inFile= " << inFile << endl;
  cout << "hist_filler outFile= " << outFile << endl;

  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");
  gSystem->Load("libhistfillertask");
	
  FairLogger::GetLogger()->SetLogToFile(kFALSE);
  FairRunAna* fRun = new FairRunAna();

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  fRun->SetInputFile(inFile);	
  fRun->SetOutputFile(outFile);
	
  HistFillerTask *hft = new HistFillerTask();
  fRun->AddTask(hft);
	
  fRun->Init(); 
  fRun->Run(0,0);

}


// ./_pandaroot-jan14/build/macro/PhysListHadrG4
// ./pandaroot-apr13-on-ext-apr13/pandaroot/buildpanda/macro/PhysListHadrG4
// ./pandaroot-jan14-on-ext-dec13/build/macro/PhysListHadrG4

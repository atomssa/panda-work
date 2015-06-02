void ana(int iplab = 0, int itype = 0, int brem = 1, int fid=0, int nevts = 10000) {

  gSystem->Load("libanatda");
  gSystem->ListLibraries();

  //double plab[3] = {5.513, 8., 12.};
  //const char *tt[2] = {"pip_pim","jpsi"};
  //const char *cbrem[2] = {"raw","brem"};
  //TString inPidFile = Form("output/pid/pbar_p_%s_pi0_plab%3.1f_%d.root", tt[itype], plab[iplab], fid);
  //TString inParFile = Form("output/par/pbar_p_%s_pi0_plab%3.1f_%d.root", tt[itype], plab[iplab], fid);
  //TString outFile = Form("test/ana_%s_%s_plab%3.1f_%d.root", tt[itype], cbrem[brem], plab[iplab], fid);
  //cout << "inPidFile= " << inPidFile << endl;
  //cout << "inParFile= " << inParFile << endl;
  //cout << " tt= " << tt[itype] <<  " Outfile= " << outFile << endl;

  const char *tt[2] = {"pi0pipm_dpm","pi0jpsi"};
  const char *cbrem[2] = {"raw","brem"};
  const char *dir = Form("/vol0/panda/work/jpsi_pi0/grid.out/tda/%s/runall.%d.%d", tt[itype], iplab, fid);
  TString inPidFile = Form("%s/pid_complete.root",dir);
  TString inParFile = Form("%s/simparams.root",dir);
  TString outFile = Form("%s/ana_%s.root",dir,cbrem[brem]);
  cout << "inPidFile= " << inPidFile << endl;
  cout << "inParFile= " << inParFile << endl;
  cout << "outFile= " << outFile << endl;

  // *** initialization
  FairLogger::GetLogger()->SetLogToFile(kFALSE);
  FairRunAna* fRun = new FairRunAna();
  fRun->SetInputFile(inPidFile);

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parIO = new FairParRootFileIo();
  parIO->open(inParFile);
  FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
  parIOPid->open((TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par").Data(),"in");
  rtdb->setFirstInput(parIO);
  rtdb->setSecondInput(parIOPid);
  rtdb->setOutput(parIO);

  AnaTda *atda = new AnaTda(iplab, itype, brem);
  atda->set_verbosity(0);
  fRun->AddTask(atda);

  fRun->SetOutputFile(outFile);

  fRun->Init();

  fRun->Run(0,nevts);

  cout << "Done with: " << outFile << endl;

}

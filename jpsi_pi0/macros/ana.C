void ana(int type = 0, int fid=0, int nevts = 10000) {

  int verb = 0;
  bool test_run = true;
  int brem = 1;

  //gSystem->Load("libanatda");
  gSystem->Load("libanatda");

  // *** the files coming from the simulation
  TString inPidFile;
  TString inParFile;

  if (type==0) {
    inPidFile  = Form("output/pid/pbar_p_pip_pim_pi0_%d.root",fid);
    inParFile  = Form("output/par/pbar_p_pip_pim_pi0_%d.root",fid);
  } else {
    inPidFile  = "output/pid/pbar_p_jpsi_pi0_100.root";
    inParFile  = "output/par/pbar_p_jpsi_pi0_100.root";
  }

  // *** initialization
  FairLogger::GetLogger()->SetLogToFile(kFALSE);
  FairRunAna* fRun = new FairRunAna();
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  fRun->SetInputFile(inPidFile);

  // *** setup parameter database
  FairParRootFileIo* parIO = new FairParRootFileIo();
  parIO->open(inParFile);
  FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
  parIOPid->open((TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par").Data(),"in");

  rtdb->setFirstInput(parIO);
  rtdb->setSecondInput(parIOPid);
  rtdb->setOutput(parIO);

  AnaTda *atda = new AnaTda(brem, type  == 0);
  atda->set_verbosity(verb);
  fRun->AddTask(atda);

  //if (type==0)
  fRun->SetOutputFile(Form("%s/ana_%s_%s_%d.root", (test_run?"test":"hists"), (type==0?"bg":"jpsi"), (brem==0?"raw":"brem"), fid) );
  //else
  //fRun->SetOutputFile(Form("%s/ana_%s_%s.root", (test_run?"test":"hists"), (type==0?"bg":"jpsi"), (brem==0?"raw":"brem")) );
  fRun->Init();

  fRun->Run(0,nevts);

}

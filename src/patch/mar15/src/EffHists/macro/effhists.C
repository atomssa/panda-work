void effhists(int sp) {

  int eid_param = 1;

  gSystem->Load("libeffhists");
  gSystem->ListLibraries();

  const TString basedir="/vol0/panda/work/jpsi_pi0/grid.out/jacek/";
  const TString subdir = Form("/runall.%d/",ifile);
  TString inPidFile = Form(basedir+EffHists::s_pid[sp]+subdir+"pid_complete.root");
  TString inParFile = Form(basedir+EffHists::s_pid[sp]+subdir+"simparams.root");
  cout << "inPidFile= " << inPidFile << endl;
  cout << "inParFile= " << inParFile << endl;

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

  EffHists *eh = new EffHists(sp);
  eh->set_verbosity(0);
  eh->set_output_name("effhists.root");
  eh->set_probcut(EffHists::iel, 0.5);
  eh->set_probcut(EffHists::imu, 0.5);
  eh->set_probcut(EffHists::ipi, 0.5);
  eh->set_probcut(EffHists::ik, 0.5);
  eh->set_probcut(EffHists::iprot, 0.5);
  fRun->AddTask(eh);

  //fRun->SetOutputFile(outFile);
  fRun->Init();
  fRun->Run(0,nevts);

}

void effhists(int batch, int ifile) {
  // 0 -> /projet/panda/Ermias/tda/pim_flat
  // 1 -> /projet/panda/Ermias/tda/pip_flat
  // 2 -> /projet/panda/Ermias/tda/elec_flat
  // 3 -> /projet/panda/Ermias/tda/posit_flat

  //int eid_param = 1;
  int sp;
  TString basedir;
  switch(batch) {
  case 0:
    sp = EffHists::ipionm;
    basedir = ifile>=200&&ifile<=299?"/vol0/panda/work/jpsi_pi0/grid.out/tda/pim_flat":"/projet/panda/Ermias/tda/pim_flat";
    break;
  case 1:
    sp = EffHists::ipionp;
    basedir = ifile>=200&&ifile<=299?"/vol0/panda/work/jpsi_pi0/grid.out/tda/pip_flat":"/projet/panda/Ermias/tda/pip_flat";
    break;
  case 2:
    sp = EffHists::ielec;
    basedir = "/projet/panda/Ermias/tda/elec_flat";
    break;
  case 3:
    sp = EffHists::iposit;
    basedir = "/projet/panda/Ermias/tda/posit_flat";
    break;
  default:
    cout << "arguments not recognized" <<endl;
    return;
  }

  gSystem->Load("libeffhists");
  gSystem->ListLibraries();

  //const TString basedir="/vol0/panda/work/jpsi_pi0/grid.out/tda/pip_flat";
  const TString subdir = Form("/runall.%d/",ifile);
  //TString inPidFile = basedir+EffHists::s_spc[sp]+subdir+"pid_complete.root";
  //TString inParFile = basedir+EffHists::s_spc[sp]+subdir+"simparams.root";
  TString inPidFile = basedir+subdir+"pid_complete.root";
  TString inParFile = basedir+subdir+"simparams.root";
  TString outFile = basedir+subdir+"effhists.root";

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
  //eh->set_prob_cut(EffHists::iel, 0.5);
  //eh->set_prob_cut(EffHists::imu, 0.5);
  //eh->set_prob_cut(EffHists::ipi, 0.5);
  //eh->set_prob_cut(EffHists::ik, 0.5);
  //eh->set_prob_cut(EffHists::iprot, 0.5);
  fRun->AddTask(eh);

  fRun->SetOutputFile(outFile);
  fRun->Init();
  fRun->Run(0,0);
  cout << "Done with: " << outFile << endl;
}

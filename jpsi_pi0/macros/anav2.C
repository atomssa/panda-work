void anav2(int iplab = 0, int itype = 0, int brem = 1, int fid=0, int nevts = 0) {

  int eid_param = 1;

  gSystem->Load("libanatda");
  gSystem->ListLibraries();

  double plab[3] = {5.513, 8., 12.};
  const char *tt[2] = {"pip_pim","jpsi"};
  TString inPidFile = Form("output/pid/pbar_p_%s_pi0_plab%3.1f_%d.root", tt[itype], plab[iplab], fid);
  TString inParFile = Form("output/par/pbar_p_%s_pi0_plab%3.1f_%d.root", tt[itype], plab[iplab], fid);
  //TString inPidFile = Form("output/.scrut14/pid/pbar_p_%s_pi0_%d.root", tt[itype], plab[iplab], fid);
  //TString inParFile = Form("output/.scrut14/par/pbar_p_%s_pi0_%d.root", tt[itype], plab[iplab], fid);
  cout << "inPidFile= " << inPidFile << endl;
  cout << "inParFile= " << inParFile << endl;

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

  AnaTdav2 *atda = new AnaTdav2(iplab, itype, brem, eid_param);
  atda->set_verbosity(0);
  atda->set_eff_file_name("eff/epem.root");
  atda->set_eff_hist_name("heffelec", false);
  atda->do_apply_pi0evsoa_cut(false);
  atda->do_apply_pi0m_cut(true);
  atda->do_apply_mtot_cut(false);
  atda->do_apply_dth_dph_cut(false);
  fRun->AddTask(atda);

  const char *cbrem[2] = {"raw","brem"};
  TString outFile = Form("test/anav2_%s_%s_plab%3.1f_%d.root", tt[itype], cbrem[brem], plab[iplab], fid);
  cout << " tt= " << tt[itype] <<  " Outfile= " << outFile << endl;
  fRun->SetOutputFile(outFile);

  fRun->Init();

  fRun->Run(0,nevts);

}

void ana(int plab_id = 0, int itype = 0, int fid=0, int nevts = 10000) {

  int verb = 0;
  bool test_run = true;
  int brem = 1;

  gSystem->Load("libanatda");
  gSystem->ListLibraries();

  double plab[3] = {5.513, 8., 12.};
  const char *tt[2] = {"pipm","jpsi"};
  TString inPidFile = Form("output/pid/pbar_p_%s_pi0_plab%3.1f_%d.root", tt[itype], fid , plab[plab_id]);
  cout << "inPidFile= " << inPidFile << endl;

  // *** initialization
  FairLogger::GetLogger()->SetLogToFile(kFALSE);
  FairRunAna* fRun = new FairRunAna();
  fRun->SetInputFile(inPidFile);

  AnaTda *atda = new AnaTda(brem, itype  == 0);
  //AnaTdav2 *atda = new AnaTdav2();
  atda->set_verbosity(verb);
  fRun->AddTask(atda);

  const char *cbrem[2] = {"raw","brem"};
  TString outFile;
  if (test_run)
    outFile = Form("test/ana_%s_%s_plab%3.1f_%d.root", tt[itype], cbrem[brem], plab[plab_id], fid);
  else
    outFile = Form("hists/ana_%s_%s_plab%3.1f_%d.root", tt[itype], cbrem[brem], plab[plab_id], fid);
  cout << " tt= " << tt[itype] <<  " Outfile= " << outFile << endl;
  fRun->SetOutputFile(outFile);

  fRun->Init();

  fRun->Run(0,nevts);

}

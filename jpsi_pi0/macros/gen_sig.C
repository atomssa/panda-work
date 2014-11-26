void gen_sig(int batch)
{

  double pbar_mom[3] = {5.513, 8.0, 12.0};
  Float_t mom = pbar_mom[batch-1];

  // TODO - there are two places where this is filled (here and in decay description file). WHY?
  Int_t nEvents = 100000;

  gSystem->Load("libEvtGenSA");
  gSystem->ListLibraries();

  string decfile, base = gSystem->Getenv("PNDDIR");
  cout << "base= \"" << base << "\"" << endl;
  if (base == "") {
    cout << "Make sure PNDDIR is set to the root PANDA work directory exiting" << endl;
    return;
  } else {
    decfile = Form("%s/work/jpsi_pi0/macros/jpsi_pi0_%3.1f.dec",base.c_str(),mom);
    std::ifstream inf(decfile.c_str());
    if (!inf.good()) {
      cout << "Decay description file " << decfile << " does not exist. exiting" << endl;
      return;
    }
    inf.close();
  }

  PndEvtGenStandAlone *gen = new PndEvtGenStandAlone("pbarpSystem", decfile, mom);
  gen->Initialize();

  gen->SetVerbose(0);

  for (int i=0; i<nEvents; ++i) {
    gen->Exec();
  }

  gen->Finalize();

}

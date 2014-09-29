void gen_sig()
{

  Float_t mom = 5.513;
  // TODO - there are two places where this is filled. WHY?
  Int_t nEvents = 1000000;

  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");
  gSystem->Load("libEvtGenSA");
  gSystem->ListLibraries();

  PndEvtGenStandAlone *gen = new PndEvtGenStandAlone("pbarpSystem", "/Users/tujuba/panda/work/jpsi_pi0/macros/jpsi_pi0.dec", mom);
  gen->Initialize();

  gen->SetVerbose(0);

  for (int i=0; i<nEvents; ++i) {
    gen->Exec();
  }

  gen->Finalize();

}

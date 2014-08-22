void gen_sig()
{

  Float_t mom = 6.5; // TODO: will hav to be fixed
  Int_t nEvents = 10000;
  
  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");
  gSystem->Load("libEvtGenSA");
  gSystem->ListLibraries();

  PndEvtGenStandAlone *gen = new PndEvtGenStandAlone("pbarpSystem", "/vol0/panda/work/jpsi_pi0/macros/jpsi_pi0.dec", mom);
  gen->Initialize();

  gen->SetVerbose(0);
  
  for (int i=0; i<nEvents; ++i) {
    gen->Exec();
  }

  gen->Finalize();

}  
  

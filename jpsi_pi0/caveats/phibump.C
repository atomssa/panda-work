void phibump(TString path="./") {

  TFile *f = TFile::Open(path+"/digi_complete.root");

  TTree *t = (TTree*)f->Get("cbmsim");

  TClonesArray* phiBumpArray = new TClonesArray("PndEmcBump");
  t->SetBranchAddress("EmcPhiBump", &phiBumpArray);

  ofstream outf;
  outf.open(path+"/phibump.log");
  int nent = t->GetEntries();
  for (int ient=0; ient<1000; ++ient) {
    t->GetEntry(ient);
    cout << "Event " << ient << endl;
    int nBump= phiBumpArray->GetEntriesFast();
    for (int iBump=0; iBump<nBump; ++iBump) {
      PndEmcBump* bump = (PndEmcBump*) phiBumpArray->At(iBump);
      outf << "Event " << ient << " bump " << iBump << "/" << nBump << " E= " << bump->energy() << " emcId= " << bump->GetClusterIndex() << endl;
    }
  }
  outf.close();
}

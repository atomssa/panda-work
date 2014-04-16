{

  string data_file_name = "output/data_QGSP_BERT_EMV_MOM_1.5_pip.root";
  
  int debug = 5;
  
  TFile *data_file = TFile::Open(data_file_name.c_str());
  TTree *t = (TTree*) data_file->Get("cbmsim");

  TClonesArray *geotrack_array = new TClonesArray("TGeoTrack");
  t->SetBranchAddress("GeoTracks",&geotrack_array);

  int nevt = t->GetEntriesFast();
  for (int ievt=0; ievt<nevt; ++ievt) {
    t->GetEntry(ievt);
    if (ievt>100) break;
    if (debug>1) cout << "=============== Event " << ievt << "  ====================" << endl;

    //-------- geo tracks -----------//
    const int ngt = geotrack_array->GetEntriesFast();
    if (debug>1) cout << "ngt= " << ngt << endl;
    for (int igt=0; igt<ngt; ++igt) {
      TGeoTrack *tgt = (TGeoTrack*) geotrack_array->At(igt);
      cout << "tgt = " << tgt << endl;
      TObject *track = (TObject*) tgt->GetParticle();
      track->Dump();
    }
    
  }

}

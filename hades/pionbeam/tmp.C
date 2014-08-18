{

  TH1F *h = new TH1F("h","h",100,-10,10);
  for (int i=0; i<10000; ++i) {
    double rnd = (double) ( (int) (gRandom->Rndm() * 13.0) - 6 );
    h->Fill(rnd);
    cout << rnd << endl;
  }

  h->Draw();
}

{

  double max = 1.166666666666666;
  TH1F *h0 = new TH1F("h0","h0",100,-max*1.1,max*1.1);
  TH1F *h = new TH1F("h","h",100,-10,10);
  TF1 *f = new TF1("pol2","[0]+x*x*[1]",-max,max);
  f->SetParameter(0,1);
  f->SetParameter(1,100);
  for (int i=0; i<100000; ++i) {
    //double rnd = (double) ( (int) (gRandom->Rndm() * 13.0) - 6 );
    //double rnd = (double) ( (int) (f->GetRandom() * 13.0) - 6 );
    double rnd0 = f->GetRandom();
    double rnd = (double) ( (int) (rnd0 * 6) );
    h->Fill(rnd);
    h0->Fill(rnd0);
    cout << rnd << endl;
  }

  TCanvas *c = new TCanvas("c","c");
  c->cd();
  h->Draw();
  TCanvas *c0 = new TCanvas("c0","c0");
  c0->cd();
  h0->Draw("same");
}

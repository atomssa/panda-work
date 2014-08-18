{

  double pitch = 0.5;
  
  TH1F *h_in = new TH1F("h_in","h_in",100,-6,6);
  TH1F *h_out = new TH1F("h_out","h_out",100,-6,6);
  for (int i=0; i<100000; ++i) {
    double rnd = 10* ( gRandom->Rndm() - 0.5);
    h_in->Fill(rnd);
    //double out = ( (double) ( floor( rnd/pitch + 0.5) ) * pitch ) + pitch/2;
    double out = ( (double) ( floor( rnd/pitch ) ) * pitch ) + pitch/2 + ((gRandom->Rndm()-0.5)*pitch);
    h_out->Fill(out);
  }

  TCanvas *tc = new TCanvas();
  tc->Divide(1,2);
  tc->cd(1);
  h_in->Draw();
  tc->cd(2);
  h_out->Draw();  
  
}

{
  gROOT->SetStyle("plain");
  
  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(2,2);

  h_phi_the->SetTitle("Track #theta vs. #phi After Acc. cut;#phi;#theta");
  h_phi_the_b->SetTitle("Track #theta vs. #phi Before Acc. cut;#phi;#theta");

  tc->cd(1);
  h_phi_the_b->Draw("colz");
  
  tc->cd(2);
  h_phi_the->Draw("colz");


  tc->cd(3);
  d_dphi->Draw();
  d_dphi->SetLineWidth(2);
  d_dphi->SetTitle("d#phi = #phi_{track} - #phi_{reco_hit}; d#phi; Y(d#phi)");
  d_dphi->Draw();

  tc->cd(4);
  d_dthe->SetLineWidth(2);
  d_dthe->SetTitle("d#theta = #theta_{track} - #theta_{reco_hit}; d#theta; Y(d#theta)");
  d_dthe->Draw();

  
}

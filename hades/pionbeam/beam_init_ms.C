void beam_init_ms() {

  TCanvas *tcx = new TCanvas("tcx","tcx");
  tcx->cd();
  gPad->SetLogy();
  t->Draw("xms[0]>>hx0(200,-15,15)","","");
  t->Draw("xms[1]>>hx1(200,-15,15)","","same");
  t->Draw("xms[2]>>hx2(200,-15,15)","","same");
  hx0->SetLineColor(1);
  hx1->SetLineColor(2);
  hx2->SetLineColor(4);

  TCanvas *tcy = new TCanvas("tcy","tcy");
  tcy->cd();
  gPad->SetLogy();
  t->Draw("yms[0]>>hy0(200,-25,25)","","");
  t->Draw("yms[1]>>hy1(200,-25,25)","","same");
  t->Draw("yms[2]>>hy2(200,-25,25)","","same");
  hy0->SetLineColor(1);
  hy1->SetLineColor(2);
  hy2->SetLineColor(4);

  TCanvas *tcr = new TCanvas("tcr","tcr");
  tcr->cd();
  gPad->SetLogy();
  t->Draw("sqrt(xms[0]*xms[0]+yms[0]*yms[0])>>hr0(200,0,25)","","");
  t->Draw("sqrt(xms[1]*xms[1]+yms[1]*yms[1])>>hr1(200,0,25)","","same");
  t->Draw("sqrt(xms[2]*xms[2]+yms[2]*yms[2])>>hr2(200,0,25)","","same");
  hr0->Fit("gaus","+");
  hr1->Fit("gaus","+");
  hr2->Fit("gaus","+");
  hr0->SetLineColor(1);
  hr1->SetLineColor(2);
  hr2->SetLineColor(4);
  hr0->Draw();
  hr1->Draw("same");
  hr2->Draw("same");

}

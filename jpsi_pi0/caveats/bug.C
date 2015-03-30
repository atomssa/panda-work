void bug() {

  gStyle->SetOptStat(0);

  TFile *fD = TFile::Open("esim_trunk/default/output_ana.root");
  TFile *fP = TFile::Open("esim_trunk/patch1/output_ana.root");

  TLegend *tlD = new TLegend(0.1,0.89,0.9,1.0); //0.15,0.8,0.8,0);
  tlD->SetHeader("trunk 26841");
  tlD->SetTextSize(0.08);
  TH1F *hDall = (TH1F*)fD->Get("eth_all");
  hDall->SetLineColor(2);
  hDall->SetLineWidth(2);
  hDall->SetTitle(";#theta_{lab}[deg];num. electron");
  tlD->AddEntry(hDall,"e^{-} (ALL)","l");
  TH1F *hDeid = (TH1F*)fD->Get("eth_eid");
  hDeid->SetLineColor(2);
  hDeid->SetLineWidth(2);
  hDeid->SetLineStyle(9);
  tlD->AddEntry(hDeid,"e^{-} (EID)","l");

  TLegend *tlP = new TLegend(0.1,0.89,0.9,1.0);
  tlP->SetTextSize(0.08);
  tlP->SetHeader("trunk 26841 + patch");
  TH1F *hPall = (TH1F*)fP->Get("eth_all");
  hPall->SetLineColor(4);
  hPall->SetLineWidth(2);
  hPall->SetTitle(";#theta_{lab}[deg];num. electron");
  tlP->AddEntry(hPall,"e^{-} (ALL)","l");
  TH1F *hPeid = (TH1F*)fP->Get("eth_eid");
  hPeid->SetLineColor(4);
  hPeid->SetLineWidth(2);
  hPeid->SetLineStyle(9);
  tlP->AddEntry(hPeid,"e^{-} (EID)","l");

  TLegend *tlE = new TLegend(0.1,0.89,0.9,1.0);
  tlE->SetHeader("Efficiency");
  tlE->SetTextSize(0.08);
  TEfficiency *effD = new TEfficiency(*hDeid, *hDall);
  //effD->GetXaxis()->SetRangeUser(12,30);
  effD->SetMarkerStyle(20);
  effD->SetMarkerSize(1);
  effD->SetMarkerColor(2);
  effD->SetLineColor(2);
  effD->SetTitle(";#theta_{lab}[deg];Efficiency");
  tlE->AddEntry(effD, "trunk 26841", "pl");

  TEfficiency *effP = new TEfficiency(*hPeid, *hPall);
  //effP->GetXaxis()->SetRangeUser(12,30);
  effP->SetMarkerStyle(20);
  effP->SetMarkerSize(1);
  effP->SetMarkerColor(4);
  effP->SetLineColor(4);
  tlE->AddEntry(effP, "trunk 26841 + patch", "pl");

  TCanvas *tc = new TCanvas("tc","tc",1000,700);
  tc->Divide(2);

  tc->cd(1);
  hDall->Draw();
  hDeid->Draw("same");
  tlD->Draw();

  tc->cd(2);
  hPall->Draw();
  hPeid->Draw("same");
  tlP->Draw();

  //tc->cd(3);
  //effP->Draw();
  //effD->Draw("same");
  //tlE->Draw();

}

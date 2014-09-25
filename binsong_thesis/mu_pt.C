{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Sep  5 17:42:53 2014) by ROOT version5.34/05
   TCanvas *c1 = new TCanvas("c1", "c1",258,103,958,827);
   c1->Range(-0.2507,-571.725,0.2563,5145.525);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SetTopMargin(0.06);
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);

   TH1F *h_mu_pt_nocor = new TH1F("h_mu_pt_nocor"," ",30,0,3);
   for(Int_t i = 0;i<30;i++) h_mu_pt_nocor->SetBinContent(i,-0.1);

   h_mu_pt_nocor->SetBinContent(6,0.0055);
   h_mu_pt_nocor->SetBinContent(11,0.00744);
   h_mu_pt_nocor->SetBinContent(16,0.00849);
   h_mu_pt_nocor->SetBinContent(21,0.00847);
   h_mu_pt_nocor->SetMinimum(-0.01);
   h_mu_pt_nocor->SetMaximum(0.01);
   h_mu_pt_nocor->SetEntries(4);
   h_mu_pt_nocor->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   h_mu_pt_nocor->SetLineColor(ci);
   h_mu_pt_nocor->SetMarkerStyle(3);
   h_mu_pt_nocor->SetMarkerSize(1.5);
   h_mu_pt_nocor->GetXaxis()->SetNdivisions(505);
   h_mu_pt_nocor->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_mu_pt_nocor->GetXaxis()->SetLabelFont(42);
   h_mu_pt_nocor->GetXaxis()->SetLabelSize(0.05);
   h_mu_pt_nocor->GetXaxis()->SetTitleSize(0.06);
   h_mu_pt_nocor->GetXaxis()->SetTitleFont(42);
   h_mu_pt_nocor->GetYaxis()->SetNdivisions(505);
   h_mu_pt_nocor->GetYaxis()->SetTitle("Peak position");
   h_mu_pt_nocor->GetYaxis()->SetTitleOffset(1.2);
   h_mu_pt_nocor->GetYaxis()->SetLabelFont(42);
   h_mu_pt_nocor->GetYaxis()->SetLabelSize(0.05);
   h_mu_pt_nocor->GetYaxis()->SetTitleSize(0.06);
   h_mu_pt_nocor->GetYaxis()->SetTitleFont(42);
   h_mu_pt_nocor->GetZaxis()->SetLabelFont(42);
   h_mu_pt_nocor->GetZaxis()->SetLabelSize(0.05);
   h_mu_pt_nocor->GetZaxis()->SetTitleSize(0.06);
   h_mu_pt_nocor->GetZaxis()->SetTitleFont(42);
   h_mu_pt_nocor->SetMarkerStyle(20);
   h_mu_pt_nocor->SetMarkerSize(2);
   h_mu_pt_nocor->Draw("P");

   TH1F *h_mu_pt_cor = new TH1F("h_mu_pt_cor"," ",30,0,3);
   for(Int_t j = 0;j<30;j++) h_mu_pt_cor->SetBinContent(j,-0.1);
   h_mu_pt_cor->SetBinContent(6,-0.00051);
   h_mu_pt_cor->SetBinContent(11,-0.000687);
   h_mu_pt_cor->SetBinContent(16,0.00134);
   h_mu_pt_cor->SetBinContent(21,0.00305);
   h_mu_pt_cor->SetEntries(4);

   ci = TColor::GetColor("#000099");
   h_mu_pt_cor->SetLineColor(ci);
   h_mu_pt_cor->SetMarkerColor(2);
   h_mu_pt_cor->SetMarkerStyle(3);
   h_mu_pt_cor->SetMarkerSize(1.5);
   h_mu_pt_cor->GetXaxis()->SetLabelFont(42);
   h_mu_pt_cor->GetXaxis()->SetLabelSize(0.035);
   h_mu_pt_cor->GetXaxis()->SetTitleSize(0.035);
   h_mu_pt_cor->GetXaxis()->SetTitleFont(42);
   h_mu_pt_cor->GetYaxis()->SetLabelFont(42);
   h_mu_pt_cor->GetYaxis()->SetLabelSize(0.035);
   h_mu_pt_cor->GetYaxis()->SetTitleSize(0.035);
   h_mu_pt_cor->GetYaxis()->SetTitleFont(42);
   h_mu_pt_cor->GetZaxis()->SetLabelFont(42);
   h_mu_pt_cor->GetZaxis()->SetLabelSize(0.035);
   h_mu_pt_cor->GetZaxis()->SetTitleSize(0.035);
   h_mu_pt_cor->GetZaxis()->SetTitleFont(42);
   h_mu_pt_cor->SetMarkerStyle(20);
   h_mu_pt_cor->SetMarkerSize(2);
   h_mu_pt_cor->Draw("PSAME");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->Print("mu_pt.png");
}

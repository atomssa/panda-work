{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Sep  5 17:29:30 2014) by ROOT version5.34/05
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

   TH1F *h_pro_pt_cor = new TH1F("h_pro_pt_cor"," ",30,0,3);
   h_pro_pt_cor->SetBinContent(6,0.7357);
   h_pro_pt_cor->SetBinContent(11,0.7787);
   h_pro_pt_cor->SetBinContent(16,0.7971);
   h_pro_pt_cor->SetBinContent(21,0.7578);
   h_pro_pt_cor->SetEntries(4);
   h_pro_pt_cor->SetMaximum(1.0);
   h_pro_pt_cor->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   h_pro_pt_cor->SetLineColor(ci);
   h_pro_pt_cor->SetMarkerColor(2);
   h_pro_pt_cor->SetMarkerStyle(3);
   h_pro_pt_cor->SetMarkerSize(1.7);
   h_pro_pt_cor->GetXaxis()->SetNdivisions(505);
   h_pro_pt_cor->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_pro_pt_cor->GetXaxis()->SetLabelFont(42);
   h_pro_pt_cor->GetXaxis()->SetLabelSize(0.05);
   h_pro_pt_cor->GetXaxis()->SetTitleSize(0.06);
   h_pro_pt_cor->GetXaxis()->SetTitleFont(42);
   h_pro_pt_cor->GetYaxis()->SetNdivisions(505);
   h_pro_pt_cor->GetYaxis()->SetTitle("N_{2#sigma}/N");
   h_pro_pt_cor->GetYaxis()->SetLabelFont(42);
   h_pro_pt_cor->GetYaxis()->SetLabelSize(0.05);
   h_pro_pt_cor->GetYaxis()->SetTitleSize(0.06);
   h_pro_pt_cor->GetYaxis()->SetTitleFont(42);
   h_pro_pt_cor->GetZaxis()->SetLabelFont(42);
   h_pro_pt_cor->GetZaxis()->SetLabelSize(0.05);
   h_pro_pt_cor->GetZaxis()->SetTitleSize(0.06);
   h_pro_pt_cor->GetZaxis()->SetTitleFont(42);
   h_pro_pt_cor->SetMarkerStyle(20);
   h_pro_pt_cor->SetMarkerSize(2);
   h_pro_pt_cor->Draw("P");

   TH1F *h_pro_pt_nocor = new TH1F("h_pro_pt_nocor"," ",30,0,3);
   h_pro_pt_nocor->SetBinContent(6,0.4463);
   h_pro_pt_nocor->SetBinContent(11,0.4752);
   h_pro_pt_nocor->SetBinContent(16,0.4956);
   h_pro_pt_nocor->SetBinContent(21,0.4974);
   h_pro_pt_nocor->SetEntries(4);

   ci = TColor::GetColor("#000099");
   h_pro_pt_nocor->SetLineColor(ci);
   h_pro_pt_nocor->SetLineStyle(2);
   h_pro_pt_nocor->SetMarkerStyle(3);
   h_pro_pt_nocor->SetMarkerSize(1.7);
   h_pro_pt_nocor->GetXaxis()->SetLabelFont(42);
   h_pro_pt_nocor->GetXaxis()->SetLabelSize(0.035);
   h_pro_pt_nocor->GetXaxis()->SetTitleSize(0.035);
   h_pro_pt_nocor->GetXaxis()->SetTitleFont(42);
   h_pro_pt_nocor->GetYaxis()->SetLabelFont(42);
   h_pro_pt_nocor->GetYaxis()->SetLabelSize(0.035);
   h_pro_pt_nocor->GetYaxis()->SetTitleSize(0.035);
   h_pro_pt_nocor->GetYaxis()->SetTitleFont(42);
   h_pro_pt_nocor->GetZaxis()->SetLabelFont(42);
   h_pro_pt_nocor->GetZaxis()->SetLabelSize(0.035);
   h_pro_pt_nocor->GetZaxis()->SetTitleSize(0.035);
   h_pro_pt_nocor->GetZaxis()->SetTitleFont(42);
   h_pro_pt_nocor->SetMarkerStyle(20);
   h_pro_pt_nocor->SetMarkerSize(2);
   h_pro_pt_nocor->Draw("PSAME");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->Print("pro_pt.png");
}

{

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  int iconfig = 8;

  Double_t config[22][3] = {
    {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
    {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
    {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
    {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
    {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
    {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);

  int _type = 0;

  if (_type==0) {
    TCanvas *tc_wtd = new TCanvas("tc_wtd","tc_wtd",1000,1000);
    tc_wtd->cd();
    h_sto->SetTitle(";(p_{MC} - p)/p_{MC}");
    h_cor->Draw();
    tt->DrawLatex(0.57,0.82,Form("Fully corr.",config[iconfig][0]));
    tt->DrawLatex(0.63,0.72,Form("p_{T}= %4.1f",config[iconfig][0]));
    tt->DrawLatex(0.57,0.62,Form("%3.0f < #theta < %3.0f",config[iconfig][1],config[iconfig][2]));
  }

  else if (_type==10) {
    TCanvas *tc_wtd = new TCanvas("tc_wtd","tc_wtd",1000,1000);
    tc_wtd->cd();
    h_cor->SetTitle(";(p_{MC} - p)/p_{MC}");
    h_cor->Draw();
    h_rec->Draw("same");
    h_sto->SetLineColor(4);
    h_sto->Draw("same");
    TLegend* tl_wtd = new TLegend(0.46,0.63,0.95,0.83);
    tl_wtd->SetFillStyle(0);
    tl_wtd->SetBorderSize(0);
    tl_wtd->SetTextSize(0.07);
    tl_wtd->AddEntry(h_rec,"Uncorrected", "pl");
    tl_wtd->AddEntry(h_sto,"Neut. cands", "pl");
    tl_wtd->AddEntry(h_cor,"Bumps", "pl");
    tl_wtd->Draw();
    tt->DrawLatex(0.53,0.52,Form("p_{T}= %4.1f",config[iconfig][0]));
    tt->DrawLatex(0.54,0.42,Form("%3.0f < #theta < %3.0f",config[iconfig][1],config[iconfig][2]));
    gPad->SetGridx();
  }


  else if (_type==1) {
    TCanvas *tc_nmcb = new TCanvas("tc_nmcb","tc_nmcb");
    tc_nmcb->cd();
    h_nmcb_gt1mev->Draw();
    for (int ii = 1; ii <= h_nmcb_gt1mev->GetXaxis()->GetNbins(); ++ii) {
      cout << "Perc. of track (nphot= " << ii-1 << ")= " << 100.0*h_nmcb_gt1mev->GetBinContent(ii)/h_nmcb_gt1mev->GetEntries() << endl;
    }
  }

  else if (_type==2){
    TCanvas *tc_fermi = new TCanvas("tc_fermi","tc_fermi",1000,1000);
    tc_fermi->cd();
    TF1* fun = new TF1("fermi","1.0/(1.+exp((x-[0])/[1]))",0,42);
    fun->SetParameter(0,21);
    fun->SetParameter(1,5);
    fun->SetLineWidth(3);
    fun->SetTitle(";R_{Calc};W(R_{Calc})");
    fun->Draw();
    fun->GetXaxis()->SetTitleSize(0.06);
    fun->GetXaxis()->SetTitleFont(62);
    fun->GetXaxis()->SetLabelSize(0.05);
    fun->GetYaxis()->SetTitleSize(0.06);
    fun->GetYaxis()->SetTitleFont(62);
    fun->GetYaxis()->SetLabelSize(0.05);

    TLatex *ttt = new TLatex();
    ttt->SetTextSize(0.07);
    ttt->SetNDC(kTRUE);
    ttt->DrawLatex(0.36,0.82,Form("W(R) = #frac{1}{1+e^{(R-R0)/R1}}",config[iconfig][0]));
    //ttt->DrawLatex(0.43,0.75,Form("W(R) = #frac{1}{1+e^{#frac{R-p0}{p1}}}",config[iconfig][0]));
    TLatex *tttt = new TLatex();
    tttt->SetTextSize(0.045);
    tttt->SetNDC(kTRUE);
    tttt->DrawLatex(0.470,0.70,"R0 = 21 cm (= R^{ext}_{TRK}/2)");
    tttt->DrawLatex(0.470,0.61,"R1 = 5.0 cm (#approx R^{ext}_{TRK}/10)");
    tc_fermi->Print("fermi.pdf");
  }

  else if (_type==3) {
    TCanvas *tc_rcor = new TCanvas("tc_rcor","tc_rcor",1000,1000);
    tc_rcor->cd();
    h_rad_calc_vs_true->Draw("colz");
    tt->DrawLatex(0.2,0.8,Form("p_{T} = %3.1f GeV/c",config[iconfig][0]));
    tc_rcor->Print(Form("tc_rcor_config%d.pdf",iconfig));
  }

  else if (_type==4) {
    TCanvas *tc_pcor = new TCanvas("tc_pcor","tc_pcor");
    tc_pcor->cd(4);
    //h_cor_vs_rec->Draw("colz");
    prec_vs_rec->Draw("colz");
  }

  else if (_type==5) {
    TCanvas *tc_res = new TCanvas("tc_res","tc_res");
    tc_res->Divide(2,2);
    tc_res->cd(1);
    h_cor->Draw();
    h_sep->Draw("same");
    h_rec->Draw("same");
    gPad->Update();
    TLegend* tl1 = new TLegend(0.6,0.4,0.9,0.6);
    tl1->SetFillStyle(0);
    tl1->SetBorderSize(0);
    tl1->SetTextSize(0.07);
    tl1->AddEntry(h_cor,"Full corr.", "pl");
    tl1->AddEntry(h_sep,"Sep. only", "pl");
    tl1->AddEntry(h_rec,"Reco", "pl");
    tl1->Draw();
    tt->DrawLatex(0.53,0.85,Form("p_{T}= %4.1f",config[iconfig][0]));
    tt->DrawLatex(0.53,0.75,Form("%3.0f < #theta < %3.0f",config[iconfig][1],config[iconfig][2]));

    tc_res->cd(2);
    h_wtd->Draw();
    h_cor->Draw("same");
    gPad->Update();
    TLegend* tl2 = new TLegend(0.6,0.4,0.6,0.6);
    tl2->SetFillStyle(0);
    tl2->SetBorderSize(0);
    tl2->SetTextSize(0.07);
    tl2->AddEntry(h_cor,"Full corr.", "pl");
    tl2->AddEntry(h_wtd,"Full corr. (wtd)", "pl");
    tl2->Draw();
    gPad->SetGridx();

    //tc_res->cd(3);
    //h_cor->Draw();
    //h_mrg->Draw("same");
    //h_sep->Draw("same");
    //TLegend* tl3 = new TLegend(0.6,0.4,0.9,0.6);
    //tl3->SetFillStyle(0);
    //tl3->SetBorderSize(0);
    //tl3->SetTextSize(0.07);
    //tl3->AddEntry(h_cor,"Full corr.", "pl");
    //tl3->AddEntry(h_sep,"Sep. only", "pl");
    //tl3->AddEntry(h_mrg,"Mrg. only", "pl");
    //tl3->Draw();

    tc_res->cd(3);
    h_cor->Draw();
    //h_cor->Fit("gaus","","",-0.05,0.05);
    gPad->SetGridx();

    tc_res->cd(4);
    h_wtd->Draw();
    //h_wtd->Fit("gaus","","",-0.05,0.05);
    gPad->SetGridx();
  }

  else if (_type==6) {
    TCanvas *tc_wtd = new TCanvas("tc_wtd","tc_wtd",1000,1000);
    tc_wtd->cd();
    TString title = h_cor->GetXaxis()->GetTitle();
    h_wtd->GetXaxis()->SetTitle(title);
    h_wtd->Draw();
    h_cor->Draw("same");
    TLegend* tl_wtd = new TLegend(0.46,0.63,0.95,0.83);
    tl_wtd->SetFillStyle(0);
    tl_wtd->SetBorderSize(0);
    tl_wtd->SetTextSize(0.07);
    tl_wtd->AddEntry(h_cor,"Full corr.", "pl");
    tl_wtd->AddEntry(h_wtd,"Full corr. (wtd)", "pl");
    tl_wtd->Draw();
    tt->DrawLatex(0.53,0.52,Form("p_{T}= %4.1f",config[iconfig][0]));
    tt->DrawLatex(0.53,0.42,Form("%3.0f < #theta < %3.0f",config[iconfig][1],config[iconfig][2]));
    gPad->SetGridx();
  }

  else if (_type==7) {
    TCanvas *tc_dphi_dthe = new TCanvas("tc_dphi_dthe","tc_dphi_dthe");
    tc_dphi_dthe->Divide(2,2);
    tc_dphi_dthe->cd(1);
    h_dphi_all->Draw();
    h_dphi_brem->Draw("same");
    tc_dphi_dthe->cd(2);
    h_dthe_all->Draw();
    h_dthe_brem->Draw("same");
    tc_dphi_dthe->cd(3);
    h_dphi_dthe_all->Draw("colz");
    gPad->SetLogz();
    tc_dphi_dthe->cd(4);
    h_dphi_dthe_brem->Draw("colz");
  }

}

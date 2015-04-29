void wtd(int grp) {
  switch(grp) {
  case 0:
    /*
      G1
      |  0 | 0.15 |   30.0 |   45.0 |
      |  1 | 0.25 |   30.0 |   45.0 |
      |  2 |  0.5 |   30.0 |   45.0 |
      |  8 |  1.0 |   30.0 |   45.0 |
      | 12 |  1.5 |   30.0 |   45.0 |
      | 14 |  2.0 |   30.0 |   45.0 |
    */
    const int nc = 6;
    int ic[nc] = {0,1,2,8,12,14};
    int col[nc] = {1,2,5,6,9,7};
    TString fname = "~/grive/physics/PANDA/collab_meeting/2015.03.17/figs//energy_dep_barrel_30to45deg.pdf";
    break;
  case 1:
    /*
      |  2 |  0.5 |   30.0 |   45.0 |
      |  3 |  0.5 |   45.0 |   60.0 |
      |  4 |  0.5 |   60.0 |   75.0 |
      |  5 |  0.5 |   75.0 |   90.0 |
      |  6 |  0.5 |   90.0 |  105.0 |
      |  7 |  0.5 |  105.0 |  120.0 |
    */
    const int nc = 6;
    int ic[nc] = {2,3,4,5,6,7};
    int col[nc] = {1,2,5,6,9,7};
    TString fname = "~/grive/physics/PANDA/collab_meeting/2015.03.17/figs//angle_dep_barrel_500mev.pdf";
    break;
  case 2:
    /*
      |  8 |  1.0 |   30.0 |   45.0 |
      |  9 |  1.0 |   45.0 |   60.0 |
      | 10 |  1.0 |   60.0 |   75.0 |
      | 11 |  1.0 |   75.0 |   90.0 |
    */
    const int nc = 4;
    int ic[nc] = {8,9,10,11};
    int col[nc] = {1,2,5,6};
    TString fname = "~/grive/physics/PANDA/collab_meeting/2015.03.17/figs//angle_dep_barrel_1gev.pdf";
    break;
  case 3:
    /*
      | 15 | 0.15 |   10.0 |   20.0 |
      | 16 | 0.25 |   10.0 |   20.0 |
      | 17 |  0.5 |   10.0 |   20.0 |
      | 18 |  1.0 |   10.0 |   20.0 |
      | 19 |  1.5 |   10.0 |   20.0 |
      | 20 |  2.0 |   10.0 |   20.0 |
      | 21 |  2.5 |   10.0 |   20.0 |
    */
    const int nc = 7;
    int ic[nc] = {15,16,17,18,19,20,21};
    int col[nc] = {1,2,5,6,9,7,42};
    TString fname = "~/grive/physics/PANDA/collab_meeting/2015.03.17/figs//energy_dep_fwd_10to20deg.pdf";
    break;
  default:
    cout << "default" << endl;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  Double_t config[22][3] = {
    {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
    {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
    {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
    {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
    {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
    {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

  TH1F* h_wtd[nc];
  TH1F* h_cor[nc];
  TH1F* h_rec[nc];

  TLegend *tl1 = new TLegend(0.2,0.55,0.5,0.75);
  tl1->SetFillStyle(0);
  tl1->SetBorderSize(0);
  tl1->SetTextSize(0.07);

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);

  TCanvas *tc_wtd = new TCanvas("tc_wtd","tc_wtd");
  //tc_wtd->Divide(3,2);

  TH1F *hdummy[6];
  TPad *pads[6];

  for (int ii = 0; ii < 6; ++ii) {

    TFile *f = new TFile(Form("../grid.out/psim_oct14_binsong_configs/bremcorr.%d_hists.root",ic[ii]));
    h_wtd[ii] = (TH1F*) f->Get("h_wtd");
    h_cor[ii] = (TH1F*) f->Get("h_cor");
    h_rec[ii] = (TH1F*) f->Get("h_rec");

    h_wtd[ii]->GetXaxis()->SetTitle(h_cor[ii]->GetXaxis()->GetTitle());
    //h_wtd[ii]->SetLineWidth(3);
    //h_cor[ii]->SetLineWidth(3);
    //h_rec[ii]->SetLineWidth(3);

    TString xtitle = h_wtd[ii]->GetXaxis()->GetTitle();
    h_wtd[ii]->SetTitle("");
    h_wtd[ii]->GetXaxis()->SetTitle(xtitle);

    double a= 0.15;
    double b= 0.05;
    double _a = a/(1-a);
    double _b = b/(1-b);
    double w = 1/(3+_a+_b);
    double yl = ii<3?0.53:0.0;
    double yh = ii<3?1.0:0.53;
    double xl[3] = {0.0, (_a+1)*w, (_a+2)*w-0.003};
    double xh[3] = {(_a+1)*w, (_a+2)*w, (_a+_b+3)*w};
    cout << "p " << ii << "(" << xl[ii%3] << ", " << yl << ", " << xh[ii%3] << ", " << yh << ")" << " w= " << xh[ii%3] - xl[ii%3] << endl;

    //pads[ii] = new TPad(Form("pads%d",ii),Form("pads%d",ii), 0.33*(ii%3), ii<3?0.0:0.5,  0.3*((ii%3)+1), ii<3?0.5:1.0);
    pads[ii] = new TPad(Form("pads%d",ii),Form("pads%d",ii), xl[ii%3], yl,  xh[ii%3], yh);
    pads[ii]->SetBorderSize(0.1);
    tc_wtd->cd(0);

    pads[ii]->Draw();
    double epsilon=1e-9;
    if ((ii%3)==0) {pads[ii]->SetRightMargin(epsilon); pads[ii]->SetLeftMargin(0.2);}
    if ((ii%3)==1) {pads[ii]->SetLeftMargin(epsilon); pads[ii]->SetRightMargin(epsilon); }
    if ((ii%3)==2) {pads[ii]->SetLeftMargin(epsilon); pads[ii]->SetRightMargin(0.1);}

    if (ii<3) {
      pads[ii]->SetTopMargin(0.05);
      pads[ii]->SetBottomMargin(epsilon);
    } else {
      pads[ii]->SetTopMargin(epsilon);
      pads[ii]->SetBottomMargin(0.15);
    }
    pads[ii]->SetTicks(0,1);
    pads[ii]->cd();

    h_wtd[ii]->SetMaximum(h_wtd[ii]->GetMaximum()*1.3);
    h_wtd[ii]->Draw();
    h_cor[ii]->Draw("same");
    //h_rec[ii]->Draw("same");

    //hdummy[ii] = tmg_yield[ii]->GetHistogram();
    //for (int ibin=0; ibin<5; ++ibin){
    //  int iibin = 0;
    //  if (ibin==0) iibin = 5;
    //  else if (ibin==1) iibin = 32;
    //  else if (ibin==2) iibin = 39;
    //  else iibin = 36+(ibin-2)*30;
    //  hdummy[ii]->GetXaxis()->SetBinLabel(iibin, bin_label[ibin]);
    //}
    //hdummy[ii]->SetLabelSize(0.05,"Y");
    //hdummy[ii]->SetLabelSize(ii==0?0.065:0.08,"X");
    //hdummy[ii]->SetLabelOffset(0.005,"Y");
    //hdummy[ii]->SetTitleSize(0.06,"Y");
    //hdummy[ii]->SetTitleOffset(1.7,"Y");
    //hdummy[ii]->SetTitle(";;counts in 2.96<M_{e^{+}e^{-}}[GeV/c^{2}]<3.22");
    //hdummy[ii]->SetMinimum(0);
    //hdummy[ii]->SetMaximum(max_yield*1.1);

    if (ii==0) {
      //tl1->AddEntry(h_rec[ii],"Uncorrected", "pl");
      tl1->AddEntry(h_wtd[ii],"No weight", "pl");
      tl1->AddEntry(h_cor[ii],"W/ weight", "pl");
      tl1->Draw();
    }
    tt->SetTextColor(1);
    tt->SetTextSize(0.07);
    tt->DrawLatex((ii%3==0)?0.23:0.05,0.87+(ii<3?0.0:0.05),Form("p_{T}= %4.2f GeV/c",config[ic[ii]][0]));
    tt->DrawLatex((ii%3==0)?0.23:0.05,0.77+(ii<3?0.0:0.05),Form("%3.0f < #theta < %3.0f",config[ic[ii]][1],config[ic[ii]][2]));
    tt->SetTextColor(2);
    tt->SetTextSize(0.07);
    if (config[ic[ii]][1]>15) {
      tt->DrawLatex(0.65,0.67+(ii<3?0.0:0.05),"BARREL");
    } else {
      tt->DrawLatex((ii%3==2)?0.7:((ii%3==0)?0.78:0.73),0.67,"FWD.");
      tt->DrawLatex((ii%3==2)?0.6:((ii%3==0)?0.71:0.65),0.60,"ENDCAP");
    }
    gPad->SetGridx();
  }

  tc_wtd->Print(fname);
}

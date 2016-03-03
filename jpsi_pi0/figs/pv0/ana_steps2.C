void set_style(TH1F* h) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetLineWidth(2);
  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  h->GetXaxis()->SetNdivisions(610);
}

void ana_steps2(int ip=0 ){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  // id - pcut - pi0kin - pi0m
  // 18 - 0.9  -   t    -  t
  // 19 - 0    -   t    -  t
  // 20 -
  // 21 - 0.7  -   t    -  t  // TEMP //For BG this was processed wrongly as 0.7-f-t or something like that
  // 22 - 0.9  -   f    -  f  // TEMP //For BG this was processed wrongly as 0.9-t-f
  // 23 - 0.9  -   f    -  t

  // no cut ref: (acceptance only)         hmep_18_0
  // eid cut only at prob 0.9 -            hmep_18_2
  // eid cut only at prob 0.7 -            hmep_21_2
  // eid cut prob 0.9 + one pi0nocut       hmep_22_3
  // eid cut prob 0.9 + one pi0mcut        hmep_23_3
  // eid cut prob 0.9 + one pi0m+kincut    hmep_18_3
  // eid cut prob 0.9 + one pi0m+kin+kfit  hmep_18_8

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";
  //TFile *fsim,*fbg;
  //if (ip==0){
  //  fsim = TFile::Open(Form("%s/hists/note.aug.2015//anav2_jpsi_brem_plab5.5.root",bdir));
  //  fbg = TFile::Open(Form("%s/hists/note.aug.2015//anav2_pip_pim_brem_plab5.5.root",bdir));
  //} else if (ip==1) {
  //  fsim = TFile::Open(Form("%s/hists/note.aug.2015//anav2_jpsi_brem_plab8.0.root",bdir));
  //  fbg = TFile::Open(Form("%s/hists/note.aug.2015//anav2_pip_pim_brem_plab8.0.root",bdir));
  //} else {
  //  fsim = TFile::Open(Form("%s/hists/note.aug.2015//anav2_jpsi_brem_plab12.0.root",bdir));
  //  fbg = TFile::Open(Form("%s/hists/note.aug.2015//anav2_pip_pim_brem_plab12.0.root",bdir));
  //}

  TFile *fsim[5],*fbg[5],*fbg2[5];
  int npass[5] = {18,19,21,22,23};
  for (int ipass=0; ipass < 5; ++ipass) {
    cout << "Openning " << Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_brem4comp_p%d_pass%d.root",bdir,ip,npass[ipass]) << endl;
    fsim[ipass] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_brem4comp_p%d_pass%d.root",bdir,ip,npass[ipass]));
    fbg[ipass] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pipm_brem_p%d_pass%d.root",bdir, ip, npass[ipass]));
    fbg2[ipass] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pi0jpsi_brem4comp_p%d_pass%d.root",bdir, ip, npass[ipass]));
  }

  double nevt[3] = {992780., 990142., 991860.};
  double nevt_sg_xsect[3] = {32780.0, 50142.0, 51860.0};

  double nevt_bg_bis[3] = {814794.0,888292.0,898721.0};
  double nevt_bg[3] = {814794.0,888292.0,898721.0};
  double nevt_bg_xsect[3] = {4.0e11, 1e11, 2e10};

  double nevt_bg2_bis[3] = {200000.0,200000.0,200000.0};
  double nevt_bg2[3] = {9.9e5, 9.9e5, 9.9e5};
  double nevt_bg2_xsect[3] = {4.0e11, 1e11, 2e10};

  double f_valid[3] = {24561.88/32780.0, 33327.48/50142.0, 23742.89/51860.0};
  double nevt_valid[3] = {0.0};
  for (int ii=0; ii < 3; ++ii) nevt_valid[ii] = nevt[ii]*f_valid[ii];

  double sg_scale[3] = {0.0}, bg_scale[3] = {0.0}, bg2_scale[3] = {0.0};
  for (int ii=0; ii < 3; ++ii) {
    sg_scale[ii] = nevt_sg_xsect[ii]/nevt[ii];
    bg_scale[ii] = nevt_bg_bis[ii]/nevt_bg[ii];
    bg2_scale[ii] = nevt_bg2_bis[ii]/nevt_bg2[ii];
  }

  cout << "Getting histograms" << endl;
  TH1F* hmep_sg[5][9];
  TH1F* hmep_bg[5][9];
  TH1F* hmep_bg2[5][9];
  for (int ipass=0; ipass < 5; ++ipass) {
    for (int istep=0; istep<9; ++istep) {

      hmep_sg[ipass][istep] = (TH1F*)fsim[ipass]->Get(Form("hmep_%s%d",(istep>=5?"valid_":""),istep))->Clone(Form("hmep_%s%d_pass%d",(istep>=5?"valid_":""),istep,npass[ipass]));
      set_style(hmep_sg[ipass][istep]);
      hmep_sg[ipass][istep]->Scale(sg_scale[ip]);
      hmep_sg[ipass][istep]->SetTitle("");

      hmep_bg[ipass][istep] = (TH1F*)fbg[ipass]->Get(Form("hmep_%s%d",(istep>=5?"valid_":""),istep))->Clone(Form("hmep_%s%d_pass%d",(istep>=5?"valid_":""),istep,npass[ipass]));
      set_style(hmep_bg[ipass][istep]);
      hmep_bg[ipass][istep]->Scale(bg_scale[ip]);
      hmep_bg[ipass][istep]->SetTitle("");

      hmep_bg2[ipass][istep] = (TH1F*)fbg2[ipass]->Get(Form("hmep_%s%d",(istep>=5?"valid_":""),istep))->Clone(Form("hmep_%s%d_pass%d",(istep>=5?"valid_":""),istep,npass[ipass]));
      set_style(hmep_bg2[ipass][istep]);
      hmep_bg2[ipass][istep]->Scale(bg2_scale[ip]);
      hmep_bg2[ipass][istep]->SetTitle("");

    }
  }
  //cout << "Done getting histograms" << endl;

  //TCanvas *tc_steps = new TCanvas("tc_steps","tc_steps",800,1600);
  //tc_steps->Divide(1,2);
  TCanvas *tc_steps = new TCanvas("tc_steps","tc_steps",1600,700);
  tc_steps->Divide(2,1);
  tc_steps->GetPad(1)->SetRightMargin(0.05);
  tc_steps->GetPad(2)->SetRightMargin(0.05);

  int immin = hmep_sg[0][0]->GetXaxis()->FindBin(2.8);
  int immax = hmep_sg[0][0]->GetXaxis()->FindBin(3.3);

  //  int npass[5] = {18,19,21,22,23};
  // no cut ref: (acceptance only)                  hmep_18(0)_0
  // eid cut only at prob 0.9 -                     hmep_18(0)_2
  // eid cut only at prob 0.7 -                     hmep_21(2)_2
  // eid cut prob 0.9 + one pi0nocut                hmep_22(3)_3
  // eid cut prob 0.9 + one pi0mcut                 hmep_23(4)_3
  // eid cut prob 0.9 + one pi0m+kincut             hmep_18(0)_3
  // eid cut prob 0.9 + one pi0m+kincut+most-btb    hmep_18(0)_4

  // eid cut prob 0.9 + one pi0m+kincut+most-btb     hmep_18(0)_5
  // eid cut prob 0.9 + one pi0m+kin+kfit0+mostbtb   hmep_18(0)_6
  // eid cut prob 0.9 + one pi0m+kin+kfit1+mostbtb   hmep_18(0)_7
  // eid cut prob 0.9 + one pi0m+kin+kfit2+mostbtb   hmep_18(0)_8

  // TODO: why are _3 and _4 different? They shouldn't be

  //TLegend *tl = new TLegend(0.38,0.5,0.85,0.9);
  //tl->SetBorderSize(0);
  //tl->SetFillStyle(0);
  //tl->SetTextSize(0.035);

  // int col = 1;
  // tc_steps->cd(1);
  // hmep_sg[0][0]->SetTitle("#pi^{0}J/#psi, p_{#bar{p}} = 5.513 GeV/c;M_{e^{+}e^{-}}[GeV/c];counts");
  // hmep_sg[0][0]->GetXaxis()->SetRangeUser(2.6,4.6);
  // hmep_sg[0][0]->SetLineColor(col);
  // hmep_sg[0][0]->Draw();
  // tc_steps->cd(2);
  // hmep_bg[0][0]->SetTitle("#pi^{0}#pi^{+}#pi^{-}, p_{#bar{p}} = 5.513 GeV/c;M_{#pi^{+}#pi^{-}}[GeV/c];counts");
  // hmep_bg[0][0]->Scale(4e-7);
  // hmep_bg[0][0]->SetLineColor(col);
  // hmep_bg[0][0]->Draw();
  // gPad->SetLogy();
  // tl->AddEntry(hmep_sg[0][0],"In acc.","l");
  //
  // TLatex *tlat = new TLatex();
  // tlat->DrawLatexNDC(0.68,0.79,"#times 4#times10^{-7}");
  // TArrow *ta = new TArrow(3.6, 650, 3.1, 300, 0.015, ">");
  // ta->Draw();
  //
  // cout << "Nnorm = " << hmep_sg[0][0]->Integral(immin,immax) << endl;
  // cout << "Acc = " << Form("%4.1f",100*hmep_sg[0][0]->Integral(immin,immax)/nevt[ip]) << endl;
  //
  // tc_steps->cd(1);
  // tl->Draw();
  // tc_steps->Update();
  // tc_steps->Print("ana_steps_0.pdf");
  //
  // col = 3;
  // tc_steps->cd(1);
  // hmep_sg[2][2]->SetLineColor(col);
  // hmep_sg[2][2]->Draw("same");
  // tc_steps->cd(2);
  // hmep_bg[2][2]->SetLineColor(col);
  // hmep_bg[2][2]->Draw("same");
  // tl->AddEntry(hmep_sg[2][2],"Pairs after EID (0.7)","l");
  // cout << "EID(0.7) = " << Form("%4.1f",100*hmep_sg[2][2]->Integral(immin,immax)/nevt[ip]) << endl;
  //
  // tc_steps->cd(1);
  // tl->Draw();
  // tc_steps->Update();
  // tc_steps->Print("ana_steps_1.pdf");
  //
  // col = 9;
  // tc_steps->cd(1);
  // hmep_sg[0][2]->SetLineColor(col);
  // hmep_sg[0][2]->Draw("same");
  // tc_steps->cd(2);
  // hmep_bg[0][2]->SetLineColor(col);
  // hmep_bg[0][2]->Draw("same");
  // tl->AddEntry(hmep_sg[0][2],"Pairs after EID (0.9)","l");
  // cout << "EID(0.9) = " << Form("%4.1f",100*hmep_sg[0][2]->Integral(immin,immax)/nevt[ip]) << endl;
  //
  // tc_steps->cd(1);
  // tl->Draw();
  // tc_steps->Update();
  // tc_steps->Print("ana_steps_2.pdf");
  //
  // col = 41;
  // tc_steps->cd(1);
  // hmep_sg[3][4]->SetLineColor(col);
  // hmep_sg[3][4]->Draw("same");
  // tc_steps->cd(2);
  // hmep_bg[3][3]->SetLineColor(col);
  // hmep_bg[3][3]->Draw("same");
  // tl->AddEntry(hmep_sg[3][4],"EID(0.9)+ N#pi^{0}_{no id}>0","l");
  // cout << "EID(0.9)+npi0>0 = " << Form("%4.1f",100*hmep_sg[3][4]->Integral(immin,immax)/nevt[ip]) << endl;
  //
  // tc_steps->cd(1);
  // tl->Draw();
  // tc_steps->Update();
  // tc_steps->Print("ana_steps_3.pdf");
  //
  // col = 6;
  // tc_steps->cd(1);
  // hmep_sg[4][4]->SetLineColor(col);
  // hmep_sg[4][4]->Draw("same");
  // tc_steps->cd(2);
  // hmep_bg[4][3]->SetLineColor(col);
  // hmep_bg[4][3]->Draw("same");
  // tl->AddEntry(hmep_sg[4][4],"EID(0.9)+ N#pi^{0}_{A}>0","l");
  // cout << "EID(0.9)+npi0_mcut>0 = " << Form("%4.1f",100*hmep_sg[4][4]->Integral(immin,immax)/nevt[ip]) << endl;
  //
  // tc_steps->cd(1);
  // tl->Draw();
  // tc_steps->Update();
  // tc_steps->Print("ana_steps_4.pdf");
  //
  // col = 38;
  // tc_steps->cd(1);
  // hmep_sg[0][4]->SetLineColor(col);
  // hmep_sg[0][4]->Draw("same");
  // tc_steps->cd(2);
  // hmep_bg[0][3]->SetLineColor(col);
  // hmep_bg[0][3]->Draw("same");
  // tl->AddEntry(hmep_bg[0][3],"EID(0.9)+ N#pi^{0}_{B}>0","l");
  // cout << "EID(0.9)+npi0_mcut_kincut>0 = " << Form("%4.1f",100*hmep_sg[0][4]->Integral(immin,immax)/nevt[ip]) << endl;
  //
  // tc_steps->cd(1);
  // tl->Draw();
  // tc_steps->Update();
  // tc_steps->Print("ana_steps_5.pdf");
  //
  // col = 46;
  // tc_steps->cd(1);
  // hmep_sg[0][4]->SetLineColor(col);
  // hmep_sg[0][4]->Draw("same");
  // tc_steps->cd(2);
  // hmep_bg[0][4]->SetLineColor(col);
  // hmep_bg[0][4]->Draw("same");
  // tl->AddEntry(hmep_sg[0][4],"EID(0.9)+N#pi^{0}_{B}>0+BtB","l");
  // cout << "EID(0.9)+npi0_mcut_kincut>0 = " << Form("%4.1f",100*hmep_sg[0][4]->Integral(immin,immax)/nevt[ip]) << endl;
  //
  // tc_steps->cd(1);
  // tl->Draw();
  // tc_steps->Update();
  // tc_steps->Print("ana_steps_6.pdf");

  TLegend *tl = new TLegend(0.35,0.55,0.85,0.85);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  tl->SetTextSize(0.038);

  col = 1;
  tc_steps->cd(1);
  hmep_sg[0][4]->SetTitle("#pi^{0}J/#psi, p_{#bar{p}} = 5.513 GeV/c;M_{e^{+}e^{-}}[GeV/c];counts");
  hmep_sg[0][4]->GetXaxis()->SetRangeUser(2.6,4.6);
  hmep_sg[0][4]->SetLineColor(col);
  hmep_sg[0][4]->Draw("same");
  tc_steps->cd(2);
  hmep_bg2[0][4]->SetTitle("#pi^{0}#pi^{0}J/#psi, p_{#bar{p}} = 5.513 GeV/c;M_{#pi^{+}#pi^{-}}[GeV/c];counts");
  hmep_bg2[0][4]->GetXaxis()->SetRangeUser(1.5,4.6);
  hmep_bg2[0][4]->Scale(8e-2);
  hmep_bg2[0][4]->SetLineColor(col);
  hmep_bg2[0][4]->Draw("same");
  gPad->SetLogy();
  tl->SetHeader(" EID(0.9) + N#pi^{0}_{B}>0 and");
  tl->AddEntry(hmep_sg[0][4],"Most back-to-back","l");
  cout << "EID(0.9)+npi0_mcut_kincut>0 = " << Form("%4.1f",100*hmep_sg[0][4]->Integral(immin,immax)/nevt_valid[ip]) << endl;

  TLatex *tlat2 = new TLatex();
  tlat2->DrawLatexNDC(0.67,0.785,"#times 8#times10^{-2}");
  TArrow *ta2 = new TArrow(3.52, 63, 3.175, 40, 0.015, ">");
  ta2->Draw();

  tc_steps->cd(1);
  tl->Draw();
  tc_steps->Update();
  tc_steps->Print("kincut_steps_0.pdf");

  col = 41;
  tc_steps->cd(1);
  hmep_sg[0][6]->SetLineColor(col);
  hmep_sg[0][6]->Draw("same");
  tc_steps->cd(2);
  hmep_bg2[0][6]->SetLineColor(col);
  hmep_bg2[0][6]->Draw("same");
  tl->AddEntry(hmep_sg[0][6],"#chi_{2#gamma}<#chi_{2#gamma,MAX}","l");
  cout << "EID(0.9)+npi0_mcut_kincut>0+chi2_sig = " << Form("%4.1f",100*hmep_sg[0][6]->Integral(immin,immax)/nevt_valid[ip]) << endl;

  tc_steps->cd(1);
  tl->Draw();
  tc_steps->Update();
  tc_steps->Print("kincut_steps_1.pdf");

  col = 38;
  tc_steps->cd(1);
  hmep_sg[0][7]->SetLineColor(col);
  hmep_sg[0][7]->Draw("same");
  tc_steps->cd(2);
  hmep_bg2[0][7]->SetLineColor(col);
  hmep_bg2[0][7]->Draw("same");
  tl->AddEntry(hmep_sg[0][7],"#chi_{2#gamma}<#chi_{2#gamma,MAX} + #chi_{4#gamma}<#chi_{2#gamma}","l");
  cout << "EID(0.9)+npi0_mcut_kincut>0+chi2_sig+chi2_bg = " << Form("%4.1f",100*hmep_sg[0][7]->Integral(immin,immax)/nevt_valid[ip]) << endl;

  tc_steps->cd(1);
  tl->Draw();
  tc_steps->Update();
  tc_steps->Print("kincut_steps_2.pdf");

  col = 46;
  tc_steps->cd(1);
  hmep_sg[0][8]->SetLineColor(col);
  hmep_sg[0][8]->Draw("same");
  tc_steps->cd(2);
  hmep_bg2[0][8]->SetLineColor(col);
  hmep_bg2[0][8]->Draw("same");
  tl->AddEntry(hmep_sg[0][8],"#chi_{2#gamma}<#chi_{2#gamma,MAX} + #chi_{4#gamma}<#chi_{2#gamma} + N#gamma_{>20MeV}<4","l");
  cout << "int = " << hmep_sg[0][8]->Integral(immin,immax) << endl;
  cout << " nevt_valid = " << nevt_valid[ip] << endl;
  cout << "EID(0.9)+npi0_mcut_kincut>0+chi2_sig+chi2_bg+ng20<4 = " << Form("%4.1f",100*hmep_sg[0][8]->Integral(immin,immax)/nevt_valid[ip]) << endl;

  tc_steps->cd(1);
  tl->Draw();
  tc_steps->Print("kincut_steps_3.pdf");

  //hmep[ii][1]->SetLineColor(2);
  //hmep[ii][1]->Draw("same");
  ////if (ii==0) hmep[ii][1]->Scale(0.9);
  //if (ii==0) tl->AddEntry(hmep[ii][1],"Pairs after eid (wt for bg)","l");
  ////if (ii==0) cout << "Eff1. Eid = " << 100*hmep[0][1]->Integral(immin,immax)/hmep[0][0]->Integral(immin,immax) << endl;
  ////if (ii==0) cout << "Eff2. (All) = " << 100*hmep[0][2]->Integral(immin,immax)/hmep[0][0]->Integral(immin,immax) <<  " should be equal to Eff1 " << endl;
  //
  //
  //hmep[ii][3]->SetLineColor(4);
  //hmep[ii][3]->SetLineWidth(3);
  //hmep[ii][3]->Draw("same");
  //if (ii==0) tl->AddEntry(hmep[ii][3],"N_{#pi^{0}}>0 after #pi^{0} cuts","l");
  ////if (ii==0) cout << "Eff3. (Npi0>0) = " << 100*hmep[0][3]->Integral(immin,immax)/hmep[0][0]->Integral(immin,immax)<< endl;
  //if (ii==0) cout << "Eff3. (Npi0>0) = " << 100*hmep[0][3]->Integral(immin,immax)/nevt[ip]<< endl;
  //
  //hmep[ii][4]->SetLineColor(3);
  //hmep[ii][4]->Draw("same");
  //if (ii==0) tl->AddEntry(hmep[ii][4],"Most b-to-b #pi^{0}-(e^{+}e^{-}) pair","l");
  //if (ii==1) tl->Draw();

  //if (ii==0) cout << "Eff4. (MostBtoB) = "     << 100*hmep[0][4]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << " should be equal to Eff3" << endl;
  //if (ii==0) cout << "Eff5. (MostBtoB.bis) = " << 100*hmep[0][5]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << " should be equal to Eff4" << endl;
  //if (ii==0) cout << "Eff6. (chi2sig) = "      << 100*hmep[0][6]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << endl;
  //if (ii==0) cout << "Eff7. (chi2bg) = "       << 100*hmep[0][7]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << endl;
  //if (ii==0) cout << "Eff8. (ngam>0) = "       << 100*hmep[0][8]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << endl;
  //if (ii==0) cout << "Eff4. (MostBtoB) = "     << 100*hmep[0][4]->Integral(immin,immax)/nevt[ip] << " should be equal to Eff3" << endl;
  //if (ii==0) cout << "Eff5. (MostBtoB.bis) = " << 100*hmep[0][5]->Integral(immin,immax)/nevt[ip] << " should be equal to Eff4" << endl;
  //if (ii==0) cout << "Eff6. (chi2sig) = "      << 100*hmep[0][6]->Integral(immin,immax)/nevt[ip] << endl;
  //if (ii==0) cout << "Eff7. (chi2bg) = "       << 100*hmep[0][7]->Integral(immin,immax)/nevt[ip] << endl;
  //if (ii==0) cout << "Eff8. (ngam>0) = "       << 100*hmep[0][8]->Integral(immin,immax)/nevt[ip] << endl;

  //tc_steps->Print(Form("%s/figs/2015.06.08/ana_steps.pdf",bdir));

}

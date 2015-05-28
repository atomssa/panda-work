void dphimax(int cfg=-1) {

  static const int nconfig = 25;
  Double_t config[nconfig][3] = {
    /*0*/{0.15, 30.0, 45.0}, /*1*/{0.25, 30.0, 45.0}, /*2*/{0.5, 30.0, 45.0}, /*3*/{0.5, 45.0, 60.0}, /*4*/{0.5, 60.0, 75.0},
    /*5*/{0.5, 75.0, 90.0}, /*6*/{0.5, 90.0, 105.0}, /*7*/{0.5, 105.0, 120.0}, /*8*/{1.0, 30.0, 45.0}, /*9*/{1.0, 45.0, 60.0},
    /*10*/{1.0, 60.0, 75.0}, /*11*/{1.0, 75.0, 90.0}, /*12*/{1.5, 30.0, 45.0}, /*13*/{1.5, 45.0, 60.0}, /*14*/{2.0, 30.0, 45.0},
    /*15*/{0.15, 10.0, 20.0}, /*16*/{0.25, 10.0, 20.0}, /*17*/{0.5, 10.0, 20.0}, /*18*/{1.0, 10.0, 20.0}, /*19*/{1.5, 10.0, 20.0},
    /*20*/{2.0, 10.0, 20.0}, /*21*/{2.5, 10.0, 20.0}, /*22*/{0.2, 30.0, 45.0}, /*23*/{0.2, 90.0, 105.0}, /*24*/{0.2, 10.0, 20.0}};

  int icfg[25] = {
    /*30-45 (6)*/ 22, 1, 2, 8, 12, 14,
    /*10-20 (7) */ 15, 24, 16, 17, 18, 19, 20,
    /*0.5 GeV/c (7)*/ 17, 2, 3, 4, 5, 6, 7,
    /*1GeV/c (5) */ 18, 8, 9, 10, 11
  };

  TCanvas *tc = new TCanvas("tc","tc");
  if (cfg==-1) tc->Divide(5,5);

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);

  bool sep = true;
  for (int iconfig = 0; iconfig < nconfig; ++iconfig) {
    if (cfg!=-1&&iconfig!=cfg) continue;
    //TString fname = Form("../grid.out/esim_oct14_binsong_configs/all.ibs.xr-0.1_0.1/bremcorr.all.ibs.cfg.%d_hists.root",icfg[iconfig]);
    TString fname = Form("../grid.out/esim_oct14_binsong_configs/xr-0.1_0.1//bremcorr.ibs.0.cfg.%d_hists.root",icfg[iconfig]);
    //TString fname = Form("../grid.out/old.brem/psim_oct14_binsong_configs/bremcorr_merged.%d_hists.root",iconfig);
    cout << "Opening file " << fname << endl;
    TFile *f = TFile::Open(fname);

    TH1F* phi_max_rec = (TH1F*) f->Get("h_dphi_max_rec"); phi_max_rec->SetLineWidth(2);
    TH1F* phi_max_mc = (TH1F*) f->Get("h_dphi_max_mc"); phi_max_mc->SetLineWidth(2);
    TH1F* dphi_all_sep = (TH1F*) f->Get("h_dphi_all_sep"); dphi_all_sep->SetLineWidth(2);
    TH1F* dphi_brem_sep = (TH1F*) f->Get("h_dphi_brem_sep"); dphi_brem_sep->SetLineWidth(2);
    TH1F* dphi_cut_sep = (TH1F*) f->Get("h_dphi_cut_sep"); dphi_cut_sep->SetLineWidth(2);
    TH1F* dphi_all_mrg = (TH1F*) f->Get("h_dphi_all_mrg"); dphi_all_mrg->SetLineWidth(2);
    TH1F* dphi_cut_mrg = (TH1F*) f->Get("h_dphi_cut_mrg"); dphi_cut_mrg->SetLineWidth(2);

    if (cfg==-1)
      tc->cd(1+iconfig);
    else
      tc->cd();

    phi_max_mc->Draw();
    phi_max_rec->Draw("same");
    //phi_max_rec->Draw();
    //phi_max_mc->Draw("same");
    phi_max_rec->SetLineColor(2);
    phi_max_mc->SetMinimum(1);

    if (sep) {
      //dphi_all_sep->Scale(double(phi_max_mc->GetEntries())/double(dphi_all_sep->GetEntries()));
      dphi_all_sep->SetLineColor(4);
      dphi_all_sep->Draw("same");
      //dphi_all_sep->Draw();
      //dphi_brem_sep->Scale(double(phi_max_mc->GetEntries())/double(dphi_all_sep->GetEntries()));
      dphi_brem_sep->SetLineColor(1);
      dphi_brem_sep->Draw("same");
      //dphi_cut_sep->Scale(double(phi_max_mc->GetEntries())/double(dphi_all_sep->GetEntries()));
      dphi_cut_sep->SetLineColor(3);
      dphi_cut_sep->Draw("same");
    } else {
      //dphi_all_mrg->Scale(double(phi_max_mc->GetEntries())/double(dphi_all_mrg->GetEntries()));
      dphi_all_mrg->SetLineColor(4);
      dphi_all_mrg->Draw("same");
      //dphi_all_mrg->Draw();
      //dphi_cut_mrg->Scale(double(phi_max_mc->GetEntries())/double(dphi_all_mrg->GetEntries()));
      dphi_cut_mrg->Scale(0.8);
      dphi_cut_mrg->SetLineColor(1);
      dphi_cut_mrg->Draw("same");
    }

    gPad->SetLogy();
    tt->SetTextColor(1);
    tt->SetTextSize(0.07);
    tt->DrawLatex(0.5,0.82,Form("p_{T}= %4.2f GeV/c",config[icfg[iconfig]][0]));
    tt->DrawLatex(0.5,0.72,Form("%3.0f < #theta < %3.0f",config[icfg[iconfig]][1],config[icfg[iconfig]][2]));
  }

}

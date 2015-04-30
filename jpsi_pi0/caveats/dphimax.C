void dphimax(int cfg=-1) {

  static const int nconfig = 22;
  Double_t config[nconfig][3] = {
    {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
    {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
    {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
    {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
    {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
    {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

  TCanvas *tc = new TCanvas("tc","tc");
  if (cfg==-1) tc->Divide(6,4);

  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);

  bool sep = true;
  for (int iconfig = 0; iconfig < nconfig; ++iconfig) {
    if (cfg!=-1&&iconfig!=cfg) continue;
    TString fname = Form("../grid.out/psim_oct14_binsong_configs/bremcorr_merged.%d_hists.root",iconfig);
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
    phi_max_rec->SetLineColor(2);
    phi_max_mc->SetMinimum(1);

    if (sep) {
      dphi_all_sep->Scale(double(phi_max_mc->GetEntries())/double(dphi_all_sep->GetEntries()));
      dphi_all_sep->SetLineColor(4);
      dphi_all_sep->Draw("same");
      dphi_brem_sep->Scale(double(phi_max_mc->GetEntries())/double(dphi_brem_sep->GetEntries()));
      dphi_brem_sep->SetLineColor(1);
      dphi_brem_sep->Draw("same");
      dphi_cut_sep->Scale(double(phi_max_mc->GetEntries())/double(dphi_cut_sep->GetEntries()));
      dphi_cut_sep->SetLineColor(3);
      dphi_cut_sep->Draw("same");
    } else {
      dphi_all_mrg->Scale(double(phi_max_mc->GetEntries())/double(dphi_all_mrg->GetEntries()));
      dphi_all_mrg->SetLineColor(4);
      dphi_all_mrg->Draw("same");
      dphi_cut_mrg->Scale(double(phi_max_mc->GetEntries())/double(dphi_cut_mrg->GetEntries()));
      dphi_cut_mrg->SetLineColor(1);
      dphi_cut_mrg->Draw("same");
    }

    gPad->SetLogy();
    tt->SetTextColor(1);
    tt->SetTextSize(0.07);
    tt->DrawLatex(0.5,0.82,Form("p_{T}= %3.1f GeV/c",config[iconfig][0]));
    tt->DrawLatex(0.5,0.72,Form("%3.0f < #theta < %3.0f",config[iconfig][1],config[iconfig][2]));
  }

}

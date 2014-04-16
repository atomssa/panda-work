void plot_g3_3vsD() {

  gStyle->SetOptStat(0);
  
  vector<string> tags;
  tags.push_back("G3_HADR3_NUCRIN");
  tags.push_back("QGSP_BIC_EMV");
  tags.push_back("QGSP_BERT_EMV");
  tags.push_back("QGSP_BERT_OPTICAL");
  tags.push_back("QGSP_BERT_EMV_OPTICAL");
  tags.push_back("QGSC_BERT_EMV");
  tags.push_back("QGSP_BERT_HP");
  tags.push_back("G3_HADR6_GCALOR");
  tags.push_back("G3_HADR5_MICAP");

  
  static const int nmom = 6;
  static const int npl = 2;
  double mom_array[nmom] = {0.5, 0.8, 1.0, 1.5, 3.0, 5.0};

  //int pl[npl] = {0,7,7};
  //int col[npl] = {1,2,4};

  int pl[npl] = {0,7};
  int col[npl] = {1,2};
  
  TH1F* h_pip[npl][nmom];
  TFile* f_pip[npl][nmom];

  TH1F* h_pim[npl][nmom];
  TFile* f_pim[npl][nmom];
  
  TCanvas *tc_3vs4_pip = new TCanvas("tc_3vs4_pip","tc_3vs4_pip",1400,1000);
  tc_3vs4_pip->Divide(3,2);

  TCanvas *tc_3vs4_pim = new TCanvas("tc_3vs4_pim","tc_3vs4_pim",1400,1000);
  tc_3vs4_pim->Divide(3,2);

  TLegend* tl_pip[nmom];
  TLegend* tl_pim[nmom];

  for (int imom=0; imom<nmom; ++imom) {

    tl_pip[imom] = new TLegend(0.12,0.12,0.75,0.32);
    tl_pip[imom]->SetFillStyle(0);
    tl_pip[imom]->SetFillColor(0);
    //tl_pip[imom]->SetBorderSize(0);
    tl_pip[imom]->SetHeader("      #color[1]{E_{reco}/E_{true}}, G3 vs. G4");

    tl_pim[imom] = new TLegend(0.12,0.12,0.75,0.32);
    tl_pim[imom]->SetFillStyle(0);
    tl_pim[imom]->SetFillColor(0);
    //tl_pim[imom]->SetBorderSize(0);
    tl_pim[imom]->SetHeader("      #color[1]{E_{reco}/E_{true}}, G3 vs. G4");

    for (int ipl=0; ipl<npl; ++ipl) {
      string fname_pip = Form("hists/hist_%s_MOM_%3.1f_pip.root",tags[pl[ipl]].c_str(),mom_array[imom]);
      cout << "Opening hists file " << fname_pip << endl;
      f_pip[ipl][imom] = TFile::Open(fname_pip.c_str());
      h_pip[ipl][imom] = (TH1F*) f_pip[ipl][imom]->Get("h_energy");
      h_pip[ipl][imom]->SetTitle(Form("#pi^{+}, p = %3.1f GeV/c; E_{reco}/E_{true}; Yield",mom_array[imom]));
      h_pip[ipl][imom]->SetLineWidth(ipl==0?6:2);
      h_pip[ipl][imom]->SetLineColor(col[ipl]);
      h_pip[ipl][imom]->GetXaxis()->SetRangeUser(0.0,1.05);
      //if (ipl==0) h_pip[ipl][imom]->Scale(1.3);
      tl_pip[imom]->AddEntry(h_pip[ipl][imom],Form("%s",tags[pl[ipl]].c_str()),"l");
    }

    for (int ipl=0; ipl<npl; ++ipl) {
      string fname_pim = Form("hists/hist_%s_MOM_%3.1f_pim.root",tags[pl[ipl]].c_str(),mom_array[imom]);
      cout << "Opening hists file " << fname_pim << endl;
      f_pim[ipl][imom] = TFile::Open(fname_pim.c_str());
      h_pim[ipl][imom] = (TH1F*) f_pim[ipl][imom]->Get("h_energy");
      h_pim[ipl][imom]->SetTitle(Form("#pi^{-}, p = %3.1f GeV/c; E_{reco}/E_{true}; Yield",mom_array[imom]));
      h_pim[ipl][imom]->SetLineWidth(ipl==0?6:2);
      h_pim[ipl][imom]->SetLineColor(col[ipl]);
      h_pim[ipl][imom]->GetXaxis()->SetRangeUser(0.0,1.05);
      //if (ipl==0) h_pim[ipl][imom]->Scale(1.3);
      tl_pim[imom]->AddEntry(h_pim[ipl][imom],Form("%s",tags[pl[ipl]].c_str()),"l");
    }
    
    tc_3vs4_pip->cd(1+imom);
    gPad->SetLogy();
    for (int ipl=0; ipl<npl; ++ipl) h_pip[ipl][imom]->Draw(ipl==0?"":"same");
    tl_pip[imom]->Draw();

    tc_3vs4_pim->cd(1+imom);
    gPad->SetLogy();
    for (int ipl=0; ipl<npl; ++ipl) h_pim[ipl][imom]->Draw(ipl==0?"":"same");
    tl_pim[imom]->Draw();

  }

  //tc_3vs4_pip->Print("figs/g4_emv_opt_pip.eps");
  //tc_3vs4_pim->Print("figs/g4_emv_opt_pim.eps");

}

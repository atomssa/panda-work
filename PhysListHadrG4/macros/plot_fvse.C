void plot_fvse() {

  gStyle->SetOptStat(0);
  
  vector<string> tags;
  tags.push_back("G3_HADR3_NUCRIN");
  tags.push_back("QGSP_BIC_EMV");
  tags.push_back("QGSP_BERT_EMV");
  tags.push_back("QGSP_BERT_OPTICAL");
  tags.push_back("QGSP_BERT_EMV_OPTICAL");
  tags.push_back("QGSC_BERT_EMV");
  tags.push_back("QGSP_BERT_HP");

  vector<string> taggs;
  taggs.push_back("Geant3");
  taggs.push_back("QGSP_BIC_EMV");
  taggs.push_back("Geant4");
  taggs.push_back("QGSP_BERT_OPTICAL");
  taggs.push_back("QGSP_BERT_EMV_OPTICAL");
  taggs.push_back("QGSC_BERT_EMV");
  taggs.push_back("QGSP_BERT_HP");

  
  static const int nsetup = 2;
  static const int nmom = 6;
  static const int npl = 2;
  double mom_array[nmom] = {0.5, 0.8, 1.0, 1.5, 3.0, 5.0};
  int pl[npl] = {0,2};
  int colp[npl] = {2,4};
  int colm[npl] = {3,9};
  
  TH1F* h_pip[npl][nsetup];
  TFile* f_pip[npl][nsetup];

  TH1F* h_pim[npl][nsetup];
  TFile* f_pim[npl][nsetup];

  TCanvas *tc_fvse = new TCanvas("tc_fvse","tc_fvse",1400,1000);
  tc_fvse->Divide(2,1);

  TLegend* tl[nsetup];
  for (int ipl=0; ipl<npl; ++ipl) {
    tl[ipl] = new TLegend(0.12,0.12,0.55,0.42);
    tl[ipl]->SetFillStyle(0);
    tl[ipl]->SetFillColor(0);
    tl[ipl]->SetHeader(Form("      %s",taggs[pl[ipl]].c_str()));
  }
  for (int is=0; is<nsetup; ++is) {

    
    for (int ipl=0; ipl<npl; ++ipl) {
      string fname_pip = Form("hists/hist%s_%s_MOM_1.0_pip.root",(is==0?"":"_FULLPANDA"),tags[pl[ipl]].c_str());
      cout << "Opening hists file " << fname_pip << endl;
      f_pip[ipl][is] = TFile::Open(fname_pip.c_str());
      h_pip[ipl][is] = (TH1F*) f_pip[ipl][is]->Get("h_energy");
      h_pip[ipl][is]->SetTitle(Form("%s, p = 1.0 GeV/c %s; E_{reco}/E_{true}; Yield", (ipl==0?"G3":"G4")," -- EMCal Only vs. Full PANDA"));
      h_pip[ipl][is]->SetLineWidth(2);
      h_pip[ipl][is]->SetLineColor(colp[0]);
      if (is==1) h_pip[ipl][is]->SetLineStyle(7);
      h_pip[ipl][is]->GetXaxis()->SetRangeUser(0.0,1.05);
      //tl[ipl]->AddEntry(h_pip[ipl][is],Form("#pi^{+} %s",(is==0?"EMCal Only":"Full PANDA")),"l");
    }

    for (int ipl=0; ipl<npl; ++ipl) {
      string fname_pim = Form("hists/hist%s_%s_MOM_1.0_pim.root",(is==0?"":"_FULLPANDA"),tags[pl[ipl]].c_str());
      cout << "Opening hists file " << fname_pim << endl;
      f_pim[ipl][is] = TFile::Open(fname_pim.c_str());
      h_pim[ipl][is] = (TH1F*) f_pim[ipl][is]->Get("h_energy");
      h_pim[ipl][is]->SetTitle(Form("%s, p = 1.0 GeV/c %s; E_{reco}/E_{true}; Yield", (ipl==0?"G3":"G4"), " -- EMCal Only vs. Full PANDA"));
      h_pim[ipl][is]->SetLineWidth(2);
      h_pim[ipl][is]->SetLineColor(4);
      if (is==1) h_pim[ipl][is]->SetLineStyle(7);
      h_pim[ipl][is]->GetXaxis()->SetRangeUser(0.0,1.05);
      //tl[ipl]->AddEntry(h_pim[ipl][is],Form("#pi^{-} %s",(is==0?"EMCal Only":"Full PANDA")),"l");
    }
  }


  for (int ipl=0; ipl<npl; ++ipl) {
    tc_fvse->cd(1+ipl);

    for (int is=0; is<nsetup; ++is) {
      tl[ipl]->AddEntry(h_pim[ipl][is],Form("#pi^{+} %s",(is==0?"EMCal Only":"Full PANDA")),"l");
    }
    for (int is=0; is<nsetup; ++is) {
      tl[ipl]->AddEntry(h_pip[ipl][is],Form("#pi^{-} %s",(is==0?"EMCal Only":"Full PANDA")),"l");
    }


    
    gPad->SetLogy();
    for (int is=0; is<nsetup; ++is) {
      h_pip[ipl][is]->Draw(is==0?"":"SAME");
      h_pim[ipl][is]->Draw("SAME");
    }
    tl[ipl]->Draw();
  }
  
  tc_fvse->Print("figs/fvse.eps");

}

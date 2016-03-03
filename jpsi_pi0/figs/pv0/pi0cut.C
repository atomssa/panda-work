void pi0cut() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";
  gROOT->LoadMacro(Form("%s/figs/pv0/ananote.C",bdir));

  TFile *fsig[nplab], *fbg[nplab];
  TH2F* sig_avg_true_mc[nplab], *sig_avg_true_rec[nplab], *sig_avg_all_rec[nplab], *sig_avg_all_cut[nplab];
  TH2F* bg_avg_true_mc[nplab], *bg_avg_true_rec[nplab], *bg_avg_all_rec[nplab], *bg_avg_all_cut[nplab];
  TH1F* sig_all_minv[nplab], *sig_all_avg_cut_minv[nplab];
  TH1F* sig_true_minv[nplab], *sig_true_avg_cut_minv[nplab];
  TH1F* bg_all_minv[nplab], *bg_all_avg_cut_minv[nplab];
  TH1F* bg_true_minv[nplab], *bg_true_avg_cut_minv[nplab];
  TCanvas *tc_sig_minv[nplab], *tc_sig_avg[nplab];
  TLatex *tl[6][nplab];
  TLegend *legend = new TLegend(0.15,0.35,0.6,0.6);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  TCanvas *tc_sig_avg[nplab],*tc_bg_avg[nplab];
  TCanvas *tc_sig_minv_allplab = new TCanvas("sig_minv_allplab","sig_minv_allplab",1000,1400);
  tc_sig_minv_allplab->Divide(2,3);
  TCanvas *tc_bg_minv_allplab = new TCanvas("bg_minv_allplab","bg_minv_allplab",1000,1400);
  tc_bg_minv_allplab->Divide(2,3);
  ibrem = 1;
  bool bidims = true;
  for (int iplab=0; iplab<nplab; ++iplab) {

    for (int i=0; i<6; ++i) {
      tl[i][iplab] =new TLatex();
      tl[i][iplab]->SetNDC(true);
      tl[i][iplab]->SetLineColor(i==0?2:1);
      tl[i][iplab]->SetTextSize(((i==2)?1.4:(i<4?1.5:1.2))*tl[i][iplab]->GetTextSize());
    }

    //fsig[iplab] = TFile::Open(Form("%s/test/ana/ana_jpsi_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"),plab[iplab]));
    //fbg[iplab] = TFile::Open(Form("%s/test/ana/ana_pip_pim_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"),plab[iplab]));

    fsig[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/ana_pi0jpsi_%s_p%d.root",bdir,(ibrem==0?"raw":"brem"),iplab));
    fbg[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/ana_pi0jpsi10cm_%s_p%d.root",bdir,(ibrem==0?"raw":"brem"),iplab));
    //fbg[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/ana_pi0pipm_%s_p%d.root",bdir,(ibrem==0?"raw":"brem"),iplab));

    // SIG 2D HISTS
    sig_avg_true_mc[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_mc")->Clone(Form("sig_avg_true_mc_p%d", iplab));
    sig_avg_true_rec[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_rec")->Clone(Form("sig_avg_true_rec_p%d", iplab));
    sig_avg_all_rec[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_all_rec")->Clone(Form("sig_avg_all_rec_p%d", iplab));
    sig_avg_all_cut[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_ana")->Clone(Form("sig_avg_all_cut_p%d", iplab));
    // SIG 1D HISTS
    sig_all_minv[iplab] = (TH1F*) fsig[iplab]->Get("gg/h_m_gg_all_rec")->Clone(Form("sig_all_minv_p%d",iplab));
    sig_all_avg_cut_minv[iplab] = (TH1F*) fsig[iplab]->Get("pi0_ana/h_m_gg_ana_0")->Clone(Form("sig_all_minv_p%d",iplab));
    sig_true_minv[iplab] = (TH1F*) fsig[iplab]->Get("gg/h_m_gg_truepi0_rec")->Clone(Form("sig_true_minv_p%d",iplab));
    sig_true_avg_cut_minv[iplab] = (TH1F*) fsig[iplab]->Get("gg/h_m_gg_pm_ana")->Clone(Form("sig_true_avg_cut_minv_p%d",iplab));
    // SIG 2D HISTS
    bg_avg_true_mc[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_mc")->Clone(Form("sig_avg_true_mc_p%d", iplab));
    bg_avg_true_rec[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_rec")->Clone(Form("sig_avg_true_rec_p%d", iplab));
    bg_avg_all_rec[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_all_rec")->Clone(Form("sig_avg_all_rec_p%d", iplab));
    bg_avg_all_cut[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_ana")->Clone(Form("sig_avg_all_cut_p%d", iplab));
    // SIG 1D HISTS
    bg_all_minv[iplab] = (TH1F*) fbg[iplab]->Get("gg/h_m_gg_all_rec")->Clone(Form("sig_all_minv_p%d",iplab));
    bg_all_avg_cut_minv[iplab] = (TH1F*) fbg[iplab]->Get("pi0_ana/h_m_gg_ana_0")->Clone(Form("sig_all_minv_p%d",iplab));
    bg_true_minv[iplab] = (TH1F*) fbg[iplab]->Get("gg/h_m_gg_truepi0_rec")->Clone(Form("bg_true_minv_p%d",iplab));
    bg_true_avg_cut_minv[iplab] = (TH1F*) fbg[iplab]->Get("gg/h_m_gg_pm_ana")->Clone(Form("bg_true_avg_cut_minv_p%d",iplab));

    if (bidims) {
      tc_sig_avg[iplab] = new TCanvas(Form("sig_avg_p%d",iplab),Form("sig_avg_p%d",iplab),1000,1000);
      tc_sig_avg[iplab]->Divide(2,2);
      tc_sig_avg[iplab]->cd(1);
      sig_avg_true_mc[iplab]->SetTitle(";OA[rad];Avg(E_{#gamma 1},E_{#gamma 2})[GeV]");
      sig_avg_true_mc[iplab]->Draw("colz");
      tl[4][iplab]->DrawLatex(0.15,0.84,Form("Avg(E_{#gamma 1},E_{#gamma 2}) vs OA"));
      tl[5][iplab]->DrawLatex(0.15,0.91,Form("#gamma-#gamma pairs from true #pi^{0}(MC)"));
      tc_sig_avg[iplab]->cd(2);
      sig_avg_true_rec[iplab]->SetTitle(";OA[rad];Avg(E_{#gamma 1},E_{#gamma 2})[GeV]");
      sig_avg_true_rec[iplab]->Draw("colz");
      tl[5][iplab]->DrawLatex(0.15,0.91,Form("#gamma-#gamma pairs from true #pi^{0}(RECO)"));
      tc_sig_avg[iplab]->cd(3);
      sig_avg_all_rec[iplab]->SetTitle(";OA[rad];Avg(E_{#gamma 1},E_{#gamma 2})[GeV]");
      sig_avg_all_rec[iplab]->Draw("colz");
      tl[5][iplab]->DrawLatex(0.22,0.91,Form("All #gamma-#gamma pairs(RECO)"));
      gPad->SetLogz();
      tc_sig_avg[iplab]->cd(4);
      sig_avg_all_cut[iplab]->SetTitle(";OA[rad];Avg(E_{#gamma 1},E_{#gamma 2})[GeV]");
      sig_avg_all_cut[iplab]->Draw("colz");
      tl[5][iplab]->DrawLatex(0.2,0.91,Form("#gamma-#gamma pairs after cut(RECO)"));

      //tc_sig_avg[iplab]->Print(Form("%s/figs/2015.09.15/pi0cut_%s.png",bdir,tc_sig_avg[iplab]->GetName()));
      //gSystem->Exec(Form("convert figs/2015.09.15/pi0cut_%s.png figs/2015.09.15/pi0cut_%s.pdf", tc_sig_avg[iplab]->GetName(), tc_sig_avg[iplab]->GetName()));

      tc_bg_avg[iplab] = new TCanvas(Form("bg_avg_p%d",iplab),Form("bg_avg_p%d",iplab),1000,1000);
      tc_bg_avg[iplab]->Divide(2,2);
      tc_bg_avg[iplab]->cd(1);
      bg_avg_true_mc[iplab]->Draw("colz");
      tc_bg_avg[iplab]->cd(2);
      bg_avg_true_rec[iplab]->Draw("colz");
      tc_bg_avg[iplab]->cd(3);
      gPad->SetLogz();
      bg_avg_all_rec[iplab]->Draw("colz");
      tc_bg_avg[iplab]->cd(4);
      bg_avg_all_cut[iplab]->Draw("colz");
      //tc_bg_avg[iplab]->Print(Form("%s/figs/2015.09.15/pi0cut_%s.png",bdir,tc_bg_avg[iplab]->GetName()));
      //gSystem->Exec(Form("convert %s/figs/2015.09.15/pi0cut_%s.png %s/figs/2015.09.15/pi0cut_%s.pdf", bdir, tc_bg_avg[iplab]->GetName(), bdir, tc_bg_avg[iplab]->GetName()));

    }

    set_style(sig_all_minv[iplab], 4);
    set_style(sig_all_avg_cut_minv[iplab], 2);
    set_style(sig_true_minv[iplab], 4);
    set_style(sig_true_avg_cut_minv[iplab], 2);
    sig_all_minv[iplab]->SetMinimum(0.0);
    sig_all_minv[iplab]->GetXaxis()->SetRangeUser(0.08,0.17);
    sig_true_minv[iplab]->SetMinimum(0.0);
    sig_true_minv[iplab]->GetXaxis()->SetRangeUser(0.08,0.17);
    tc_sig_minv_allplab->cd(1+2*iplab);
    sig_all_minv[iplab]->Draw();
    sig_all_avg_cut_minv[iplab]->Draw("same");
    tl[0][iplab]->DrawLatex(iplab==0?0.15:0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[1][iplab]->DrawLatex(0.15,0.18,Form("All #gamma-#gamma pairs (#bar{p}p#rightarrow#pi^{0}J/#psi)"));
    tc_sig_minv_allplab->cd(1+2*iplab+1);
    tc_sig_minv_allplab->cd(2+2*iplab);
    sig_true_minv[iplab]->Draw();
    sig_true_avg_cut_minv[iplab]->Draw("same");
    tl[2][iplab]->DrawLatex(0.15,0.8,Form("Full Truth Matched #pi^{0}s"));
    tl[3][iplab]->DrawLatex(0.25,0.7,Form("#bar{p}p#rightarrow#pi^{0}J/#psi"));
    if (iplab==0) {
      legend->AddEntry(sig_true_minv[iplab],"Before Cut");
      legend->AddEntry(sig_true_avg_cut_minv[iplab],"After Cut");
      legend->Draw();
    }

    set_style(bg_all_minv[iplab], 4);
    set_style(bg_all_avg_cut_minv[iplab], 2);
    set_style(bg_true_minv[iplab], 4);
    set_style(bg_true_avg_cut_minv[iplab], 2);
    bg_all_minv[iplab]->SetMinimum(0.0);
    bg_all_minv[iplab]->GetXaxis()->SetRangeUser(0.08,0.17);
    bg_true_minv[iplab]->SetMinimum(0.0);
    bg_true_minv[iplab]->GetXaxis()->SetRangeUser(0.08,0.17);
    tc_bg_minv_allplab->cd(1+2*iplab);
    bg_all_minv[iplab]->Draw();
    bg_all_avg_cut_minv[iplab]->Draw("same");
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[1][iplab]->DrawLatex(0.15,0.18,Form("All #gamma-#gamma pairs (#bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-})"));
    tc_bg_minv_allplab->cd(1+2*iplab+1);
    tc_bg_minv_allplab->cd(2+2*iplab);
    bg_true_minv[iplab]->Draw();
    bg_true_avg_cut_minv[iplab]->Draw("same");
    tl[2][iplab]->DrawLatex(0.15,0.8,Form("Full Truth Matched #pi^{0}s"));
    tl[3][iplab]->DrawLatex(0.25,0.7,Form("#bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-}"));
    if (iplab==0) {
      legend->Draw();
    }
  }

  //tc_sig_minv_allplab->Print(Form("%s/figs/2015.09.15/pi0cut_%s.pdf",bdir,tc_sig_minv_allplab->GetName()));
  ///tc_bg_minv_allplab->Print(Form("%s/figs/2015.09.15/pi0cut_%s.pdf",bdir,tc_bg_minv_allplab->GetName()));

}

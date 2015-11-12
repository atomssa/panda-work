void gen_dists() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  TFile *fsig[nplab], *fbg[nplab];
  TH1F* h_ttrue_sig[nplab], *h_ttrue_bg[nplab];
  TH1F* h_trec_sig[nplab], *h_trec_bg[nplab];

  TH1F* h_cthcm_true_sig[nplab], *h_cthcm_true_bg[nplab];
  TH1F* h_thlab_true_sig[nplab], *h_thlab_true_bg[nplab];

  TLatex *tl[6][nplab];
  for (int iplab=0; iplab<nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("%s/hists/note.aug.2015/anav2_jpsi_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));
    fbg[iplab] = TFile::Open(Form("%s/hists/note.aug.2015/anav2_pip_pim_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));
    h_ttrue_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/httrumc")->Clone(Form("httruemc_sig_p%d",iplab));
    h_ttrue_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/httrumc")->Clone(Form("httruemc_bg_p%d",iplab));
    h_trec_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrecgg")->Clone(Form("htrecgg_sig_p%d",iplab));
    h_trec_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrecgg")->Clone(Form("htrecgg_bg_p%d",iplab));
    set_style(h_ttrue_sig[iplab], 1, 5);
    set_style(h_ttrue_bg[iplab], 1, 5);
    set_style(h_trec_sig[iplab], 2, 5);
    set_style(h_trec_bg[iplab], 2, 5);

    h_cthcm_true_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrupi0costhcm")->Clone(Form("h_cthcm_sig_p%d",iplab));
    h_thlab_true_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrupi0thlab")->Clone(Form("h_thlab_sig_p%d",iplab));
    h_cthcm_true_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrupi0costhcm_mc")->Clone(Form("h_cthcm_bg_p%d",iplab));
    h_thlab_true_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrupi0thlab_mc")->Clone(Form("h_thlab_bg_p%d",iplab));
    set_style(h_cthcm_true_sig[iplab], 1, 5);
    set_style(h_thlab_true_sig[iplab], 1, 5);
    set_style(h_cthcm_true_bg[iplab], 1, 5);
    set_style(h_thlab_true_bg[iplab], 1, 5);

    h_cthcm_true_sig[iplab]->SetTitle(";cos{#theta_{CM}}");
    h_thlab_true_sig[iplab]->SetTitle(";#theta_{LAB}");
    h_cthcm_true_bg[iplab]->SetTitle(";cos{#theta_{CM}}");
    h_thlab_true_bg[iplab]->SetTitle(";#theta_{LAB}");

    for (int i=0; i<6; ++i) {
      tl[i][iplab] =new TLatex();
      tl[i][iplab]->SetNDC(true);
      tl[i][iplab]->SetLineColor(i==0?2:1);
      tl[i][iplab]->SetTextSize(((i==2)?1.4:1.5)*tl[i][iplab]->GetTextSize());
    }
  }

  TCanvas *tc_gen_t_sig = new TCanvas("gen_t_dist_sig","gen_t_dist_sig",1600,800);
  tc_gen_t_sig->Divide(3);
  TCanvas *tc_gen_t_bg = new TCanvas("gen_t_dist_bg","gen_t_dist_bg",1600,800);
  tc_gen_t_bg->Divide(3);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tc_gen_t_sig->cd(1+iplab);
    h_ttrue_sig[iplab]->SetTitle(";t[GeV^{2}]");
    h_ttrue_sig[iplab]->Draw();
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tc_gen_t_bg->cd(1+iplab);
    h_ttrue_bg[iplab]->SetTitle(";t[GeV^{2}]");
    h_ttrue_bg[iplab]->Draw();
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  }
  tc_gen_t_sig->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_sig->GetName()));
  tc_gen_t_bg->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_bg->GetName()));


  TCanvas *tc_gen_cthcm_sig = new TCanvas("gen_cthcm_dist_sig","gen_cthcm_dist_sig",1600,800);
  tc_gen_cthcm_sig->Divide(3);
  TCanvas *tc_gen_cthcm_bg = new TCanvas("gen_cthcm_dist_bg","gen_cthcm_dist_bg",1600,800);
  tc_gen_cthcm_bg->Divide(3);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tc_gen_cthcm_sig->cd(1+iplab);
    h_cthcm_true_sig[iplab]->SetTitle(";cos(#theta_{CM})");
    h_cthcm_true_sig[iplab]->Draw();
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tc_gen_cthcm_bg->cd(1+iplab);
    h_cthcm_true_bg[iplab]->SetTitle(";cos(#theta_{CM})");
    h_cthcm_true_bg[iplab]->Draw();
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  }
  tc_gen_cthcm_sig->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_cthcm_sig->GetName()));
  tc_gen_cthcm_bg->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_cthcm_bg->GetName()));

  TCanvas *tc_gen_thlab_sig = new TCanvas("gen_thlab_dist_sig","gen_thlab_dist_sig",1600,800);
  tc_gen_thlab_sig->Divide(3);
  TCanvas *tc_gen_thlab_bg = new TCanvas("gen_thlab_dist_bg","gen_thlab_dist_bg",1600,800);
  tc_gen_thlab_bg->Divide(3);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tc_gen_thlab_sig->cd(1+iplab);
    h_thlab_true_sig[iplab]->SetTitle(";#theta_{LAB}[rad]");
    h_thlab_true_sig[iplab]->Draw();
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tc_gen_thlab_bg->cd(1+iplab);
    h_thlab_true_bg[iplab]->SetTitle(";#theta_{LAB}[rad]");
    h_thlab_true_bg[iplab]->Draw();
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  }
  tc_gen_thlab_sig->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_thlab_sig->GetName()));
  tc_gen_thlab_bg->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_thlab_bg->GetName()));

  TCanvas *tc_gen_t_sig2 = new TCanvas("gen_rec_t_dist_sig","gen_rec_t_dist_sig",1600,800);
  tc_gen_t_sig2->Divide(3);
  TCanvas *tc_gen_t_bg2 = new TCanvas("gen_rec_t_dist_bg2","gen_rec_t_dist_bg",1600,800);
  tc_gen_t_bg2->Divide(3);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tc_gen_t_sig2->cd(1+iplab);
    h_ttrue_sig[iplab]->Draw();
    h_trec_sig[iplab]->Draw("same");
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tc_gen_t_bg2->cd(1+iplab);
    h_ttrue_bg[iplab]->Draw();
    h_trec_bg[iplab]->Draw("same");
    tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  }
  tc_gen_t_sig2->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_sig2->GetName()));
  tc_gen_t_bg2->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_bg2->GetName()));

  //TCanvas *tc_gen_rec_t_sig = new TCanvas("gen_rec_t_dist_sig","gen_rec_t_dist_sig");
  //tc_gen_rec_t_sig->Divide(3);
  //TCanvas *tc_gen_rec_t_bg = new TCanvas("gen_rec_t_dist_bg","gen_rec_t_dist_bg");
  //tc_gen_rec_t_sig->Divide(3);
  //for (int iplab = 0; iplab < nplab; ++iplab) {
  //  cout << iplab << endl;
  //  tc_gen_rec_t_sig->cd(1+iplab);
  //  h_ttrue_sig[iplab]->Draw();
  //  //h_trec_sig[iplab]->SetLineStyle(9);
  //  h_trec_sig[iplab]->Draw("same");
  //  tc_gen_rec_t_bg->cd(1+iplab);
  //  h_ttrue_bg[iplab]->Draw();
  //  //h_trec_bg[iplab]->SetLineStyle(9);
  //  h_trec_bg[iplab]->Draw("same");
  //}

}

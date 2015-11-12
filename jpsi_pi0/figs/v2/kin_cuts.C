void kin_cuts() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  TFile *fsig[nplab], *fbg[nplab];
  TH2F* cmdth_vs_mmtot_true_sig[nplab];
  TH2F* cmdth_vs_cmdph_true_sig[nplab];
  TH2F* cmdth_vs_mmtot_true_bg[nplab];
  TH2F* cmdth_vs_cmdph_true_bg[nplab];
  TH1F* cmdth_sig[nplab];
  TH1F* cmdph_sig[nplab];
  TH1F* mmtot_sig[nplab];
  TH1F* cmdth_bg[nplab];
  TH1F* cmdph_bg[nplab];
  TH1F* mmtot_bg[nplab];
  TLegend *tl = new TLegend(0.15,0.7,0.6,0.85);
  tl->SetTextSize(0.06);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("%s/test/ana/ana_jpsi_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"),plab[iplab]));
    cmdth_vs_mmtot_true_sig[iplab] = (TH2F*) fsig[iplab]->Get("full_sys/h_dth_vs_mass_gg_epair_truepi0jpsi_rec")->Clone(Form("cmdth_vs_mtot_true_sig_iplab%d",iplab));
    cmdth_vs_cmdph_true_sig[iplab] = (TH2F*) fsig[iplab]->Get("full_sys/h_dth_vs_dph_gg_epair_truepi0jpsi_rec")->Clone(Form("cmdth_vs_cmdph_true_sig_iplab%d",iplab));
    //set_style(cmdth_vs_mmtot_true_sig[iplab]);
    //set_style(cmdth_vs_cmdph_true_sig[iplab]);
    cmdth_sig[iplab] =(TH1F*) cmdth_vs_cmdph_true_sig[iplab]->ProjectionX()->Clone(Form("cmdth_sig_iplab%d",iplab));
    cmdph_sig[iplab] =(TH1F*) cmdth_vs_cmdph_true_sig[iplab]->ProjectionY()->Clone(Form("cmdph_sig_iplab%d",iplab));
    mmtot_sig[iplab] =(TH1F*) cmdth_vs_mmtot_true_sig[iplab]->ProjectionX()->Clone(Form("mmtot_sig_iplab%d",iplab));
    set_style(cmdth_sig[iplab],1);
    set_style(cmdph_sig[iplab],1);
    set_style(mmtot_sig[iplab],1);

    fbg[iplab] = TFile::Open(Form("%s/test/ana/ana_pip_pim_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"),plab[iplab]));
    cmdth_vs_mmtot_true_bg[iplab] = (TH2F*) fbg[iplab]->Get("full_sys/h_dth_vs_mass_gg_epair_truepi0jpsi_rec")->Clone(Form("cmdth_vs_mtot_true_bg_iplab%d",iplab));
    cmdth_vs_cmdph_true_bg[iplab] = (TH2F*) fbg[iplab]->Get("full_sys/h_dth_vs_dph_gg_epair_truepi0jpsi_rec")->Clone(Form("cmdth_vs_cmdph_true_bg_iplab%d",iplab));
    //set_style(cmdth_vs_mmtot_true_bg[iplab]);
    //set_style(cmdth_vs_cmdph_true_bg[iplab]);
    cmdth_bg[iplab] =(TH1F*) cmdth_vs_cmdph_true_bg[iplab]->ProjectionX()->Clone(Form("cmdth_bg_iplab%d",iplab));
    cmdph_bg[iplab] =(TH1F*) cmdth_vs_cmdph_true_bg[iplab]->ProjectionY()->Clone(Form("cmdph_bg_iplab%d",iplab));
    mmtot_bg[iplab] =(TH1F*) cmdth_vs_mmtot_true_bg[iplab]->ProjectionX()->Clone(Form("mmtot_bg_iplab%d",iplab));
    set_style(cmdth_bg[iplab],2);
    set_style(cmdph_bg[iplab],2);
    set_style(mmtot_bg[iplab],2);

    if (iplab==0) {
      tl->AddEntry(cmdth_sig[iplab],"Signal");
      tl->AddEntry(cmdth_bg[iplab],"Background");
    }
  }

  double mmin = 3.0;
  double dth_min = 2.4;
  double dth_max = 3.7;
  double dph_min = 2.6;
  double dph_max = 3.7;

  TCanvas *tc1  = new TCanvas("cm_kin1","cm_kin1",1400,800);
  tc1->Divide(3,1);
  tc1->cd(1);
  mmtot_sig[0]->Draw();
  //mmtot_bg[0]->Draw("same");
  //tl->Draw();
  TLine *tline = new TLine();
  tline->SetLineColor(2);
  tline->SetLineWidth(2);
  tline->DrawLine(mmin,0,mmin,100);
  tc1->cd(2);
  cmdth_sig[0]->Draw();
  tline->DrawLine(dth_min,0,dth_min,200);
  tline->DrawLine(dth_max,0,dth_max,200);
  //cmdth_bg[0]->Draw("same");
  tc1->cd(3);
  cmdph_sig[0]->Draw();
  tline->DrawLine(dph_min,0,dph_min,200);
  tline->DrawLine(dph_max,0,dph_max,200);
  //cmdph_bg[0]->Draw("same");
  tc1->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc1->GetName()));

}

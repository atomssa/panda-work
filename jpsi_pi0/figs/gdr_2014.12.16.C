#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1F.h"

static int ibrem = 1;
static const int nplab = 3;
static const double plab[nplab] = {5.513, 8., 12.};
// these are numbers derived from analytical formula for
// values of t at cos_theta_cm = {-1, 0, 1}
static double limits[nplab][3] =
  { {-1.50406, -0.443789, 0.616486},
    {-5.96462, -2.75368, 0.457248},
    {-13.2926, -6.48861, 0.31538} };
static double tmin[nplab]={-0.443789, -0.5, -0.5};
static double tmax[nplab]={0.616486, 0.457248, 0.31538};

void set_style(TH1* h, int col=4, int rebin=1, bool sumw2=false);
void set_style_ana(TH1* h, int col=4, int rebin=1, bool sumw2=false);
void set_style(TGraph *g, int col=4, int sty=20, int s=2, int w=2);

void fig_num_evt() {
  double sig[nplab] = {28040.,23619,15460};
  double bg[nplab] = {36000, 9000, 4500};

  TLatex *tl[6][nplab];
  TLegend *legend = new TLegend(0.3,0.6,0.99,0.8);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    for (int i=0; i<6; ++i) {
      tl[i][iplab] =new TLatex();
      //tl[i][iplab]->SetNDC(true);
      tl[i][iplab]->SetLineColor(i==0?2:1);
      //tl[i][iplab]->SetTextSize((i==1?1.1:1.4)*tl[i][iplab]->GetTextSize());
      //tl[i][iplab]->SetTextSize(((i==2)?1.4:1.5)*tl[i][iplab]->GetTextSize());
    }
  }

  TH1F* hdummy = new TH1F("my_dummy","my_dummy",10,0,15);
  set_style(hdummy,1);

  TGraph *tgsig = new TGraph(nplab, plab, sig);
  tgsig->SetTitle("Number of events in full signal MC;p_{lab}[GeV/c];num. evt");
  tgsig->SetMarkerStyle(20);
  tgsig->SetMarkerSize(2);
  tgsig->SetMarkerColor(2);
  tgsig->SetMarkerColor(2);
  TCanvas *tcsig = new TCanvas("num_evt_sig","num_evt_sig");
  tcsig->cd();
  hdummy->Draw();
  tgsig->Draw("ap");
  tgsig->SetMinimum(0);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tl[0][iplab]->DrawLatex(plab[iplab]+(iplab==2?-2.6:0.3),sig[iplab]*0.97,Form("%4.2f<|t/u|<%4.2f",tmin[iplab],tmax[iplab]));
  }
  tcsig->Print(Form("figs/2014.12.16/%s.pdf", tcsig->GetName()));

  TGraph *tgbg = new TGraph(nplab, plab, bg);
  tgbg->SetTitle("Number of events in full background MC;p_{lab}[GeV/c];num. evt");
  tgbg->SetMarkerStyle(24);
  tgbg->SetMarkerSize(2);
  tgbg->SetMarkerColor(2);
  TCanvas *tcbg = new TCanvas("num_evt_bg","num_evt_bg");
  tcbg->cd();
  hdummy->Draw();
  tgbg->Draw("ap");
  tgbg->SetMinimum(2);
  //gPad->SetLogy();
  //TGraph *tmg =  new TMultiGraph();
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tl[1][iplab]->DrawLatex(plab[iplab]+(iplab==2||iplab==1?-2.6:0.3),bg[iplab]*0.97,Form("%4.2f<|t/u|<%4.2f",tmin[iplab],tmax[iplab]));
  }
  tcbg->Print(Form("figs/2014.12.16/%s.pdf", tcbg->GetName()));

}

void fig_kin_cuts() {
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
    fsig[iplab] = TFile::Open(Form("test/ana/ana_jpsi_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
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

    fbg[iplab] = TFile::Open(Form("test/ana/ana_pip_pim_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
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
  tc1->Print(Form("figs/2014.12.16/%s.pdf",tc1->GetName()));

}

//void fig_kin_cuts2() {
//  TFile *fsig[nplab], *fbg[nplab];
//  for (int iplab = 0; iplab < nplab; ++iplab) {
//    fsig[iplab] = TFile::Open(Form("test/anav2_jpsi_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
//    fbg[iplab] = TFile::Open(Form("test/anav2_pip_pim_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
//    for (int istep = 0; istep < nstep; ++istep) {
//      sig_nep[iplab][istep] = (TH1F*) fsig[iplab]->Get(Form("hnep_%d",istep))->Clone(Form("sig_hnep_s%d_p%d",istep,iplab));
//      bg_nep[iplab][istep] = (TH1F*) fbg[iplab]->Get(Form("hnep_%d",istep))->Clone(Form("bg_hnep_s%d_p%d",istep,iplab));
//    }
//  }
//}

void fig_npair() {
  static const int nstep = 5;
  TH1F* sig_nep[nplab][nstep], *sig_ngg[nplab][nstep];
  TH1F* bg_nep[nplab][nstep], *bg_ngg[nplab][nstep];
  TFile *fsig[nplab], *fbg[nplab];
  for (int iplab = 0; iplab < nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("test/anav2_jpsi_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
    fbg[iplab] = TFile::Open(Form("test/anav2_pip_pim_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
    for (int istep = 0; istep < nstep; ++istep) {
      sig_nep[iplab][istep] = (TH1F*) fsig[iplab]->Get(Form("hnep_%d",istep))->Clone(Form("sig_hnep_s%d_p%d",istep,iplab));
      bg_nep[iplab][istep] = (TH1F*) fbg[iplab]->Get(Form("hnep_%d",istep))->Clone(Form("bg_hnep_s%d_p%d",istep,iplab));
    }
  }

  TArrow *ar[nplab];
  TLatex *tl[6][nplab];
  TLegend *legend = new TLegend(0.3,0.6,0.99,0.8);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    for (int i=0; i<6; ++i) {
      tl[i][iplab] =new TLatex();
      tl[i][iplab]->SetNDC(true);
      tl[i][iplab]->SetLineColor(i==0?2:1);
      tl[i][iplab]->SetTextSize((i==1?1.1:1.4)*tl[i][iplab]->GetTextSize());
      tl[i][iplab]->SetTextSize(((i==2)?1.4:1.5)*tl[i][iplab]->GetTextSize());
    }
    set_style(sig_nep[iplab][0],1);
    set_style(sig_nep[iplab][1],2);
    sig_nep[iplab][iplab==0?0:1]->SetTitle(";N_{e^{+}e^{-}}/Event");
  }

  TCanvas *sig_s0_s1 = new TCanvas("sig_s0_s1","sig_s0_s1");
  sig_s0_s1->Divide(3,1);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    sig_s0_s1->cd(1+iplab);
    sig_nep[iplab][iplab==0?0:1]->Draw();
    sig_nep[iplab][iplab==0?1:0]->Draw("same");
    sig_nep[iplab][1]->Draw("same");
    if (true) {
      tl[0][iplab]->DrawLatex(0.45,0.35,"N_{evt}(N_{e^{+}e^{-}}>1)");
      ar[iplab] = new TArrow(3.5,iplab==0?4500:(iplab==1?8000:9000),2.5,500,0.02);
      ar[iplab]->SetLineWidth(2);
      ar[iplab]->Draw();
    }
    tl[1][iplab]->DrawLatex(0.25,0.93,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    if (iplab==0) {
      legend->SetTextSize(0.09);
      legend->AddEntry(sig_nep[iplab][0],"#color[1]{Before EID}","");
      legend->AddEntry(sig_nep[iplab][1],"#color[2]{After EID}","");
      legend->Draw();
    }
  }

  TCanvas *sig_s0_s1_noarrow = new TCanvas("sig_s0_s1_noarrow","sig_s0_s1_noarrow");
  sig_s0_s1_noarrow->Divide(3);
  for (int iplab = 0; iplab < nplab; ++iplab) {
    sig_s0_s1_noarrow->cd(1+iplab);
    sig_nep[iplab][iplab==0?0:1]->Draw();
    sig_nep[iplab][iplab==0?1:0]->Draw("same");
    sig_nep[iplab][1]->Draw("same");
    tl[1][iplab]->DrawLatex(0.25,0.93,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    if (iplab==0) {
      legend->Draw();
    }
  }
  sig_s0_s1->Print("figs/2014.12.16/num_ep_pair_vs_cut_sig.pdf");
  sig_s0_s1_noarrow->Print("figs/2014.12.16/num_ep_pair_vs_cut_sig_noarrow.pdf");
}

void fig_pi0cut() {

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

    fsig[iplab] = TFile::Open(Form("test/ana/ana_jpsi_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
    fbg[iplab] = TFile::Open(Form("test/ana/ana_pip_pim_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"),plab[iplab]));
    sig_avg_true_mc[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_mc")->Clone(Form("sig_avg_true_mc_p%d", iplab));
    sig_avg_true_rec[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_rec")->Clone(Form("sig_avg_true_rec_p%d", iplab));
    sig_avg_all_rec[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_all_rec")->Clone(Form("sig_avg_all_rec_p%d", iplab));
    sig_avg_all_cut[iplab]= (TH2F*) fsig[iplab]->Get("gg/h_oa_gg_avg_e_g_ana")->Clone(Form("sig_avg_all_cut_p%d", iplab));
    sig_all_minv[iplab] = (TH1F*) fsig[iplab]->Get("gg/h_m_gg_all_rec")->Clone(Form("sig_all_minv_p%d",iplab));
    sig_all_avg_cut_minv[iplab] = (TH1F*) fsig[iplab]->Get("pi0_ana/h_m_gg_ana_0")->Clone(Form("sig_all_minv_p%d",iplab));
    sig_true_minv[iplab] = (TH1F*) fsig[iplab]->Get("gg/h_m_gg_truepi0_rec")->Clone(Form("sig_true_minv_p%d",iplab));
    sig_true_avg_cut_minv[iplab] = (TH1F*) fsig[iplab]->Get("gg/h_m_gg_pm_ana")->Clone(Form("sig_true_avg_cut_minv_p%d",iplab));

    bg_avg_true_mc[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_mc")->Clone(Form("sig_avg_true_mc_p%d", iplab));
    bg_avg_true_rec[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_truepi0_rec")->Clone(Form("sig_avg_true_rec_p%d", iplab));
    bg_avg_all_rec[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_all_rec")->Clone(Form("sig_avg_all_rec_p%d", iplab));
    bg_avg_all_cut[iplab]= (TH2F*) fbg[iplab]->Get("gg/h_oa_gg_avg_e_g_ana")->Clone(Form("sig_avg_all_cut_p%d", iplab));
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
      tc_sig_avg[iplab]->Print(Form("figs/2014.12.16/pi0cut_%s.png",tc_sig_avg[iplab]->GetName()));
      gSystem->Exec(Form("convert figs/2014.12.16/pi0cut_%s.png figs/2014.12.16/pi0cut_%s.pdf", tc_sig_avg[iplab]->GetName(), tc_sig_avg[iplab]->GetName()));

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
      tc_bg_avg[iplab]->Print(Form("figs/2014.12.16/pi0cut_%s.png",tc_bg_avg[iplab]->GetName()));
      gSystem->Exec(Form("convert figs/2014.12.16/pi0cut_%s.png figs/2014.12.16/pi0cut_%s.pdf", tc_bg_avg[iplab]->GetName(), tc_bg_avg[iplab]->GetName()));

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

  tc_sig_minv_allplab->Print(Form("figs/2014.12.16/pi0cut_%s.pdf",tc_sig_minv_allplab->GetName()));
  tc_bg_minv_allplab->Print(Form("figs/2014.12.16/pi0cut_%s.pdf",tc_bg_minv_allplab->GetName()));

}

void draw_fg_sig_bg(TCanvas *canv, int ipad, TH1F*fg, TH1F*sig, TH1F*bg, bool legend=false) {
  canv->cd(ipad);
  fg->Draw("hist");
  sig->Draw("same");
  bg->Draw("same");
  if (legend) {
    TLegend *tl = new TLegend(0.2,0.6, 0.45,0.88);
    tl->AddEntry(fg,"Foreground","l");
    tl->AddEntry(sig,"Signal","pl");
    tl->AddEntry(bg,"Background","pl");
    tl->Draw();
  }
}

void fetch_fg_bg_sig(TFile*fsig, TFile *fbg, const char*dir, const char *stag, const char*dtag, TH1F* fg, TH1F* sig, TH1F*bg){
  fg = (TH1F*) fsig->Get(Form("%shmep%s",dir,stag))->Clone(Form("hmep_fg_%s",dtag));
  sig = (TH1F*) fsig->Get(Form("%shmep%s",dir,stag))->Clone(Form("hmep_sig_%s",dtag));
  bg = (TH1F*) fbg->Get(Form("%shmep%s",dir,stag))->Clone(Form("hmep_bg_%s",dtag));
  fg->Add(bg);
  set_style(fg, 1);
  set_style(sig, 2);
  set_style(bg, 4);
}

void fig_gen_dists() {

  TFile *fsig[nplab], *fbg[nplab];
  TH1F* h_ttrue_sig[nplab], *h_ttrue_bg[nplab];
  TH1F* h_trec_sig[nplab], *h_trec_bg[nplab];

  TH1F* h_cthcm_true_sig[nplab], *h_cthcm_true_bg[nplab];
  TH1F* h_thlab_true_sig[nplab], *h_thlab_true_bg[nplab];

  TLatex *tl[6][nplab];
  for (int iplab=0; iplab<nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("test/anav2_jpsi_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"), plab[iplab]));
    fbg[iplab] = TFile::Open(Form("test/anav2_pip_pim_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"), plab[iplab]));
    h_ttrue_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/httrumc")->Clone(Form("httruemc_sig_p%d",iplab));
    h_ttrue_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/httrumc")->Clone(Form("httruemc_bg_p%d",iplab));
    h_trec_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrecgg")->Clone(Form("htrecgg_sig_p%d",iplab));
    h_trec_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrecgg")->Clone(Form("htrecgg_bg_p%d",iplab));
    set_style(h_ttrue_sig[iplab], 1);
    set_style(h_ttrue_bg[iplab], 1);
    set_style(h_trec_sig[iplab], 2);
    set_style(h_trec_bg[iplab], 2);

    h_cthcm_true_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrupi0costhcm")->Clone(Form("h_cthcm_sig_p%d",iplab));
    h_thlab_true_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrupi0thlab")->Clone(Form("h_thlab_sig_p%d",iplab));
    h_cthcm_true_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrupi0costhcm")->Clone(Form("h_cthcm_bg_p%d",iplab));
    h_thlab_true_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrupi0thlab")->Clone(Form("h_thlab_bg_p%d",iplab));
    set_style(h_cthcm_true_sig[iplab],1);
    set_style(h_thlab_true_sig[iplab],1);
    set_style(h_cthcm_true_bg[iplab],1);
    set_style(h_thlab_true_bg[iplab],1);

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
  tc_gen_t_sig->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_t_sig->GetName()));
  tc_gen_t_bg->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_t_bg->GetName()));


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
  tc_gen_cthcm_sig->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_cthcm_sig->GetName()));
  tc_gen_cthcm_bg->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_cthcm_bg->GetName()));

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
  tc_gen_thlab_sig->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_thlab_sig->GetName()));
  tc_gen_thlab_bg->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_thlab_bg->GetName()));

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
  tc_gen_t_sig2->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_t_sig2->GetName()));
  tc_gen_t_bg2->Print(Form("figs/2014.12.16/%s.pdf", tc_gen_t_bg2->GetName()));

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

void fig_ana() {

  const int nstep = 5;
  const int ntbin = 12;

  double mmin = 2.96, mmax = 3.22;
  int immin = 0, immax = 0;
  double x[nstep] = {1, 2 ,3, 4, 5};

  double yield_sig[nplab][nstep]= {{0.0}}, yield_bg[nplab][nstep]= {{0.0}};
  double eff_sig[nplab][nstep]= {{0.0}}, eff_bg[nplab][nstep]= {{0.0}};
  double stob[nplab][nstep]= {{0.0}};

  TGraph *tg_sig_eff[nplab], *tg_bg_eff[nplab];
  TGraph *tg_sig_yield[nplab], *tg_bg_yield[nplab];
  TGraph *tg_stob[nplab];

  TMultiGraph *tmg_yield[nplab], *tmg_eff[nplab];

  TCanvas *tc_mep_all_steps[nplab];
  TCanvas *tc_mep_steps[nplab][nstep];

  TFile *fsig[nplab], *fbg[nplab];
  TH1F* hmsg[nplab][nstep], *hmbg[nplab][nstep], *hmfg[nplab][nstep];
  TH1F* hmsg_tb[nplab][ntbin], *hmbg_tb[nplab][ntbin], *hmfg_tb[nplab][ntbin];

  TLatex *tl[6][nplab];
  for (int ii = 0; ii < 6; ++ii) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tl[ii][iplab] = new TLatex();
      tl[ii][iplab]->SetNDC(true);
      //tl[ii][iplab]->SetTextColor(ii==0?2:1);
      if (ii==0) tl[ii][iplab]->SetTextSize((iplab==0?1.5:2.0)*tl[ii][iplab]->GetTextSize());
      if (ii==1) tl[ii][iplab]->SetTextSize(0.8*(iplab==0?1.5:2.0)*tl[ii][iplab]->GetTextSize());
      if (ii==2) tl[ii][iplab]->SetTextSize(1.8*tl[ii][iplab]->GetTextSize());
      if (ii==3) tl[ii][iplab]->SetTextSize(0.8*1.8*tl[ii][iplab]->GetTextSize());
    }
  }

  TLegend *leg[3];
  leg[0] = new TLegend(0.2,0.5,0.8,0.7);
  leg[1] = new TLegend(0.15,0.65,0.7,0.75);
  leg[2] = new TLegend(0.2,0.25,0.8,0.45);
  for (int jj = 0; jj < 3; ++jj) {
    leg[jj]->SetFillStyle(0);
    leg[jj]->SetBorderSize(0);
    leg[jj]->SetTextSize(0.07);
  }

  for (int iplab=0; iplab<nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("test/anav2_jpsi_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"), plab[iplab]));
    fbg[iplab] = TFile::Open(Form("test/anav2_pip_pim_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"), plab[iplab]));

    for (int istep=0; istep<nstep; ++istep) {

      hmfg[iplab][istep] = (TH1F*) fsig[iplab]->Get(Form("hmep_%d",istep))->Clone(Form("hmep_fg_p%d_s%d",iplab,istep));
      hmsg[iplab][istep] = (TH1F*) fsig[iplab]->Get(Form("hmep_%d",istep))->Clone(Form("hmep_sig_p%d_s%d",iplab,istep));
      hmbg[iplab][istep] = (TH1F*) fbg[iplab]->Get(Form("hmep_%d",istep))->Clone(Form("hmep_bg_p%d_s%d",iplab,istep));
      hmfg[iplab][istep]->Add(hmbg[iplab][istep]);
      hmfg[iplab][istep]->SetTitle(";M_{inv}[GeV/c^{2}]");
      if (immin==0||immax==0) {
	immin = hmsg[iplab][istep]->GetXaxis()->FindBin(mmin);
	immax = hmsg[iplab][istep]->GetXaxis()->FindBin(mmax);
      }

      // Hackery to combine the pid steps for SIG and BG
      if (istep==0) {
	yield_sig[iplab][istep] = hmsg[iplab][istep]->Integral(immin, immax);
	yield_bg[iplab][istep] = hmbg[iplab][istep]->Integral(immin, immax);
	eff_sig[iplab][istep] = istep==0?100.:100*yield_sig[iplab][istep]/yield_sig[iplab][0];
	eff_bg[iplab][istep] = istep==0?100.:100*yield_bg[iplab][istep]/yield_bg[iplab][0];
	stob[iplab][istep] = yield_sig[iplab][istep]/yield_bg[iplab][istep];
      } else if (istep>1) {
	yield_sig[iplab][istep-1] = hmsg[iplab][istep]->Integral(immin, immax);
	yield_bg[iplab][istep-1] = hmbg[iplab][istep]->Integral(immin, immax);
	eff_sig[iplab][istep-1] = istep==0?100.:100*yield_sig[iplab][istep-1]/yield_sig[iplab][0];
	eff_bg[iplab][istep-1] = istep==0?100.:100*yield_bg[iplab][istep-1]/yield_bg[iplab][0];
	stob[iplab][istep-1] = yield_sig[iplab][istep-1]/yield_bg[iplab][istep-1];
      }

      set_style(hmfg[iplab][istep], 1, 2, true);
      set_style(hmsg[iplab][istep], 2, 2, true);
      set_style(hmbg[iplab][istep], 4, 2, true);

      hmfg[iplab][istep]->GetXaxis()->SetRangeUser(0,4.0);
      hmsg[iplab][istep]->GetXaxis()->SetRangeUser(0,4.0);
      hmbg[iplab][istep]->GetXaxis()->SetRangeUser(0,4.0);

    }

    tc_mep_all_steps[iplab] = new TCanvas(Form("tc_mep_all_steps_p%d",iplab),Form("tc_mep_all_steps_p%d",iplab));
    tc_mep_all_steps[iplab]->Divide(3,2);
    for (int istep=0; istep<nstep; ++istep) {
      tc_mep_steps[iplab][istep] = new TCanvas(Form("tc_mep_steps_p%d_s%d",iplab,istep),Form("tc_mep_steps_p%d_s%d",iplab,istep));
      draw_fg_sig_bg(tc_mep_all_steps[iplab], istep+1, hmfg[iplab][istep], hmsg[iplab][istep], hmbg[iplab][istep], true);
      draw_fg_sig_bg(tc_mep_steps[iplab][istep], 0, hmfg[iplab][istep], hmsg[iplab][istep], hmbg[iplab][istep], true);
      //tc_mep_steps[iplab][istep]->Print(Form("figs/2014.12.16/%s.pdf",tc_mep_steps[iplab][istep]->GetName()));
    }
    //tc_mep_all_steps[iplab]->Print(Form("figs/2014.12.16/%s.pdf",tc_mep_all_steps[iplab]->GetName()));

    tmg_yield[iplab] = new TMultiGraph(Form("tmg_yield_p%d",iplab),Form("tmg_yield_p%d",iplab));
    tg_sig_yield[iplab] = new TGraph(nstep-1,x,yield_sig[iplab]);
    set_style(tg_sig_yield[iplab], 1, 20, 3);
    //set_style(tg_sig_yield[iplab], iplab==0?1:(iplab==1?2:4), 20, 1);
    tg_sig_yield[iplab]->SetMinimum(0);
    tmg_yield[iplab]->Add(tg_sig_yield[iplab],"p");
    tg_bg_yield[iplab] = new TGraph(nstep-1,x,yield_bg[iplab]);
    //set_style(tg_bg_yield[iplab], iplab==0?1:(iplab==1?2:4), 24, 1);
    set_style(tg_bg_yield[iplab], 1, 24, 3);
    tg_bg_yield[iplab]->SetMinimum(0);
    tmg_yield[iplab]->Add(tg_bg_yield[iplab],"p");
    for (int jj = 0; jj < 3; ++jj) {
      if (iplab==0) {
	leg[jj]->SetHeader(jj==0?"Raw yields":(jj==1?"Yields relative to no cuts":"S/B"));
	leg[jj]->AddEntry(tg_sig_yield[iplab],"Signal","p");
	leg[jj]->AddEntry(tg_bg_yield[iplab],"Background","p");
      }
    }

    tmg_eff[iplab] = new TMultiGraph(Form("tmg_eff_p%d",iplab),Form("tmg_eff_p%d",iplab));
    tg_sig_eff[iplab] = new TGraph(nstep-1,x,eff_sig[iplab]);
    //set_style(tg_sig_eff[iplab], iplab==0?1:(iplab==1?2:4), 20, 1);
    set_style(tg_sig_eff[iplab], 1, 20, 2);
    tg_sig_eff[iplab]->SetMinimum(0);
    tmg_eff[iplab]->Add(tg_sig_eff[iplab],"p");
    tg_bg_eff[iplab] = new TGraph(nstep-1,x,eff_bg[iplab]);
    //set_style(tg_bg_eff[iplab], iplab==0?1:(iplab==1?2:4), 24, 1);
    set_style(tg_bg_eff[iplab], 1, 24, 2);
    tg_bg_eff[iplab]->SetMinimum(0);
    tmg_eff[iplab]->Add(tg_bg_eff[iplab],"p");

    tg_stob[iplab] = new TGraph(nstep-1,x,stob[iplab]);
    //set_style(tg_stob[iplab], iplab==0?1:(iplab==1?2:4), 20, 1);
    set_style(tg_stob[iplab], 1, 20, 3);
    tg_stob[iplab]->SetMinimum(0);

  }

  const char* bin_label[5] = {"EID(BG)","EID(SIG)","& N_{e^{+}e^{-}}=1  ","N_{#pi^{0}}>0","BtoB"};

  TCanvas *yield_step = new TCanvas("yield_step","yield_step",1800,1000);
  TPad *pad_yield[nplab];
  TH1F *hdummy_yield[nplab];
  double max_yield = 0;
  for (int iplab = 0; iplab < nplab; ++iplab) { if( yield_sig[iplab][0] > max_yield) max_yield = yield_sig[iplab][0];  }
  for (int iplab=0; iplab<nplab; ++iplab) {
    double xl = (iplab==0?0:0.1)+iplab*0.3;
    double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
    cout << "xl= " << xl << " xh= " << xh << endl;
    pad_yield[iplab] = new TPad(Form("pad_yield%d",iplab),Form("pad_yield%d",iplab),xl,0.0,xh,1.0);
    yield_step->cd(0);
    pad_yield[iplab]->Draw();
    double epsilon=1e-9;
    if (iplab==0) {pad_yield[iplab]->SetRightMargin(epsilon); pad_yield[iplab]->SetLeftMargin(0.2);}
    if (iplab==1) {pad_yield[iplab]->SetLeftMargin(epsilon); pad_yield[iplab]->SetRightMargin(epsilon); }
    if (iplab==2) {pad_yield[iplab]->SetLeftMargin(epsilon); pad_yield[iplab]->SetRightMargin(0.1);}
    pad_yield[iplab]->SetTopMargin(0.05);
    pad_yield[iplab]->SetBottomMargin(0.15);
    pad_yield[iplab]->SetTicks(0,1);
    pad_yield[iplab]->cd();
    tmg_yield[iplab]->Draw("a");
    hdummy_yield[iplab] = tmg_yield[iplab]->GetHistogram();
    for (int ibin=0; ibin<5; ++ibin){
      int iibin = 0;
      if (ibin==0) iibin = 5;
      else if (ibin==1) iibin = 32;
      else if (ibin==2) iibin = 39;
      else iibin = 36+(ibin-2)*30;
      hdummy_yield[iplab]->GetXaxis()->SetBinLabel(iibin, bin_label[ibin]);
    }
    hdummy_yield[iplab]->SetLabelSize(0.05,"Y");
    hdummy_yield[iplab]->SetLabelSize(iplab==0?0.065:0.08,"X");
    hdummy_yield[iplab]->SetLabelOffset(0.005,"Y");
    hdummy_yield[iplab]->SetTitleSize(0.06,"Y");
    hdummy_yield[iplab]->SetTitleOffset(1.7,"Y");
    hdummy_yield[iplab]->SetTitle(";;counts in 2.96<M_{e^{+}e^{-}}[GeV/c^{2}]<3.22");
    hdummy_yield[iplab]->SetMinimum(0);
    hdummy_yield[iplab]->SetMaximum(max_yield*1.1);
    tl[0][iplab]->DrawLatex(iplab==0?0.35:(iplab==1?0.25:0.15),0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[1][iplab]->DrawLatex(iplab==0?0.40:(iplab==1?0.35:0.25),0.75,Form("(%3.2f<|u|<%3.2f)",tmin[iplab],tmax[iplab]));
    if (iplab==2) leg[0]->Draw();
  }
  //yield_step->Print("figs/2014.12.16/yield_vs_steps.pdf");
  yield_step->Print("yield_vs_steps.pdf");

  TCanvas *eff_step = new TCanvas("eff_step","eff_step");
  TPad *pad_eff[nplab];
  TH1F *hdummy_eff[nplab];
  for (int iplab=0; iplab<nplab; ++iplab) {
    double xl = (iplab==0?0:0.1)+iplab*0.3;
    double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
    cout << "xl= " << xl << " xh= " << xh << endl;
    pad_eff[iplab] = new TPad(Form("pad_eff%d",iplab),Form("pad_eff%d",iplab),xl,0.0,xh,1.0);
    eff_step->cd(0);
    pad_eff[iplab]->Draw();
    double epsilon=1e-9;
    if (iplab==0) {pad_eff[iplab]->SetRightMargin(epsilon); pad_eff[iplab]->SetLeftMargin(0.2);}
    if (iplab==1) {pad_eff[iplab]->SetLeftMargin(epsilon); pad_eff[iplab]->SetRightMargin(epsilon); }
    if (iplab==2) {pad_eff[iplab]->SetLeftMargin(epsilon); pad_eff[iplab]->SetRightMargin(0.1);}
    pad_eff[iplab]->SetTicks(0,1);
    pad_eff[iplab]->cd();
    tmg_eff[iplab]->Draw("a");
    hdummy_eff[iplab] = tmg_eff[iplab]->GetHistogram();
    hdummy_eff[iplab]->SetLabelSize(0.05,"Y");
    hdummy_eff[iplab]->SetLabelSize(0.05,"X");
    hdummy_eff[iplab]->SetLabelOffset(0.005,"Y");
    hdummy_eff[iplab]->SetTitleSize(0.06,"Y");
    hdummy_eff[iplab]->SetTitleOffset(1.1,"Y");
    hdummy_eff[iplab]->SetTitle(";step number;efficiency in 2.96<M_{e^{+}e^{-}}[GeV/c^{2}]<3.22");
    hdummy_eff[iplab]->SetMinimum(0);
    hdummy_eff[iplab]->SetMaximum(105);
    tl[0][iplab]->DrawLatex(iplab==0?0.35:(iplab==1?0.25:0.15),0.23,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[1][iplab]->DrawLatex(iplab==0?0.40:(iplab==1?0.35:0.25),0.18,Form("(%3.2f<|u|<%3.2f)",tmin[iplab],tmax[iplab]));
    if (iplab==2) leg[1]->Draw();
  }
  //eff_step->Print("figs/2014.12.16/eff_vs_steps.pdf");

  TCanvas *stob_step = new TCanvas("stob_step","stob_step",2200,1000);
  stob_step->Divide(3,1);
  TH1F *hdummy_stob[nplab];
  for (int iplab=0; iplab<nplab; ++iplab) {
    stob_step->cd(1+iplab);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    tg_stob[iplab]->Draw("ap");
    hdummy_stob[iplab] = tg_stob[iplab]->GetHistogram();
    for (int ibin=0; ibin<5; ++ibin){
      int iibin = 0;
      if (ibin==0) iibin = 5;
      else if (ibin==1) iibin = 32;
      else if (ibin==2) iibin = 39;
      else iibin = 36+(ibin-2)*30;
      hdummy_stob[iplab]->GetXaxis()->SetBinLabel(iibin, bin_label[ibin]);
    }
    hdummy_stob[iplab]->SetLabelSize(0.05,"Y");
    hdummy_stob[iplab]->SetLabelSize(0.07,"X");
    hdummy_stob[iplab]->SetLabelOffset(0.005,"Y");
    hdummy_stob[iplab]->SetTitleSize(0.06,"Y");
    hdummy_stob[iplab]->SetTitleOffset(1.1,"Y");
    hdummy_stob[iplab]->SetTitle(";;");
    hdummy_stob[iplab]->SetMinimum(0);
    //hdummy_stob[iplab]->SetMaximum(105);
    tl[2][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[3][iplab]->DrawLatex(0.35,0.75,Form("(%3.2f<|u|<%3.2f)",tmin[iplab],tmax[iplab]));
    tl[3][iplab]->DrawLatex(0.25,0.69,"Signal/Background");
    //if (iplab==0) leg[2]->Draw();
  }
  //stob_step->Print("figs/2014.12.16/stob_vs_steps.pdf");
  stob_step->Print("stob_vs_steps.pdf");
}

Double_t background2(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t background3(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

const double sqrt_2pi = TMath::Sqrt(2*TMath::Pi());
double binw = 0.0;
Double_t gaussianPeak(Double_t *x, Double_t *par) {
  double exp = (x[0]-par[1])/par[2];
  return binw*par[0] * TMath::Exp(-exp*exp/2.)/par[2]/sqrt_2pi;
}


Double_t fitFunctionPol2(Double_t *x, Double_t *par) {
  //return background(x,par) + langaufun(x,&par[3]);
  return background2(x,par) + gaussianPeak(x,&par[3]);
}

Double_t fitFunctionPol3(Double_t *x, Double_t *par) {
  //return background(x,par) + langaufun(x,&par[3]);
  return background3(x,par) + gaussianPeak(x,&par[4]);
}

void fig_ana2() {

  // Prepare to fit...
  const int ntu = 2;
  const int ntbin = 12;
  const char* toru[2] = {"t","u"};

  bool pol3 = false;

  double mmin = 2.0, mmax = 5.0;
  //double mmin = 2.9, mmax = 3.3; //-> 3sigma
  //double mmin = 2.8, mmax = 3.3;
  //double mmin = 2.0, mmax = 4.0;
  int immin = 0, immax = 0;

  TGraph *tg_sig_eff[nplab], *tg_bg_eff[nplab];
  TGraph *tg_sig_yield[nplab], *tg_bg_yield[nplab];
  TMultiGraph *tmg_yield[nplab], *tmg_eff[nplab];;

  TCanvas *tc_mep_tbins[ntu][nplab];

  TFile *fsig[nplab], *fbg[nplab];
  TH1F* hmsg_tu[ntu][nplab][ntbin], *hmbg_tu[ntu][nplab][ntbin], *hmfg_tu[ntu][nplab][ntbin];
  //TH1F* hmsg_tu_nofit[ntu][nplab][ntbin], *hmbg_tu_nofit[ntu][nplab][ntbin], *hmfg_tu_nofit[ntu][nplab][ntbin]; // DIRTY
  TF1* fmsg_tu[ntu][nplab][ntbin], *fmbg_tu[ntu][nplab][ntbin], *fmfg_tu[ntu][nplab][ntbin];

  TLatex *tl[10][nplab];
  for (int ii = 0; ii < 10; ++ii) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tl[ii][iplab] = new TLatex();
      tl[ii][iplab]->SetNDC(true);
      tl[ii][iplab]->SetTextColor(ii==0?2:1);
      tl[ii][iplab]->SetTextSize((ii==0?2.5:(iplab==0?1.5:2.0))*tl[ii][iplab]->GetTextSize());
    }
  }

  double t[ntbin]= {0.0}, t_er[ntbin] = {0.};
  double t_cnt[ntbin]= {0.0};
  for (int itbin = 0; itbin < ntbin; ++itbin) {
    t[itbin] = -0.45+0.1*itbin;
    t_cnt[itbin] = t[itbin]+0.01;
  }

  int mar[2] = {20,21};
  int mar_cnt[2] = {24,25};
  int col[3] = {1,2,4};
  double tvalid[ntu][nplab][ntbin]={{{0.0}}};
  double yield[ntu][nplab][ntbin]={{{0.0}}};
  double yield_er[ntu][nplab][ntbin]={{{0.0}}};
  double yield_cnt[ntu][nplab][ntbin]={{{0.0}}};
  double yield_cnt_er[ntu][nplab][ntbin]={{{0.0}}};
  int nptok[ntu][nplab] = {{0}};
  double chi2ndf[ntu][nplab][ntbin]={{{0.0}}};
  double prob[ntu][nplab][ntbin]={{{0.0}}};

  TGraphErrors *tg_yield[ntu][nplab];
  TGraphErrors *tg_yield_cnt[ntu][nplab];
  TGraph *tg_chi2ndf[ntu][nplab];
  TGraph *tg_prob[ntu][nplab];

  TMultiGraph *tmg_yield[ntu];
  TMultiGraph *tmg_yield_pbp[ntu][nplab];
  TMultiGraph *tmg_yield_cnt[ntu];
  TMultiGraph *tmg_yield_cnt_pbp[ntu][nplab];
  TMultiGraph *tmg_chi2ndf[ntu];
  TMultiGraph *tmg_prob[ntu];
  for (int itu = 0; itu < ntu; ++itu) {
    tmg_yield[itu] = new TMultiGraph(Form("tmg_yield_%s",toru[itu]),Form("tmg_yield_%s",toru[itu]));
    tmg_yield_cnt[itu] = new TMultiGraph(Form("tmg_yield_cnt_%s",toru[itu]),Form("tmg_yield_cnt_%s",toru[itu]));
    tmg_chi2ndf[itu] = new TMultiGraph(Form("tmg_chi2ndf_%s",toru[itu]),Form("tmg_chi2ndf_%s",toru[itu]));
    tmg_prob[itu] = new TMultiGraph(Form("tmg_prob_%s",toru[itu]),Form("tmg_prob_%s",toru[itu]));
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tmg_yield_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_pbp%s_t%d",toru[itu],iplab),Form("tmg_yield_pbp%s_p%d",toru[itu],iplab));
      tmg_yield_cnt_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_cnt_pbp%s_p%d",toru[itu],iplab),Form("tmg_yield_cnt_pbp%s_p%d",toru[itu],plab));
    }
  }

  TLine *tline = new TLine();
  tline->SetLineColor(4);
  tline->SetLineWidth(2);
  for (int iplab=0; iplab<nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("test/anav2_jpsi_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"), plab[iplab]));
    fbg[iplab] = TFile::Open(Form("test/anav2_pip_pim_%s_plab%3.1f.root",(ibrem==0?"raw":"brem"), plab[iplab]));
    for (int itu = 0; itu < ntu; ++itu) {
      tc_mep_tbins[itu][iplab] = new TCanvas(Form("fitted_mep_%sbins_p%d",toru[itu],iplab),Form("fitted_mep_%sbins_p%d",toru[itu],iplab));
      tc_mep_tbins[itu][iplab]->Divide(3,4);
      for (int itbin=0; itbin<ntbin; ++itbin) {

	hmfg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_fg_p%d_%s%db",iplab,toru[itu],itbin));
	hmsg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_sig_p%d_%s%d",iplab,toru[itu],itbin));
	hmbg_tu[itu][iplab][itbin] = (TH1F*) fbg[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_bg_p%d_%s%d",iplab,toru[itu],itbin));
	hmfg_tu[itu][iplab][itbin]->Add(hmbg_tu[itu][iplab][itbin]);
	hmfg_tu[itu][iplab][itbin]->SetTitle(";M_{inv}[GeV/c^{2}]");

	// Do counting before rebin, for more precise control
	immin = hmfg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmin);
	immax = hmfg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmax);
	double integral = hmsg_tu[itu][iplab][itbin]->Integral(immin, immax);
	yield_cnt[itu][iplab][itbin] = integral;
	yield_cnt_er[itu][iplab][itbin] = TMath::Sqrt(integral);

	set_style_ana(hmfg_tu[itu][iplab][itbin], 1, 4, true);
	set_style_ana(hmsg_tu[itu][iplab][itbin], 2, 4, true);
	set_style_ana(hmbg_tu[itu][iplab][itbin], 4, 4, true);

	hmfg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.0);
	hmsg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.0);
	hmbg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.0);

	binw = hmfg_tu[itu][iplab][itbin]->GetBinWidth(3);
	tc_mep_tbins[itu][iplab]->cd(itbin+1);
	// Define and setup fit function
	//fmfg_tu[itu][iplab][itbin] = new TF1(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),fitFunction,1.5,3.8,7);
	if (pol3) {
	  fmfg_tu[itu][iplab][itbin] = new TF1(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),fitFunctionPol3,2.3,3.8,7);
	  fmfg_tu[itu][iplab][itbin]->SetParameters(1,1,1,1,100,3.1,0.1);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(5,3.0,3.2);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(6,0.05,0.3);
	} else {
	  fmfg_tu[itu][iplab][itbin] = new TF1(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),fitFunctionPol2,2.3,3.8,6);
	  fmfg_tu[itu][iplab][itbin]->SetParameters(1,1,1,100,3.1,0.1);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(4,3.0,3.2);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(5,0.05,0.1);
	}
	fmfg_tu[itu][iplab][itbin]->SetNpx(500);
	fmfg_tu[itu][iplab][itbin]->SetLineWidth(2);
	fmfg_tu[itu][iplab][itbin]->SetLineColor(kCyan);


	if (integral>30000005) {
	  hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),"0Q");
	  fmfg_tu[itu][iplab][itbin]->SetParameter(pol3?4:3, fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3));
	  hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin), "Q+R", "ep");
	  cout << "par4 = " << fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3) << endl;
	  yield[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3);
	  yield_er[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParError(pol3?4:3);
	  tvalid[itu][iplab][nptok[itu][iplab]] = t[itbin];
	  nptok[itu][iplab]++;
	} else {
	  hmfg_tu[itu][iplab][itbin]->Draw();
	}

	//hmsg_tu[itu][iplab][itbin]->SetLineStyle(9);
	//hmbg_tu[itu][iplab][itbin]->SetLineStyle(9);
	//hmsg_tu[itu][iplab][itbin]->Draw("hist,same");
	//hmbg_tu[itu][iplab][itbin]->Draw("hist,same");

	chi2ndf[itu][iplab][itbin] = fmfg_tu[itu][iplab][itbin]->GetChisquare(); //fmfg_tu[itu][iplab][itbin]->GetNDF();
	prob[itu][iplab][itbin] = fmfg_tu[itu][iplab][itbin]->GetProb();
	tl[0][iplab]->DrawLatex(0.15,0.7,Form("%s",hmsg_tu[itu][iplab][itbin]->GetTitle()));

	//tline->DrawLine(mmin,0,mmin,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());
	//tline->DrawLine(mmax,0,mmax,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());
      }

      tc_mep_tbins[itu][iplab]->Print(Form("figs/2014.12.16/%s.pdf",tc_mep_tbins[itu][iplab]->GetName()));

      tg_yield[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],yield[itu][iplab],t_er,yield_er[itu][iplab]);
      tg_yield[itu][iplab]->SetMarkerStyle(mar[itu]);
      tg_yield[itu][iplab]->SetMarkerSize(1);
      tg_yield[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield[itu][iplab]->SetLineWidth(2);
      tmg_yield[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_pbp[itu][iplab]->Add(tg_yield[itu][iplab],"p");

      tg_yield_cnt[itu][iplab] = new TGraphErrors(ntbin-1,t_cnt,yield_cnt[itu][iplab],t_er,yield_cnt_er[itu][iplab]);
      tg_yield_cnt[itu][iplab]->SetMarkerStyle(mar_cnt[itu]);
      tg_yield_cnt[itu][iplab]->SetMarkerSize(1);
      tg_yield_cnt[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield_cnt[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield_cnt[itu][iplab]->SetLineWidth(2);
      tmg_yield_cnt[itu]->Add(tg_yield_cnt[itu][iplab],"p");
      tmg_yield_cnt[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_cnt_pbp[itu][iplab]->Add(tg_yield_cnt[itu][iplab],"p");
      tmg_yield_cnt_pbp[itu][iplab]->Add(tg_yield[itu][iplab],"p");

      tg_chi2ndf[itu][iplab] = new TGraph(ntbin-1,t,chi2ndf[itu][iplab]);
      tg_chi2ndf[itu][iplab]->SetMarkerStyle(mar[itu]);
      tg_chi2ndf[itu][iplab]->SetMarkerSize(1);
      tg_chi2ndf[itu][iplab]->SetMarkerColor(col[iplab]);
      tmg_chi2ndf[itu]->Add(tg_chi2ndf[itu][iplab],"p");

      tg_prob[itu][iplab] = new TGraph(ntbin-1,t,prob[itu][iplab]);
      tg_prob[itu][iplab]->SetMarkerStyle(mar[itu]);
      tg_prob[itu][iplab]->SetMarkerSize(1);
      tg_prob[itu][iplab]->SetMarkerColor(col[iplab]);
      tmg_prob[itu]->Add(tg_prob[itu][iplab],"p");
    }
  }

  TCanvas *tc_yield_pbp[ntu];
  TCanvas *tc_yield_cnt_pbp[ntu];
  TCanvas *tc_chi2ndf[ntu];
  TCanvas *tc_prob[ntu];

  TPad *pad_yield[ntu][nplab];
  TH1F *hdummy_yield[ntu][nplab];

  TPad *pad_yield_cnt[ntu][nplab];
  TH1F *hdummy_yield_cnt[ntu][nplab];

  double max_yield[ntu] = {0.0};
  for (int itu = 0; itu < ntu; ++itu) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double gr_max = 0;
      for (int itbin = 0; itbin < ntbin; ++itbin) {
	double top = yield_cnt[itu][iplab][itbin]+yield_cnt_er[itu][iplab][itbin];
	if (top > gr_max) { gr_max = top; }
      }
      if(gr_max  > max_yield[itu]) max_yield[itu] = gr_max;
    }
  }

  for (int itu = 0; itu < 1; ++itu) {

    tc_yield_pbp[itu] = new TCanvas(Form("fitted_yield_pbp%s",toru[itu]),Form("fitted_yield_pbp%s",toru[itu]));
    tc_yield_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double xl = (iplab==0?0:0.1)+iplab*0.3;
      double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
      pad_yield[itu][iplab] = new TPad(Form("pad_yield_%s_%d",toru[itu],iplab),Form("pad_yield_%s_%d",toru[itu],iplab),xl,0.0,xh,1.0);
      tc_yield_pbp[itu]->cd(0);
      pad_yield[itu][iplab]->Draw();
      double epsilon=1e-9;
      if (iplab==0) {pad_yield[itu][iplab]->SetRightMargin(epsilon); pad_yield[itu][iplab]->SetLeftMargin(0.2);}
      if (iplab==1) {pad_yield[itu][iplab]->SetLeftMargin(epsilon);  pad_yield[itu][iplab]->SetRightMargin(epsilon); }
      if (iplab==2) {pad_yield[itu][iplab]->SetLeftMargin(epsilon);  pad_yield[itu][iplab]->SetRightMargin(0.1);}
      pad_yield[itu][iplab]->SetTicks(0,1);
      pad_yield[itu][iplab]->cd();
      tmg_yield_pbp[itu][iplab]->Draw("a");
      hdummy_yield[itu][iplab] = (TH1F*)tmg_yield_pbp[itu][iplab]->GetHistogram();
      hdummy_yield[itu][iplab]->SetLabelSize(0.065,"Y");
      hdummy_yield[itu][iplab]->SetLabelSize(0.065,"X");
      hdummy_yield[itu][iplab]->SetLabelOffset(0.005,"Y");
      hdummy_yield[itu][iplab]->SetTitleSize(0.06,"X");
      hdummy_yield[itu][iplab]->SetTitleSize(0.06,"Y");
      hdummy_yield[itu][iplab]->SetTitleOffset(1.5,"Y");
      hdummy_yield[itu][iplab]->SetTitle(";t[GeV^{2}];yield");
      hdummy_yield[itu][iplab]->SetMinimum(0);
      hdummy_yield[itu][iplab]->SetMaximum(max_yield[itu]*1.1);
      tmg_yield_pbp[itu][iplab]->SetMinimum(0.0);
      tl[1][iplab]->DrawLatex((iplab==0?0.2:0.15),0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    }
    tc_yield_pbp[itu]->Print(Form("figs/2014.12.16/%s.pdf",tc_yield_pbp[itu]->GetName()));

    tc_yield_cnt_pbp[itu] = new TCanvas(Form("fitted_yield_cnt_pbp%s",toru[itu]),Form("fitted_yield_cnt_pbp%s",toru[itu]));
    tc_yield_cnt_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double xl = (iplab==0?0:0.1)+iplab*0.3;
      double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
      pad_yield_cnt[itu][iplab] = new TPad(Form("pad_yield_%s_%d",toru[itu],iplab),Form("pad_yield_%s_%d",toru[itu],iplab),xl,0.0,xh,1.0);
      tc_yield_cnt_pbp[itu]->cd(0);
      pad_yield_cnt[itu][iplab]->Draw();
      double epsilon=1e-9;
      if (iplab==0) {pad_yield_cnt[itu][iplab]->SetRightMargin(epsilon); pad_yield_cnt[itu][iplab]->SetLeftMargin(0.2);}
      if (iplab==1) {pad_yield_cnt[itu][iplab]->SetLeftMargin(epsilon);  pad_yield_cnt[itu][iplab]->SetRightMargin(epsilon); }
      if (iplab==2) {pad_yield_cnt[itu][iplab]->SetLeftMargin(epsilon);  pad_yield_cnt[itu][iplab]->SetRightMargin(0.1);}
      pad_yield_cnt[itu][iplab]->SetTicks(0,1);
      pad_yield_cnt[itu][iplab]->cd();
      tmg_yield_cnt_pbp[itu][iplab]->Draw("a");
      hdummy_yield_cnt[itu][iplab] = (TH1F*)tmg_yield_cnt_pbp[itu][iplab]->GetHistogram();

      hdummy_yield_cnt[itu][iplab]->SetLabelSize(0.065,"Y");
      hdummy_yield_cnt[itu][iplab]->SetLabelSize(0.065,"X");
      hdummy_yield_cnt[itu][iplab]->SetLabelOffset(0.005,"Y");
      hdummy_yield_cnt[itu][iplab]->SetTitleSize(0.06,"X");
      hdummy_yield_cnt[itu][iplab]->SetTitleSize(0.06,"Y");
      hdummy_yield_cnt[itu][iplab]->SetTitleOffset(1.5,"Y");
      hdummy_yield_cnt[itu][iplab]->SetTitle(";t[GeV^{2}];yield");
      hdummy_yield_cnt[itu][iplab]->SetMinimum(0);
      hdummy_yield_cnt[itu][iplab]->SetMaximum(max_yield[itu]*1.1);

      tmg_yield_cnt_pbp[itu][iplab]->SetMinimum(0.0);
      tl[1][iplab]->DrawLatex(iplab==0?0.2:0.15,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    }
    tc_yield_cnt_pbp[itu]->Print(Form("figs/2014.12.16/%s.pdf",tc_yield_cnt_pbp[itu]->GetName()));
    //tc_yield_cnt_pbp[itu] = new TCanvas(Form("fitted_yield_cnt_pbp%s",toru[itu]),Form("fitted_yield_cnt_pbp%s",toru[itu]));
    //tc_yield_cnt_pbp[itu]->Divide(3,1);
    //for (int iplab = 0; iplab < nplab; ++iplab) {
    //  tc_yield_cnt_pbp[itu]->cd(1+iplab);
    //  tmg_yield_cnt_pbp[itu][iplab]->Draw("a");
    //  tmg_yield_cnt_pbp[itu][iplab]->SetMinimum(0.0);
    //}

    continue;
    tc_chi2ndf[itu] = new TCanvas(Form("fit_chi2ndf_%s",toru[itu]),Form("fit_chi2ndf_%s",toru[itu]));
    tc_chi2ndf[itu]->cd();
    tmg_chi2ndf[itu]->Draw("a");
    tc_prob[itu] = new TCanvas(Form("fit_prob_%s",toru[itu]),Form("fit_prob_%s",toru[itu]));
    tc_prob[itu]->cd();
    tmg_prob[itu]->Draw("a");
    tc_yield[itu]->Print(Form("figs/2014.12.16/%s.pdf",tc_yield[itu]->GetName()));
  }

}


static const double mpi = 0.135;
static const double mj = 3.096;
static const double mp = 0.938;
double s(double _plab) {
  return 2.*mp*mp + 2.*mp*TMath::Hypot(mp,_plab);
}
double pcm(double _plab) {
  return _plab*mp/TMath::Sqrt(s(_plab));
}
double q(double _plab) {
  return TMath::Sqrt( pow(s(_plab),2) + pow(mj*mj-mpi*mpi,2) - (2*s(_plab)*(mj*mj+mpi*mpi)))/2./TMath::Sqrt(s(_plab));
}
double costh(double t, double _plab) {
  return ( t + pow(pcm(_plab),2) + pow(q(_plab),2) - pow(TMath::Hypot(q(_plab),mpi) - TMath::Hypot(pcm(_plab),mp), 2) )/2./q(_plab)/pcm(_plab);
}

void fig_tvsthcm_sig() {

  double umin[nplab] = {0.0}, umax[nplab] = {0.0};
  for (int iplab=0; iplab<nplab; ++iplab) {
    umin[iplab] = limits[iplab][0];
    umax[iplab] = limits[iplab][0]-(tmin[iplab]-tmax[iplab]);
  }

  int col[3] = {1,2,4};

  TLatex *tt[6][nplab];
  for (int iplab = 0; iplab < nplab; ++iplab) {
    for (int i=0; i<6; ++i) {
      tt[i][iplab] = new TLatex();
      tt[i][iplab]->SetNDC(true);
      tt[i][iplab]->SetTextColor(col[iplab]);
      if (tt>3) tt[i][iplab]->SetTextSize(0.7*tt[i][iplab]->GetTextSize());
      //tt[i][iplab]->SetTextColor(col[iplab]);
    }
  }

  TFile *f[3];
  TH2F* tvscost[3];
  TCanvas *tc[4];
  for (int icanv = 0; icanv < 4; ++icanv) tc[icanv] = new TCanvas(Form("tc%d",icanv),Form("tc%d",icanv),1400,1000);
  for (int iplab=nplab-1; iplab>=0; --iplab) {
    f[iplab]= TFile::Open(Form("hists/sig%d_nofilt_noeff.root",iplab));
    tvscost[iplab] = (TH2F*)f[iplab]->Get("inv_p_t_pbar_pi0_cm_cost_pi0")->Clone(Form("tvscost%d",3-iplab));
    tvscost[iplab]->SetMarkerColor(col[iplab]);
    for (int icanv = 0; icanv < 4; ++icanv) {
      tc[icanv]->cd();
      tvscost[iplab]->SetTitle(";t[GeV^{2}];cos(#theta_{#pi^{0}}^{CM})");
      tvscost[iplab]->Draw(iplab==nplab-1?"":"same");
    }
  }

  tc[3]->cd();
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tt[3][iplab]->DrawLatex(0.15,0.75-iplab*0.1,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  }
  tc[3]->Update();
  tc[3]->Print((Form("figs/2014.12.16/validity_iplab3.png")));
  gSystem->Exec(Form("convert figs/2014.12.16/validity_iplab3.png figs/2014.12.16/validity_iplab3.pdf"));

  TLine *tl_t[nplab][2][2]; // [plab][vert-horiz][lower-upper]
  TLine *tl_u[nplab][2][2]; // [plab][vert-horiz][lower-upper]
  TBox *tb_t[nplab][2];
  TBox *tb_u[nplab][2];
  for (int iplab=0; iplab<nplab; ++iplab) {
    cout << costh(tmin[iplab], plab[iplab]) << endl;
    tc[iplab]->cd();

    tl_t[iplab][0][0] = new TLine(tmin[iplab], -1.1, tmin[iplab], costh(tmin[iplab], plab[iplab]));
    tl_t[iplab][0][1] = new TLine(tmax[iplab], -1.1, tmax[iplab], costh(tmax[iplab], plab[iplab]));
    tl_t[iplab][1][0] = new TLine(-14.0, costh(tmin[iplab], plab[iplab]), tmin[iplab], costh(tmin[iplab], plab[iplab]));
    tl_t[iplab][1][1] = new TLine(-14.0, costh(tmax[iplab], plab[iplab]), tmax[iplab], costh(tmax[iplab], plab[iplab]));

    tl_u[iplab][0][0] = new TLine(umin[iplab], -1.1, umin[iplab], costh(umin[iplab], plab[iplab]));
    tl_u[iplab][0][1] = new TLine(umax[iplab], -1.1, umax[iplab], costh(umax[iplab], plab[iplab]));
    tl_u[iplab][1][0] = new TLine(-14.0, costh(umin[iplab], plab[iplab]), umin[iplab], costh(umin[iplab], plab[iplab]));
    tl_u[iplab][1][1] = new TLine(-14.0, costh(umax[iplab], plab[iplab]), umax[iplab], costh(umax[iplab], plab[iplab]));

    for (int hv=0; hv<2; ++hv) {
      for (int lu=0; lu<2; ++lu) {
	tl_t[iplab][hv][lu]->SetLineColor(col[iplab]);
	tl_t[iplab][hv][lu]->SetLineWidth(2);
	tl_t[iplab][hv][lu]->Draw();
	tl_u[iplab][hv][lu]->SetLineColor(col[iplab]);
	tl_u[iplab][hv][lu]->SetLineWidth(2);
	tl_u[iplab][hv][lu]->Draw();
      }
    }

    tb_t[iplab][0] = new TBox(-14.0, costh(tmin[iplab], plab[iplab]), tmax[iplab], costh(tmax[iplab], plab[iplab]));
    tb_t[iplab][1] = new TBox(tmin[iplab], -1.1, tmax[iplab], costh(tmax[iplab], plab[iplab]));
    for (int i=0; i<2; ++i) {
      tb_t[iplab][i]->SetFillStyle(3001);
      tb_t[iplab][i]->SetFillColor(col[iplab]);
    }
    tb_t[iplab][0]->Draw();
    tb_t[iplab][1]->Draw();

    tb_u[iplab][0] = new TBox(-14.0, costh(umin[iplab], plab[iplab]), umax[iplab], costh(umax[iplab], plab[iplab]));
    tb_u[iplab][1] = new TBox(umin[iplab], -1.1, umax[iplab], costh(umax[iplab], plab[iplab]));
    for (int i=0; i<2; ++i) {
      tb_u[iplab][i]->SetFillStyle(3003);
      tb_u[iplab][i]->SetFillColor(col[iplab]);
    }
    tb_u[iplab][0]->Draw();
    tb_u[iplab][1]->Draw();

    tt[0][iplab]->DrawLatex(0.15,0.45,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

    tt[1][iplab]->DrawLatex(iplab==0?0.4:(iplab==1?0.35:0.25),iplab==0?0.3:(iplab==1?0.22:0.20),"|u|<<Q^{2}");
    tt[4][iplab]->DrawLatex(iplab==0?0.36:(iplab==1?0.31:0.21),iplab==0?0.25:(iplab==1?0.17:0.15),Form("%3.2f<|u|<%3.2f",tmin[iplab],tmax[iplab]));

    tt[2][iplab]->DrawLatex(iplab==0?0.35:(iplab==1?0.4:0.45),iplab==0?0.71:(iplab==1?0.82:0.78),"|t|<<Q^{2}");
    tt[5][iplab]->DrawLatex(iplab==0?0.31:(iplab==1?0.36:0.41),iplab==0?0.67:(iplab==1?0.78:0.74),Form("%3.2f<|t|<%3.2f",tmin[iplab],tmax[iplab]));

    tc[iplab]->Update();

    tc[iplab]->Print(Form("figs/2014.12.16/validity_iplab%d.png",iplab));
    gSystem->Exec(Form("convert figs/2014.12.16/validity_iplab%d.png figs/2014.12.16/validity_iplab%d.pdf", iplab, iplab));
  }

}

void fig_bg_xsect() {
  ifstream inf;
  inf.open("figs/bg_xsect_data.txt");
  const int npt = 27;
  double ecm[npt], ecm_er[npt], _plab[npt], _plab_er[npt], xsect[npt], xsect_er[npt];
  int ctr = 0, ipt;
  while(true) {
    inf >> ipt;
    if (!inf.good() || ipt!= ctr || ipt > 26) break;
    inf >> ecm[ipt] >> _plab[ipt] >>  xsect[ipt] >> xsect_er[ipt];
    ecm_er[ipt] = _plab_er[ipt] = 0.0;
    cout << "ipt= " << ipt << " " << ecm[ipt] << " " << _plab[ipt] << " " << xsect[ipt] << " " << xsect_er[ipt] << endl;
    ++ctr;
  }
  TGraphErrors *tge_ecm= new TGraphErrors(npt, ecm,xsect,ecm_er,xsect_er);
  tge_ecm->SetTitle(";E_{cm}[GeV];#sigma(#bar{p}+p#rightarrow#pi^{0})[mb]");
  set_style(tge_ecm, 2);

  TGraphErrors *tge_plab= new TGraphErrors(npt, _plab,xsect,_plab_er,xsect_er);
  tge_plab->SetTitle(";p_{LAB}^{#bar{p}}[GeV];#sigma(#bar{p}+p#rightarrow#pi^{0})[mb]");
  set_style(tge_plab, 2, 20, 1, 2);

  TCanvas *tc_ecm = new TCanvas("tc_ecm","tc_ecm");
  tc_ecm->cd();
  tge_ecm->Draw("ap");
  tge_ecm->SetMaximum(10);
  gPad->SetLogy();

  TCanvas *tc_plab = new TCanvas("tc_plab","tc_plab");
  tc_plab->cd();
  tge_plab->Draw("ap");
  tge_plab->SetMaximum(10);
  tge_plab->GetHistogram()->GetXaxis()->SetRange(0, 15);
  gPad->SetLogy();

  TLine *tl[3];
  for (int iplab=0; iplab<nplab; ++iplab) {
    tl[iplab] = new TLine();
    tl[iplab]->SetLineWidth(4);
    tl[iplab]->SetLineColor(1);
    cout << "plab = " << plab[iplab] << endl;
    tl[iplab]->DrawLine(plab[iplab], tge_plab->GetMinimum(), plab[iplab], tge_plab->GetMaximum());
  }
}

void figs(int _page) {
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  switch(_page) {
  case 0: fig_tvsthcm_sig(); break;
  case 1: fig_bg_xsect(); break;
  case 2: fig_ana(); break;
  case 3: fig_pi0cut(); break;
  case 4: fig_npair(); break;
  case 5: fig_gen_dists(); break;
  case 6: fig_kin_cuts(); break;
  case 7: fig_num_evt(); break;
  case 8: fig_ana2(); break;
  default: return;
  }
}

void set_style(TH1* h, int col, int rebin, bool sumw2) {
  if (rebin>0)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void set_style_ana(TH1* h, int col, int rebin, bool sumw2) {
  if (rebin>1)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.08);
  h->GetXaxis()->SetLabelSize(0.08);
  h->GetXaxis()->SetLabelOffset(0.);
  h->GetXaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.08);
  h->GetYaxis()->SetLabelSize(0.08);
  h->GetYaxis()->SetLabelOffset(0.02);
  h->GetYaxis()->SetTitleOffset(0.8);

  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void set_style(TGraph *g, int col, int sty, int s, int w) {
  g->SetMarkerStyle(sty);
  g->SetMarkerSize(s);
  g->SetMarkerColor(col);
  g->SetLineColor(col);
  g->SetLineWidth(w);
}

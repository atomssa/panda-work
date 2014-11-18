#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1F.h"

void set_style(TH1* h, int col=4, double scale=1.) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Scale(scale);
  if (col>0) {
    h->SetLineWidth(3);
    h->SetLineColor(col);
  }
}

void set_style2(TH1* h, int col=4, double sty=1) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  if (col>0) {
    h->SetLineWidth(3);
    h->SetLineColor(col);
    h->SetLineStyle(sty);
  }
}

TLegend* set_tl(double xmin, double ymin, double xmax, double ymax, int ts=0.05) {
  TLegend *r = new TLegend(xmin,ymin,xmax,ymax);
  r->SetTextSize(0.05);
  r->SetFillStyle(0);
  r->SetBorderSize(0);
  return r;
}

TLine *set_line(double xmin, double ymin, double xmax, double ymax) {
  TLine *line;
  line = new TLine(xmin,ymin,xmax,ymax);
  line->SetLineColor(8);
  line->SetLineWidth(2);
  return line;
}

static TFile *f_brem[2];
static TCanvas *tc;
static TLegend *tl[3];
static double em[8] = {0.005,0.01,0.015,0.02,0.025,0.05,0.075,0.1};
static bool bgsig = false;
static bool savefig = true;
static int page;

void fig0() {

  // photon, pi+, pi-, e+, e-
  double ph_max[] = {16.e-3, 0.01, 0.01, 0.01, 0.01 };
  double th_max[] = {6.e-3, 0.005, 0.005, 0.005, 0.005 };
  int col[3] = {1,2,4};
  TEllipse *ellipse[5][3];
  for (int isp=0; isp<5; ++isp) {
    for (int ie=0; ie<3; ++ie) {
      ellipse[isp][ie] = new TEllipse(0,0,ph_max[isp]*(1.0+ie), th_max[isp]*(1.0+ie));
      ellipse[isp][ie]->SetFillStyle(0);
      ellipse[isp][ie]->SetLineColor(col[ie]);
    }
  }

  TH2F *phth[5];
  phth[0] = (TH2F*) f_brem[0]->Get("resid/h_resid_phth0");
  phth[1] = (TH2F*) f_brem[0]->Get("resid/h_resid_pip_phth0");
  phth[2] = (TH2F*) f_brem[0]->Get("resid/h_resid_pim_phth0");
  phth[3] = (TH2F*) f_brem[1]->Get("resid/h_resid_ep_phth0");
  phth[4] = (TH2F*) f_brem[1]->Get("resid/h_resid_em_phth0");
  TH1F *mom[4];
  mom[0] = (TH1F*) f_brem[0]->Get("resid/h_resid_pip_mom0");
  mom[1] = (TH1F*) f_brem[0]->Get("resid/h_resid_pim_mom0");
  mom[2] = (TH1F*) f_brem[1]->Get("resid/h_resid_ep_mom0");
  mom[3] = (TH1F*) f_brem[1]->Get("resid/h_resid_em_mom0");

  tc = new TCanvas("tc_resid","tc_resid");
  tc->Divide(3,3);
  for (int isp = 0; isp < 5; isp++) {
    tc->cd(isp+1);
    phth[isp]->Draw("colz");
    for (int ie = 0; ie < 3; ++ie) {
      ellipse[isp][ie]->Draw();
    }
    if (isp>0) {
      tc->cd(isp+5);
      mom[isp-1]->Draw();
    }
  }

}

void fig1() {

  tc = new TCanvas("tc_pi0_sel","tc_pi0_sel",1400,900);
  tc->Divide(3,2);

  TCanvas *tc_loc[2];
  tc_loc[0] = new TCanvas("tc_pi0_sel_bg","tc_pi0_sel_bg",1400,900);
  tc_loc[1] = new TCanvas("tc_pi0_sel_sig","tc_pi0_sel_sig",1400,900);
  tc_loc[0]->Divide(3,1);
  tc_loc[1]->Divide(3,1);
  TLegend *tl_loc[2][3];
  for (int ibs=0; ibs<2; ++ibs){
    for (int ip=0; ip<3; ++ip ) {
      tl_loc[ibs][ip] = new TLegend(0.15,0.75,1.1,0.9);
      tl_loc[ibs][ip]->SetFillStyle(0);
      tl_loc[ibs][ip]->SetBorderSize(0);
    }
  }

  TText *ttt = new TText();
  ttt->SetTextSize(ttt->GetTextSize()*1.5);

  TH1F *hgg[2][3][3];
  const char *nt[3] = {"all_rec", "truepi0_rec", "truepi0_mc"};
  const char *np[3] = {"m_gg", "oa_gg", "e_g"};
  const char *lt[3] = {"All (Reco)", " Pair Matched (Reco)", "Pair Matched (MC)"};
  const int c[3] = {1, 2, 4};
  const double s[2][3][3] = {  {{1, 1, 1}, {1, 20, 20}, {1, 1, 1} }  , {{1, 1, 1}, {1, 5, 5}, {1, 1, 1} } };

  for (int ibs=0; ibs<2; ++ibs) {
    for (int p=0; p<3; ++p) {
      for (int t=0; t<3; ++t) {
	const char *name = Form("%s/h_%s_%s",(p<2?"gg":"single"),np[p], nt[t]);
	cout << "name = " << name << endl;

	hgg[ibs][p][t] = (TH1F*) f_brem[ibs]->Get(name);
	set_style(hgg[ibs][p][t],c[t], s[ibs][p][t]);

	tc_loc[ibs]->cd(1+p);
	hgg[ibs][p][t]->Draw(t==0?"":"same");
	if (p!=1) gPad->SetLogy();

	tc->cd((ibs*3)+1+p);
	hgg[ibs][p][t]->Draw(t==0?"":"same");
	//if (p!=1)
	gPad->SetLogy();

	tl_loc[ibs][p]->AddEntry(hgg[ibs][p][t],Form("%s%s", (s[ibs][p][t]!=1?Form("(#times %2.0f)",s[ibs][p][t]):""), lt[t]),"l");
      }
    }
    tc_loc[ibs]->cd(1);
    (set_line(0.1,0,0.1,hgg[ibs][1][0]->GetMaximum()))->Draw();
    (set_line(0.16,0,0.16,hgg[ibs][1][0]->GetMaximum()))->Draw();
    tc_loc[ibs]->cd(2);
    (set_line(0.2,0,0.2,hgg[ibs][1][0]->GetMaximum()))->Draw();
    tc_loc[ibs]->cd(3);
    (set_line(em[5],0,em[5],hgg[ibs][1][0]->GetMaximum()))->Draw();

    tc->cd(ibs*3+1);
    (set_line(0.1,0,0.1,hgg[ibs][1][0]->GetMaximum()))->Draw();
    (set_line(0.16,0,0.16,hgg[ibs][1][0]->GetMaximum()))->Draw();
    tc->cd(ibs*3+2);
    (set_line(0.2,0,0.2,hgg[ibs][1][0]->GetMaximum()))->Draw();
    tc->cd(ibs*3+3);
    (set_line(em[5],0,em[5],hgg[ibs][1][0]->GetMaximum()))->Draw();

    //for (int i=0; i<8; ++i) (set_line(em[i],0,em[i],hgg[ibs][1][0]->GetMaximum()))->Draw();
    for (int p=0; p<3; ++p) {
      tc_loc[ibs]->cd(1+p);
      tl_loc[ibs][p]->Draw();
      tc->cd((ibs*3)+1+p);
      tl_loc[ibs][p]->Draw();
      tc->cd((ibs*3)+1+p);
      ttt->DrawTextNDC(0.6,0.6,(ibs==0?"Background":"Signal"));
    }
  }


}

void fig2() {
  for (int i=0; i<2; ++i) tl[i] = set_tl(0.15, 0.75, 1.1, 0.9);
  const char *nt[4] = {"all_rec", "truejpsi_rec", "truejpsi_mc", "pm_ana"};
  const char *lt[4] = {"All (Reco)", " Pair Matched (Reco)", "Pair Matched (MC)", "(2) With mass cut"};
  tc = new TCanvas("tc_jpsi","tc_cut_jpsi",1400,1000);
  tc->Divide(2,1);
  const int cc[4] = {1, 2, 4, 6};
  TH1F *hgg[2][4];
  for (int ff=0; ff<2; ++ff) {
    for (int i=0; i<3; ++i) {
      if (ff==1 && (i==2)) continue;
      const char *name = Form("epem/h_m_epem_%s",nt[i]);
      cout << "name = " << name << endl;
      hgg[ff][i] = (TH1F*) f_brem[ff]->Get(name);
      tc->cd(ff+1);
      set_style(hgg[ff][i], cc[i]);
      hgg[ff][i]->Draw((i==0?"","same"));

      //if (ff!=0||(ff==0&&i==0))
	tl[ff]->AddEntry(hgg[ff][i],Form("%s",lt[i]), "l");

    }
  }
  bgsig = true;
  for (int p=0; p<2; ++p) {
    tc->cd(1+p);
    tl[p]->Draw();
  }
}


void fig3() {

  for (int i=0; i<2; ++i) tl[i] = set_tl(0.15, 0.65, 1.0, 0.9, 0.05);

  TLegend *tl_loc[2][2];
  TCanvas *tc_loc[2][2];
  for (int itype=0; itype<2; ++itype) {
    for (int imcut=0; imcut<2; ++imcut) {

      tc_loc[itype][imcut] = new TCanvas(Form("tc_cut_dep_type%d_mcut%d",itype,imcut),Form("tc_cut_dep_type%d_mcut%d",itype,imcut),1400,1000);
      tc_loc[itype][imcut]->Divide(2,1);

      tl_loc[itype][imcut] = new TLegend(0.45,0.4,0.9,0.8);
      tl_loc[itype][imcut]->SetFillStyle(0);
      tl_loc[itype][imcut]->SetBorderSize(0);
    }
  }

  const int cc[8] = {1, 2, 4, 6, 1, 2, 4, 6};
  const int cs[8] = {1, 1, 1, 1, 8, 8, 8, 8};
  TH1F *hgg[2][8];

  int ff = 0;
  int cut_base = 0;

  double eff[2][2][8] = {{{0.}}};
  for (int itype=0; itype<2; ++itype) {
    for (int imcut=0; imcut<2; ++imcut) {
      for (int icut=0; icut<8; ++icut) {

	const char *name0 = Form("pi0_ana/h_m_gg_ana_%d", (imcut==0?0:8)+icut);
	hgg[0][icut] = (TH1F*) f_brem[itype]->Get(name0);
	set_style2(hgg[0][icut], cc[icut], cs[icut]);
	tc_loc[itype][imcut]->cd(1);
	hgg[0][icut]->Draw((icut==0?"":"same"));

	const char *name1 = Form("pi0_ana/h_m_gg_pm_ana_%d", (imcut==0?0:8)+icut);
	hgg[1][icut] = (TH1F*) f_brem[itype]->Get(name1);
	set_style2(hgg[1][icut], cc[icut], cs[icut]);
	tc_loc[itype][imcut]->cd(2);
	hgg[1][icut]->Draw((icut==0?"":"same"));

	const char *txt = Form("%s, OA%s + E> %3.0f MeV", (itype==0?"BG":"SIG"), (imcut==0?"":" + M"), em[icut]*1e3);
	tl_loc[itype][imcut]->AddEntry(hgg[0][icut], txt, "l");

	eff[itype][imcut][icut] = (100.*hgg[1][icut]->GetEntries())/(itype==0?1e5:1e4);
	//cout << "eff(E> " << em[icut] << ") = " << eff[itype][imcut][icut] << endl;

      }

      tc_loc[itype][imcut]->cd(1);
      tl_loc[itype][imcut]->Draw();

    }
  }

  //for (int p=0; p<1; ++p) {
  //  tc->cd(1+p);
  //  tl[p]->Draw();
  //  tc->cd(1+p);
  //  tl[p]->Draw();
  //}

  savefig = false;

  TMultiGraph *tmg =  new TMultiGraph("eff","Efficiency vs. Min. #gamma_{E} cut; Min(#gamma_{E}); #varepsilon_{#pi^{0}}");
  TGraphErrors *tge[2][2];
  for (int itype=0; itype<2; ++itype) {
    for (int imcut=0; imcut<2; ++imcut) {
      tge[itype][imcut] = new TGraphErrors(8, em, eff[itype][imcut], 0, 0);
      tge[itype][imcut]->SetMarkerStyle(imcut==0?20:24);
      tge[itype][imcut]->SetMarkerSize(2);
      tge[itype][imcut]->SetMarkerColor(itype==0?2:4);
      tmg->Add(tge[itype][imcut],"pl");
    }
  }

  TCanvas *tc_eff = new TCanvas();
  tc_eff->cd();
  tmg->Draw("a");
  tmg->SetMinimum(0);
  tmg->SetMaximum(100);
}

void fig4() {
  tc = new TCanvas("tc_sig_bg0","tc_sig_bg0",1400,1000);
  TH1F *hepm[3];
  for (int i=0; i<2; ++i) {
    hepm[i] = (TH1F*) f_brem[i]->Get("epem/h_m_epem_ana");
  }
  hepm[2] = (TH1F*) (TH2F*) f_brem[0]->Get("epem/h_m_epem_ana")->Clone("h_m_epem_ana_bg");
  hepm[2]->Add((TH1F*) f_brem[1]->Get("h_m_epem_ana"));
  set_style(hepm[0],4);
  set_style(hepm[1],2);
  set_style(hepm[2],1);
  tc->cd();
  hepm[1]->Draw();
  hepm[0]->Draw("same");
  hepm[2]->Draw("same");
}

void fig6() {

  TCanvas *tc_loc[2];
  tc_loc[0] = new TCanvas("tc_dth_s","tc_dth_s_ana",1400,1000);
  tc_loc[1] = new TCanvas("tc_dth_s_1d","tc_dth_s_ana_1d",1400,1000);

  TH2F *hdths[2];
  TH2F *hdths_pm[2];
  TH1F *hdth[2];
  TH1F *hdth_pm[2];
  TH1F *hs[2];
  TH1F *hs_pm[2];

  for (int ff=0; ff<2; ++ff) {
    hdths[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_mass_gg_epair_ana");
    hdths_pm[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_mass_gg_epair_pm_ana");

    hdth[ff] = (TH1F*) hdths[ff]->ProjectionY()->Clone(ff==0?"hdth_bg":"hdth_sig");
    hdth_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionY()->Clone(ff==0?"hdth_pm_bg":"hdth_pm_sig");
    hs[ff] = (TH1F*) hdths[ff]->ProjectionX()->Clone(ff==0?"hs_bg":"hs_sig");
    hs_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionX()->Clone(ff==0?"hs_pm_bg":"hs_pm_sig");
  }

  const double xmin = hdths[0]->GetXaxis()->GetBinLowEdge(1);
  const double xmax = hdths[0]->GetXaxis()->GetBinLowEdge(hdths[0]->GetXaxis()->GetNbins());
  const double ymin = hdths[0]->GetYaxis()->GetBinLowEdge(1);
  const double ymax = hdths[0]->GetYaxis()->GetBinLowEdge(hdths[0]->GetXaxis()->GetNbins());
  const double dthcut_min = 3.0;
  const double dthcut_max = 3.3;
  const double scut_min = 3.4;
  const double scut_max = 3.6;

  TBox *tb_s = new TBox(scut_min,dthcut_min,scut_max,dthcut_max);
  tb_s->SetLineWidth(4);
  tb_s->SetLineColor(2);
  tb_s->SetFillStyle(0);

  // 2D
  TText *tt = new TText();
  tt->SetTextSize(0.08);
  TLine *tlv = new TLine(xmin,TMath::Pi(),xmax,TMath::Pi());
  tlv->SetLineWidth(2);
  TLine *tlh = new TLine(3.5,ymin,3.5,ymax);
  tlh->SetLineWidth(2);
  tc_loc[0]->Divide(2,2);
  tc_loc[0]->cd(1);
  hdths[1]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Signal (All)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();
  tc_loc[0]->cd(2);
  hdths_pm[1]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Signal (Pair Matched)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();
  tc_loc[0]->cd(3);
  hdths[0]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Background (All)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();
  tc_loc[0]->cd(4);
  hdths_pm[0]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Background (Pair Matched)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();

  // 1D
  TLine *tl_s_min = new TLine(scut_min,0,scut_min,80);
  TLine *tl_s_max = new TLine(scut_max,0,scut_max,80);
  TLine *tl_th_min = new TLine(dthcut_min,0,dthcut_min,160);
  TLine *tl_th_max = new TLine(dthcut_max,0,dthcut_max,160);
  tl_s_min->SetLineWidth(2);
  tl_s_max->SetLineWidth(2);
  tl_th_min->SetLineWidth(2);
  tl_th_max->SetLineWidth(2);

  tc_loc[1]->Divide(2,2);
  set_style(hdth[1],1,4);
  set_style(hdth_pm[1],2,4);
  set_style(hdth[0],4);
  set_style(hdth_pm[0],6);
  set_style(hs[1],1,1.8);
  set_style(hs_pm[1],2,1.8);
  set_style(hs[0],4);
  set_style(hs_pm[0],6);

  tl[0] = set_tl(0.15, 0.75, 0.45, 0.9);
  tl[0]->AddEntry(hdth[1], "Signal, all #pi^{0}-J/\#psi pairs");
  tl[0]->AddEntry(hdth[0], "BG, all #pi^{0}#pi^{+}#pi^{-} pairs");

  tl[1] = set_tl(0.15, 0.75, 0.45, 0.9);
  tl[1]->AddEntry(hdth_pm[1], "Signal, MC matching #pi^{0}-J/\#psi pairs");
  tl[1]->AddEntry(hdth_pm[0], "BG, MC matching  #pi^{0}#pi^{+}#pi^{-} pairs");

  tc_loc[1]->cd(1);
  hdth[0]->Draw();
  hdth[1]->Draw("same");
  tl[0]->Draw();
  tl_th_min->Draw();
  tl_th_max->Draw();

  tc_loc[1]->cd(2);
  hdth_pm[1]->Draw();
  hdth_pm[0]->Draw("same");
  tl[1]->Draw();
  tl_th_min->Draw();
  tl_th_max->Draw();

  tc_loc[1]->cd(3);
  hs[0]->Draw();
  hs[0]->SetMinimum(0);
  hs[1]->Draw("same");
  tl_s_min->Draw();
  tl_s_max->Draw();

  tc_loc[1]->cd(4);
  hs_pm[1]->Draw();
  hs_pm[0]->Draw("same");
  tl_s_min->Draw();
  tl_s_max->Draw();

  savefig = false;
  tc_loc[0]->Print(Form("figs/2014.11.06/fig_p%d_sub1.pdf",page));
  tc_loc[1]->Print(Form("figs/2014.11.06/fig_p%d_sub2.pdf",page));

}

void fig5() {

  TCanvas *tc_loc[2];
  tc_loc[0] = new TCanvas("tc_dth_dphi","tc_dth_dphi_ana",1400,1000);
  tc_loc[1] = new TCanvas("tc_dth_dphi_1d","tc_dth_dphi_ana_1d",1400,1000);

  TH2F *hdthdph[2];
  TH2F *hdthdph_pm[2];
  TH1F *hdth[2];
  TH1F *hdth_pm[2];
  TH1F *hdph[2];
  TH1F *hdph_pm[2];
  for (int ff=0; ff<2; ++ff) {
    hdthdph[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_dph_gg_epair_ana");
    hdthdph_pm[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_dph_gg_epair_pm_ana");
    //hdthdph[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_dph_gg_epair_ana");
    //hdthdph_pm[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_dph_gg_epair_pm_ana");
    hdth[ff] = (TH1F*) hdthdph[ff]->ProjectionY()->Clone(ff==0?"hdth_bg":"hdth_sig");
    hdth_pm[ff] = (TH1F*) hdthdph_pm[ff]->ProjectionY()->Clone(ff==0?"hdth_pm_bg":"hdth_pm_sig");
    hdph[ff] = (TH1F*) hdthdph[ff]->ProjectionX()->Clone(ff==0?"hdph_bg":"hdph_sig");
    hdph_pm[ff] = (TH1F*) hdthdph_pm[ff]->ProjectionX()->Clone(ff==0?"hdph_pm_bg":"hdph_pm_sig");
  }
  const double xmin = hdthdph[0]->GetXaxis()->GetBinLowEdge(1);
  const double xmax = hdthdph[0]->GetXaxis()->GetBinLowEdge(hdthdph[0]->GetXaxis()->GetNbins());
  const double ymin = hdthdph[0]->GetYaxis()->GetBinLowEdge(1);
  const double ymax = hdthdph[0]->GetYaxis()->GetBinLowEdge(hdthdph[0]->GetXaxis()->GetNbins());

  // 2D
  TText *tt = new TText();
  tt->SetTextSize(0.08);
  TLine *tlv = new TLine(xmin,TMath::Pi(),xmax,TMath::Pi());
  tlv->SetLineWidth(2);
  TLine *tlh = new TLine(TMath::Pi(),ymin,TMath::Pi(),ymax);
  tlh->SetLineWidth(2);
  tc_loc[0]->Divide(2,2);
  tc_loc[0]->cd(1);
  hdthdph[1]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Signal (All)");
  tlv->Draw();
  tlh->Draw();
  tc_loc[0]->cd(2);
  hdthdph_pm[1]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Signal (Pair Matched)");
  tlv->Draw();
  tlh->Draw();
  tc_loc[0]->cd(3);
  hdthdph[0]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Background (All)");
  tlv->Draw();
  tlh->Draw();
  tc_loc[0]->cd(4);
  hdthdph_pm[0]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Background (Pair Matched)");
  tlv->Draw();
  tlh->Draw();

  // 1D
  tc_loc[1]->Divide(2,2);
  set_style(hdth[1],1,4);
  set_style(hdth_pm[1],2,4);
  set_style(hdth[0],4);
  set_style(hdth_pm[0],6);
  set_style(hdph[1],1,5);
  set_style(hdph_pm[1],2,5);
  set_style(hdph[0],4);
  set_style(hdph_pm[0],6);

  tl[0] = set_tl(0.15, 0.75, 0.45, 0.9);
  tl[0]->AddEntry(hdth[1], "Signal, all #pi^{0}-J/\#psi pairs");
  tl[0]->AddEntry(hdth[0], "BG, all #pi^{0}#pi^{+}#pi^{-} pairs");

  tl[1] = set_tl(0.15, 0.75, 0.45, 0.9);
  tl[1]->AddEntry(hdth_pm[1], "Signal, MC matching #pi^{0}-J/\#psi pairs");
  tl[1]->AddEntry(hdth_pm[0], "BG, MC matching  #pi^{0}#pi^{+}#pi^{-} pairs");

  tc_loc[1]->cd(1);
  hdth[0]->Draw();
  hdth[1]->Draw("same");
  tl[0]->Draw();

  tc_loc[1]->cd(2);
  hdth_pm[1]->Draw();
  hdth_pm[0]->Draw("same");
  tl[1]->Draw();

  tc_loc[1]->cd(3);
  hdph[0]->Draw();
  hdph[0]->SetMinimum(0);
  hdph[1]->Draw("same");

  tc_loc[1]->cd(4);
  hdph_pm[1]->Draw();
  hdph_pm[0]->Draw("same");

  savefig = false;
  tc_loc[0]->Print(Form("figs/2014.11.06/fig_p%d_sub1.pdf",page));
  tc_loc[1]->Print(Form("figs/2014.11.06/fig_p%d_sub2.pdf",page));
}




void fig8() {

  TCanvas *tc_loc[2];
  tc_loc[0] = new TCanvas("tc_dth_s","tc_dth_s_ana",1400,1000);
  tc_loc[1] = new TCanvas("tc_dth_s_1d","tc_dth_s_ana_1d",1400,1000);

  TH2F *hdths[2];
  TH2F *hdths_pm[2];
  TH1F *hdth[2];
  TH1F *hdth_pm[2];
  TH1F *hs[2];
  TH1F *hs_pm[2];

  for (int ff=0; ff<2; ++ff) {
    // THIS IS TOO HACKY, PLEASE FIX!!
    hdths[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_mass_gg_epair_pm_ana");
    hdths_pm[ff] = (TH2F*) f_brem[ff]->Get("full_sys/h_dth_vs_mass_gg_epair_mconst");

    hdth[ff] = (TH1F*) hdths[ff]->ProjectionY()->Clone(ff==0?"hdth_bg":"hdth_sig");
    hdth_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionY()->Clone(ff==0?"hdth_pm_bg":"hdth_pm_sig");
    hs[ff] = (TH1F*) hdths[ff]->ProjectionX()->Clone(ff==0?"hs_bg":"hs_sig");
    hs_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionX()->Clone(ff==0?"hs_pm_bg":"hs_pm_sig");
  }

  const double xmin = hdths[0]->GetXaxis()->GetBinLowEdge(1);
  const double xmax = hdths[0]->GetXaxis()->GetBinLowEdge(hdths[0]->GetXaxis()->GetNbins());
  const double ymin = hdths[0]->GetYaxis()->GetBinLowEdge(1);
  const double ymax = hdths[0]->GetYaxis()->GetBinLowEdge(hdths[0]->GetXaxis()->GetNbins());
  const double dthcut_min = 3.0;
  const double dthcut_max = 3.3;
  const double scut_min = 3.4;
  const double scut_max = 3.6;

  TBox *tb_s = new TBox(scut_min,dthcut_min,scut_max,dthcut_max);
  tb_s->SetLineWidth(4);
  tb_s->SetLineColor(2);
  tb_s->SetFillStyle(0);

  // 2D
  TText *tt = new TText();
  tt->SetTextSize(0.08);
  TLine *tlv = new TLine(xmin,TMath::Pi(),xmax,TMath::Pi());
  tlv->SetLineWidth(2);
  TLine *tlh = new TLine(3.5,ymin,3.5,ymax);
  tlh->SetLineWidth(2);
  tc_loc[0]->Divide(2,2);
  tc_loc[0]->cd(1);
  hdths[1]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Signal (Pair Matched)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();
  tc_loc[0]->cd(2);
  hdths_pm[1]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Signal (Mass Const.)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();
  tc_loc[0]->cd(3);
  hdths[0]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Background (Pair Matched)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();
  tc_loc[0]->cd(4);
  hdths_pm[0]->Draw("colz");
  tt->DrawTextNDC(0.15,0.15,"Background (Mass Const.)");
  tlv->Draw();
  tlh->Draw();
  tb_s->Draw();

  // 1D
  TLine *tl_s_min = new TLine(scut_min,0,scut_min,80);
  TLine *tl_s_max = new TLine(scut_max,0,scut_max,80);
  TLine *tl_th_min = new TLine(dthcut_min,0,dthcut_min,160);
  TLine *tl_th_max = new TLine(dthcut_max,0,dthcut_max,160);
  tl_s_min->SetLineWidth(2);
  tl_s_max->SetLineWidth(2);
  tl_th_min->SetLineWidth(2);
  tl_th_max->SetLineWidth(2);

  tc_loc[1]->Divide(2,2);
  set_style(hdth[1],1,1);
  set_style(hdth_pm[1],2,1);
  set_style(hdth[0],4);
  set_style(hdth_pm[0],6);
  set_style(hs[1],1,1);
  set_style(hs_pm[1],2,1);
  set_style(hs[0],4);
  set_style(hs_pm[0],6);

  tl[0] = set_tl(0.15, 0.75, 0.45, 0.9);
  tl[0]->AddEntry(hdth[1], "Signal, all #pi^{0}-J/\#psi pairs");
  tl[0]->AddEntry(hdth[0], "BG, all #pi^{0}#pi^{+}#pi^{-} pairs");

  tl[1] = set_tl(0.15, 0.75, 0.45, 0.9);
  tl[1]->AddEntry(hdth_pm[1], "Signal, MC matching #pi^{0}-J/\#psi pairs");
  tl[1]->AddEntry(hdth_pm[0], "BG, MC matching  #pi^{0}#pi^{+}#pi^{-} pairs");

  tc_loc[1]->cd(1);
  hdth[0]->Draw();
  hdth[1]->Draw("same");
  tl[0]->Draw();
  tl_th_min->Draw();
  tl_th_max->Draw();

  tc_loc[1]->cd(2);
  hdth_pm[1]->Draw();
  hdth_pm[0]->Draw("same");
  tl[1]->Draw();
  tl_th_min->Draw();
  tl_th_max->Draw();

  tc_loc[1]->cd(3);
  hs[0]->Draw();
  hs[0]->SetMinimum(0);
  hs[1]->Draw("same");
  tl_s_min->Draw();
  tl_s_max->Draw();

  tc_loc[1]->cd(4);
  hs_pm[1]->Draw();
  hs_pm[0]->Draw("same");
  tl_s_min->Draw();
  tl_s_max->Draw();

  savefig = false;
  tc_loc[0]->Print(Form("figs/2014.11.06/fig_p%d_sub1.pdf",page));
  tc_loc[1]->Print(Form("figs/2014.11.06/fig_p%d_sub2.pdf",page));

}

void fig7() {

  //TCanvas *tc_loc[5];
  //for (int i=0; i<5; ++i) tc_loc[i] = new TCanvas(Form("tc_eff_%d",i),Form("tc_pi0_sel_%d",i),1400,900);

  tc = new TCanvas("tc_pi0_sel","tc_pi0_sel",1400,900);
  tc->Divide(2,1);
  int ie[] = {0, 1, 2, 3};
  int col[6] = {1,2,4,7,9,6};
  TH1F* h_eff_thpi0[2][5];
  cout << "loop1" << endl;
  for (int t=0; t<2; ++t) {
    for (int e=0; e<4; ++e) {
      cout << "loop" << t << " " << e << " hbrem[t]= " << f_brem[t] << endl;
      //cout << "hist= " << ((TH1F*) f_brem[t]->Get(Form("h_eff_thpi0_%d",e))) << endl;
      h_eff_thpi0[t][e] = ((TH1F*) f_brem[t]->Get(Form("eff/h_eff_thpi0_%d",ie[e])));//->Clone(Form("h_eff_thpi0_%d_%s",e,(t=0?"bg":"sig")));
      cout << "loop" << t << " " << ie[e] << " done" << endl;
      tc->cd(t+1);
      //h_eff_thpi0[t][e]->GetXaxis()->SetRangeUser(0,0.6);
      h_eff_thpi0[t][e]->SetLineColor(col[e]);
      h_eff_thpi0[t][e]->SetLineWidth(2);
      h_eff_thpi0[t][e]->Draw(e==0?"":"same");
      if (e==0)
	h_eff_thpi0[t][e]->SetTitle("#theta_{#pi^{0}}^{lab}");
      if (t==0)
	tl[0]->AddEntry(h_eff_thpi0[t][e], (e==0?"MC Matched":h_eff_thpi0[t][e]->GetTitle()), "l");
    }
  }
  tc->cd(1);
  tl[0]->Draw();
}

void fig9() {
  tc = new TCanvas("tc_money_plot","tc_money_plot",1400,1000);
  TH1F *hepm[3];
  for (int i=0; i<2; ++i) {
    hepm[i] = (TH1F*) f_brem[i]->Get("man_kin_fit/h_m_epem_btb_ana");
  }
  hepm[2] = (TH1F*) (TH2F*) f_brem[0]->Get("man_kin_fit/h_m_epem_btb_ana")->Clone("h_m_epem_ana_bg");
  hepm[2]->Add((TH1F*) f_brem[1]->Get("man_kin_fit/h_m_epem_btb_ana"));
  set_style(hepm[0],4);
  set_style(hepm[1],2);
  set_style(hepm[2],1);
  tc->cd();
  hepm[2]->Draw();
  hepm[1]->Draw("same");
  hepm[0]->Draw("same");
  break;
}

void fig10() {
  tc = new TCanvas("tc_money_plot","tc_money_plot",1400,1000);
  TH1F *hm[4];
  hm[0] = (TH1F*) f_brem[0]->Get("epem/h_m_epem_pm_ana");
  hm[1] = (TH1F*) f_brem[0]->Get("epem/h_m_epem_mconst");
  hm[2] = (TH1F*) f_brem[1]->Get("epem/h_m_epem_pm_ana");
  hm[3] = (TH1F*) f_brem[1]->Get("epem/h_m_epem_mconst");

  tc->Divide(2,2);
  for (int i=0; i<4; ++i) {
    tc->cd(1+i);
    cout << hm[i] << endl;
    hm[i]->Draw();
  }

}

void figs(int _page) {

  page = _page;
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  // 0 - bg, 1 - Sig
  f_brem[0] = TFile::Open("test/ana_bg_brem.root");
  f_brem[1] = TFile::Open("test/ana_jpsi_brem.root");

  tl[0] = new TLegend(0.15,0.6,0.9,0.8);
  tl[0]->SetFillStyle(0);
  tl[0]->SetBorderSize(0);
  tl[1] = new TLegend(0.15,0.6,0.9,0.8);
  tl[1]->SetFillStyle(0);
  tl[1]->SetBorderSize(0);
  tl[2] = new TLegend(0.15,0.6,0.9,0.8);
  tl[2]->SetFillStyle(0);
  tl[2]->SetBorderSize(0);

  switch(page) {
  case 0: fig1(); break;
  case 1: fig1(); break;
  case 2: fig2(); break;
  case 3: fig3(); break;
  case 4: fig4(); break;
  case 5: fig5(); break;
  case 6: fig6(); break;
  case 7: fig7(); break;
  case 8: fig8(); break;
  case 9: fig9(); break;
  case 10: fig10(); break;
  case 11: fig11(); break;
  case 12: fig12(); break;
  case 13: fig13(); break;
  default: return;
  }

  if (bgsig) {
    TText *tt = new TText();
    tc->cd(1);
    tt->SetTextSize(0.08);
    tt->DrawTextNDC(0.44,0.9,"Background");
    tc->cd(2);
    tt->DrawTextNDC(0.65,0.9,"Signal");
  }

  if (savefig)
    tc->Print(Form("figs/2014.11.06/fig_p%d.pdf",page));

}

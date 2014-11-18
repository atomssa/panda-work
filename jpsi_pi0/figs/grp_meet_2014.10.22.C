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


void fig1() {
  tc = new TCanvas("tc_pi0_sel","tc_pi0_sel",1400,900);
  tc->Divide(3,1);
  TH1F *hgg[3][3];
  for (int i=0; i<3; ++i) tl[i] = set_tl(0.15, 0.75, 1.1, 0.9);
  const char *nt[3] = {"all_rec", "truepi0_rec", "truepi0_mc"};
  const char *np[3] = {"m_gg", "oa_gg", "e_g"};
  const char *lt[3] = {"All (Reco)", " Pair Matched (Reco)", "Pair Matched (MC)"};
  const int c[3] = {1, 2, 4};
  const double s[3][3] = {{1, 1, 1}, {1, 10, 10}, {1, 1, 1} };
  int f = 1;
  for (int p=0; p<3; ++p) {
    for (int t=0; t<3; ++t) {
      const char *name = Form("h_%s_%s",np[p], nt[t]);
      cout << "name = " << name << endl;
      hgg[p][t] = (TH1F*) f_brem[f]->Get(name);
      tc->cd(1+p);
      if (p==0&&t==2) continue;
      hgg[p][t]->Draw(t==0?"":"same");
      set_style(hgg[p][t],c[t], s[p][t]);
      gPad->SetLogy();
      tl[p]->AddEntry(hgg[p][t],Form("%s%s", lt[t], s[p][t]!=1?Form("(#times %2.0f)",s[p][t]):""),"l");
    }
  }
  tc->cd(1);
  (set_line(0.1,0,0.1,hgg[1][0]->GetMaximum()))->Draw();
  (set_line(0.16,0,0.16,hgg[1][0]->GetMaximum()))->Draw();
  tc->cd(2);
  (set_line(0.2,0,0.2,hgg[1][0]->GetMaximum()))->Draw();
  tc->cd(3);
  for (int i=0; i<8; ++i) (set_line(em[i],0,em[i],hgg[1][0]->GetMaximum()))->Draw();
  for (int p=0; p<3; ++p) {
    tc->cd(1+p);
    tl[p]->Draw();
  }
}

void fig2() {
  for (int i=0; i<2; ++i) tl[i] = set_tl(0.15, 0.65, 1.0, 0.9, 0.05);
  tc = new TCanvas("tc_cut_dep","tc_cut_dep",1400,1000);
  tc->Divide(2,1);
  const int cc[8] = {1, 2, 4, 6, 1, 2, 4, 6};
  const int cs[8] = {1, 1, 1, 1, 8, 8, 8, 8};
  TH1F *hgg[2][8];
  int ff = 1;
  int cut_base = 8;
  double eff[8] = {0.};
  for (int i=0; i<2; ++i) {
    tc->cd(i+1);
    for (int c=0; c<8; ++c) {
      const char *name = Form("h_m_gg_%s_%d",(i==0?"ana":"pm_ana"), cut_base+c);
      //cout << "name = " << name << endl;
      hgg[i][c] = (TH1F*) f_brem[ff]->Get(name);
      set_style2(hgg[i][c], cc[c], cs[c]);
      //if (i==0) { gPad->SetLogy();}
      hgg[i][c]->Draw((c==0?"","same"));
      //cout << "hgg= " << hgg[i][c] << endl;;
      const char *txt = Form("OA%s + E> %4.3f GeV", (cut_base==0?"":" + M"), em[c]);
      tl[i]->AddEntry(hgg[i][c], txt, "l");
      if (i==1) {
	eff[c] = (100.*hgg[i][c]->GetEntries())/1e4;
	cout << "eff(E> " << em[c] << ") = " << eff[c] << endl;
      }
    }
  }

  for (int p=0; p<1; ++p) {
    tc->cd(1+p);
    tl[p]->Draw();
    tc->cd(1+p);
    tl[p]->Draw();
  }

  savefig = false;
  //tc->Print(Form("figs/2014.10.22/fig_p%d_file%d_cut%d.pdf",page,ff,cut_base));

}

void fig3() {
  for (int i=0; i<2; ++i) tl[i] = set_tl(0.15, 0.75, 1.1, 0.9);
  const char *nt[4] = {"all_rec", "truejpsi_rec", "truejpsi_mc", "pm_ana"};
  const char *lt[4] = {"All (Reco)", " Pair Matched (Reco)", "Pair Matched (MC)", "(2) With mass cut"};
  tc = new TCanvas("tc_jpsi","tc_cut_jpsi",1400,1000);
  tc->Divide(2,1);
  const int cc[4] = {1, 2, 4, 6};
  double eff[2] = {0.};
  TH1F *hgg[2][4];
  for (int ff=0; ff<2; ++ff) {
    for (int i=0; i<4; ++i) {
      const char *name = Form("h_m_epem_%s",nt[i]);
      cout << "name = " << name << endl;
      hgg[ff][i] = (TH1F*) f_brem[ff]->Get(name);
      tc->cd(ff+1);
      set_style(hgg[ff][i], cc[i]);
      hgg[ff][i]->Draw((i==0?"","same"));
      if (ff!=0||(ff==0&&i==0)) tl[ff]->AddEntry(hgg[ff][i],Form("%s",lt[i]), "l");
    }
    if (ff==1) {
      eff[ff] = hgg[ff][3]->GetEntries()/hgg[ff][2]->GetEntries();
      cout << "eff_sig= " << eff[ff] << endl;
    }
  }
  bgsig = true;
  for (int p=0; p<2; ++p) {
    tc->cd(1+p);
    tl[p]->Draw();
  }
}

void fig4() {
  tc = new TCanvas("tc_sig_bg0","tc_sig_bg0",1400,1000);
  TH1F *hepm[3];
  for (int i=0; i<2; ++i) {
    hepm[i] = (TH1F*) f_brem[i]->Get("h_m_epem_ana");
  }
  hepm[2] = (TH1F*) f_brem[0]->Get("h_m_epem_ana")->Clone("h_m_epem_ana_bg");
  hepm[2]->Add((TH1F*) f_brem[1]->Get("h_m_epem_ana"));
  set_style(hepm[0],4);
  set_style(hepm[1],2);
  set_style(hepm[2],1);
  tc->cd();
  hepm[2]->Draw();
  hepm[1]->Draw("same");
  hepm[0]->Draw("same");
}

void fig5() {
  tc = new TCanvas("tc_dth_ana","tc_dth_ana",1400,1000);
  TH2F *hdths[2];
  TH2F *hdths_pm[2];
  TH1F *hdth[2];
  TH1F *hdth_pm[2];
  TH1F *hs[2];
  TH1F *hs_pm[2];
  for (int ff=0; ff<2; ++ff) {
    hdths[ff] = (TH2F*) f_brem[ff]->Get("h_dth_vs_mass_gg_epair_ana");
    hdths_pm[ff] = (TH2F*) f_brem[ff]->Get("h_dth_vs_mass_gg_epair_pm_ana");
    hdth[ff] = (TH1F*) hdths[ff]->ProjectionY()->Clone(ff==0?"hdth_bg":"hdth_sig");
    hdth_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionY()->Clone(ff==0?"hdth_pm_bg":"hdth_pm_sig");
    hs[ff] = (TH1F*) hdths[ff]->ProjectionX()->Clone(ff==0?"hs_bg":"hs_sig");
    hs_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionX()->Clone(ff==0?"hs_pm_bg":"hs_pm_sig");
  }
  int sub = 1;
  if (sub==0) {
    tc->Divide(3,1);
    tc->cd(1);
    hdths[1]->Draw("colz");
    cout << "sig = " << hdths[1]->GetEntries() << endl;;
    tc->cd(2);
    hdths_pm[1]->Draw("colz");
    cout << "sig_pm = " << hdths_pm[1]->GetEntries() << endl;;
    tc->cd(3);
    hdths[0]->Draw("colz");
    cout << "bg = " << hdths[0]->GetEntries() << endl;;
    TText *tt = new TText();
    tc->cd(1);
    tt->SetTextSize(0.08);
    tt->DrawTextNDC(0.15,0.15,"Signal (All)");
    tc->cd(2);
    tt->DrawTextNDC(0.15,0.15,"Signal (MC Match)");
    tc->cd(3);
    tt->DrawTextNDC(0.15,0.15,"Background (MC Match)");
  } else if (sub==1) {
    tc->Divide(2,1);
    tc->cd(1);
    set_style(hdth[1],1,2.6);
    set_style(hdth_pm[1],2,6.5);
    set_style(hdth[0],4);
    hdth[0]->Draw();
    hdth[1]->Draw("same");
    hdth_pm[1]->Draw("same");
    tl[0] = set_tl(0.15, 0.75, 1.1, 0.9);
    tl[0]->AddEntry(hdth[1], "Signal, all #pi^{0}-J\#psi pairs");
    tl[0]->AddEntry(hdth_pm[1], "Signal, MC matching #pi^{0}-J\#psi pairs");
    tl[0]->AddEntry(hdth[0], "BG, all #pi^{0}-J\#psi pairs");
    tl[0]->Draw();
    tc->cd(2);
    set_style(hs[1],1,2.6);
    set_style(hs_pm[1],2,6.5);
    set_style(hs[0],4);
    hs[0]->Draw();
    hs[0]->SetMinimum(0);
    hs[1]->Draw("same");
    hs_pm[1]->Draw("same");
  }
  savefig = false;
  tc->Print(Form("figs/2014.10.22/fig_p%d_sub%d.pdf",page,sub));
}

void fig6() {
  tc = new TCanvas("tc_dth_ana","tc_dth_ana",1400,1000);
  TH2F *hdths[2];
  TH2F *hdths_pm[2];
  TH1F *hdth[2];
  TH1F *hdth_pm[2];
  TH1F *hs[2];
  TH1F *hs_pm[2];
  for (int ff=0; ff<2; ++ff) {
    hdths[ff] = (TH2F*) f_brem[ff]->Get("h_dth_vs_mass_gg_epair_btb_ana");
    hdths_pm[ff] = (TH2F*) f_brem[ff]->Get("h_dth_vs_mass_gg_epair_btb_pm_ana");
    hdth[ff] = (TH1F*) hdths[ff]->ProjectionY()->Clone(ff==0?"hdth_bg":"hdth_sig");
    hdth_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionY()->Clone(ff==0?"hdth_pm_bg":"hdth_pm_sig");
    hs[ff] = (TH1F*) hdths[ff]->ProjectionX()->Clone(ff==0?"hs_bg":"hs_sig");
    hs_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionX()->Clone(ff==0?"hs_pm_bg":"hs_pm_sig");
  }
  int sub = 1;
  if (sub==0) {
    tc->Divide(3,1);
    tc->cd(1);
    hdths[1]->Draw("colz");
    cout << "sig = " << hdths[1]->GetEntries() << endl;;
    tc->cd(2);
    hdths_pm[1]->Draw("colz");
    cout << "sig_pm = " << hdths_pm[1]->GetEntries() << endl;;
    tc->cd(3);
    hdths[0]->Draw("colz");
    cout << "bg = " << hdths[0]->GetEntries() << endl;;
    TText *tt = new TText();
    tc->cd(1);
    tt->SetTextSize(0.08);
    tt->DrawTextNDC(0.15,0.15,"Signal (All)");
    tc->cd(2);
    tt->DrawTextNDC(0.15,0.15,"Signal (MC Match)");
    tc->cd(3);
    tt->DrawTextNDC(0.15,0.15,"Background (MC Match)");
  } else if (sub==1) {
    tc->Divide(2,1);
    tc->cd(1);
    set_style(hdth[1],1,2);
    set_style(hdth_pm[1],2,2);
    set_style(hdth[0],4);
    hdth[0]->Draw();
    hdth[1]->Draw("same");
    hdth_pm[1]->Draw("same");
    tl[0] = set_tl(0.15, 0.75, 1.1, 0.9);
    tl[0]->AddEntry(hdth[1], "Signal, all #pi^{0}-J\#psi pairs");
    tl[0]->AddEntry(hdth_pm[1], "Signal, MC matching #pi^{0}-J\#psi pairs");
    tl[0]->AddEntry(hdth[0], "BG, all #pi^{0}-J\#psi pairs");
    tl[0]->Draw();
    tc->cd(2);
    set_style(hs[1],1,1.7);
    set_style(hs_pm[1],2,1.7);
    set_style(hs[0],4);
    hs[0]->Draw();
    hs[0]->SetMinimum(0);
    hs[1]->Draw("same");
    hs_pm[1]->Draw("same");
  }

  savefig = false;
  tc->Print(Form("figs/2014.10.22/fig_p%d_sub%d.pdf",page,sub));
}

void fig7() {
  tc = new TCanvas("tc_dth_ana","tc_dth_ana",1400,1000);
  TH2F *hdths[2];
  TH2F *hdths_pm[2];
  TH1F *hdth[2];
  TH1F *hdth_pm[2];
  TH1F *hs[2];
  TH1F *hs_pm[2];
  for (int ff=0; ff<2; ++ff) {
    hdths[ff] = (TH2F*) f_brem[ff]->Get("h_dth_vs_mass_gg_epair_cts_ana");
    hdths_pm[ff] = (TH2F*) f_brem[ff]->Get("h_dth_vs_mass_gg_epair_cts_pm_ana");
    hdth[ff] = (TH1F*) hdths[ff]->ProjectionY()->Clone(ff==0?"hdth_bg":"hdth_sig");
    hdth_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionY()->Clone(ff==0?"hdth_pm_bg":"hdth_pm_sig");
    hs[ff] = (TH1F*) hdths[ff]->ProjectionX()->Clone(ff==0?"hs_bg":"hs_sig");
    hs_pm[ff] = (TH1F*) hdths_pm[ff]->ProjectionX()->Clone(ff==0?"hs_pm_bg":"hs_pm_sig");
  }
  int sub = 1;
  if (sub==0) {
    tc->Divide(3,1);
    tc->cd(1);
    hdths[1]->Draw("colz");
    cout << "sig = " << hdths[1]->GetEntries() << endl;;
    tc->cd(2);
    hdths_pm[1]->Draw("colz");
    cout << "sig_pm = " << hdths_pm[1]->GetEntries() << endl;;
    tc->cd(3);
    hdths[0]->Draw("colz");
    cout << "bg = " << hdths[0]->GetEntries() << endl;;
    TText *tt = new TText();
    tc->cd(1);
    tt->SetTextSize(0.08);
    tt->DrawTextNDC(0.15,0.15,"Signal (All)");
    tc->cd(2);
    tt->DrawTextNDC(0.15,0.15,"Signal (MC Match)");
    tc->cd(3);
    tt->DrawTextNDC(0.15,0.15,"Background (MC Match)");
  } else if (sub==1) {
    tc->Divide(2,1);
    tc->cd(1);
    set_style(hdth[1],1,2);
    set_style(hdth_pm[1],2,2);
    set_style(hdth[0],4);
    hdth[0]->Draw();
    hdth[1]->Draw("same");
    hdth_pm[1]->Draw("same");
    tl[0] = set_tl(0.15, 0.75, 1.1, 0.9);
    tl[0]->AddEntry(hdth[1], "Signal, all #pi^{0}-J\#psi pairs");
    tl[0]->AddEntry(hdth_pm[1], "Signal, MC matching #pi^{0}-J\#psi pairs");
    tl[0]->AddEntry(hdth[0], "BG, all #pi^{0}-J\#psi pairs");
    tl[0]->Draw();
    tc->cd(2);
    set_style(hs[1],1,1.7);
    set_style(hs_pm[1],2,1.7);
    set_style(hs[0],4);
    hs[0]->Draw();
    hs[0]->SetMinimum(0);
    hs[1]->Draw("same");
    hs_pm[1]->Draw("same");
  }
  savefig = false;
  tc->Print(Form("figs/2014.10.22/fig_p%d_sub%d.pdf",page,sub));
}

void fig8() {
  tc = new TCanvas("tc_money_plot","tc_money_plot",1400,1000);
  TH1F *hepm[3];
  for (int i=0; i<2; ++i) {
    hepm[i] = (TH1F*) f_brem[i]->Get("h_m_epem_btb_ana");
  }
  hepm[2] = (TH1F*) f_brem[0]->Get("h_m_epem_btb_ana")->Clone("h_m_epem_ana_bg");
  hepm[2]->Add((TH1F*) f_brem[1]->Get("h_m_epem_btb_ana"));
  set_style(hepm[0],4);
  set_style(hepm[1],2);
  set_style(hepm[2],1);
  tc->cd();
  hepm[2]->Draw();
  hepm[1]->Draw("same");
  hepm[0]->Draw("same");
  break;
}

void fig9() {
  tc = new TCanvas("tc_prob_4c_epem","tc_prob_4c_epem",1400,1000);
  tc->Divide(2,1);
  TH1F *h_prob_4c_cts[2];
  for (int i=0; i<2; ++i) {
    h_prob_4c_cts[i] = (TH1F*) f_brem[i]->Get("h_4c_prob_vs_m_cts_epem_4c");
    tc->cd(i+1);
    h_prob_4c_cts[i]->Draw("colz");
  }
}

void fig10() {
  tc = new TCanvas("tc_prob_4c_epem","tc_prob_4c_epem",1400,1000);
  tc->Divide(2,1);
  TH1F *h_prob_4c_cts[2];
  for (int i=0; i<2; ++i) {
    h_prob_4c_cts[i] = (TH1F*) f_brem[i]->Get("h_4c_prob_cts_epempi0");
    tc->cd(i+1);
    h_prob_4c_cts[i]->Draw("colz");
    gPad->SetLogy();
  }
}

void figs(int page) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  // 0 - bg, 1 - Sig
  f_brem[0] = TFile::Open("hists/no_mcut/ana_bg_brem.root");
  f_brem[1] = TFile::Open("hists/no_mcut/ana_jpsi_brem.root");

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
    tc->Print(Form("figs/2014.10.22/fig_p%d.pdf",page));

}

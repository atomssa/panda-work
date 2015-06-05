#include "TEfficiency.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TGAxis.h"

#include <iostream>

void set_style(TH2* h) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);

  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);

  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
}

void set_style(TH1* h, int col) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);

  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);

  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}


double avg_eff(TEfficiency *eff, TH1* wt) {
  int ndim = eff->GetDimension();
  const TH1* href = eff->GetTotalHistogram();
  int nbinx = href->GetXaxis()->GetNbins();
  if (wt->GetXaxis()->GetNbins()!=nbinx) {
    cout << "ERROR: eff hist and weight hist not same binning" << endl;
  }
  double num = 0;
  double den = 0;
  for (int ibinx=0; ibinx < nbinx; ++ibinx) {
    double _eff = eff->GetEfficiency(ibinx);
    double _wt = wt->GetBinContent(ibinx);
    num += _eff*_wt;
    den += _wt;
    //num += _eff;
    //den += 1.0;
  }
  return num/den;
}

double avg_eff2d(TEfficiency *eff, TH2* _wt, int rebin) {
  if (eff->GetDimension()!=2) {
    cout << "2d average not possible with 1d eff" << endl;
    return 0;
  }
  TH2F* wt = (TH2F*) _wt->Clone("_wt");
  int ndim = eff->GetDimension();
  const TH1* href = eff->GetTotalHistogram();
  int nbinx = href->GetXaxis()->GetNbins();
  int nbiny = href->GetYaxis()->GetNbins();

  if (wt->GetXaxis()->GetNbins()!=nbinx) {
    cout << "ERROR: eff hist and weight hist not same binning" << endl;
    cout << "nbinx= " << nbinx << " nbiny= " << nbiny << endl;
    cout << "nbinxWt= " << wt->GetXaxis()->GetNbins() << " nbiny= " << wt->GetYaxis()->GetNbins() << endl;
  }
  if (wt->GetYaxis()->GetNbins()!=nbiny) {
    cout << "ERROR: eff hist and weight hist not same binning" << endl;
    cout << "nbinx= " << nbinx << " nbiny= " << nbiny << endl;
    cout << "nbinxWt= " << wt->GetXaxis()->GetNbins() << " nbiny= " << wt->GetYaxis()->GetNbins() << endl;
  }

  if (rebin>=2) {
    cout << "rebinning wt" << endl;
    wt->RebinX(rebin);
    wt->RebinY(rebin);
    TH2F* tmpN = (TH2F*) eff->GetPassedHistogram()->Clone("tmpN");
    TH2F* tmpD = (TH2F*) eff->GetTotalHistogram()->Clone("tmpD");
    cout << "rebinning eff" << endl;
    tmpN->RebinX(rebin);
    tmpN->RebinY(rebin);
    tmpD->RebinX(rebin);
    tmpD->RebinY(rebin);
    eff = new TEfficiency(*tmpN,*tmpD);
  }

  nbinx = href->GetXaxis()->GetNbins();
  nbiny = href->GetYaxis()->GetNbins();
  double num = 0;
  double den = 0;
  for (int ibinx=0; ibinx < nbinx; ++ibinx) {
    for (int ibiny=0; ibiny < nbiny; ++ibiny) {
      double _eff = eff->GetEfficiency(eff->GetGlobalBin(ibinx,ibiny));
      double _wt = wt->GetBinContent(ibinx,ibiny);
      num += _eff*_wt;
      den += _wt;
      //num += _eff;
      //den += 1.0;
    }
  }
  return num/den;
}

void dist() {

  gStyle->SetOptStat(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  const int np=3;
  const int nt=2;
  const char *tt[2] = {"pi0pipm_dpm","pi0jpsi"};
  const char *rxn_tex[2] = {"#bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-}", "#bar{p}p#rightarrow #pi^{0}J/#psi(e^{+}e^{-})"};
  const char *plab_tex[3] = {"p_{#bar{p}}= 5.513 GeV/c", "p_{#bar{p}}= 8 GeV/c", "p_{#bar{p}}= 12 GeV/c"};

  TFile *f;
  TH2F *h_wt[nt][np];
  TCanvas *tc_wt = new TCanvas("tc_wt","tc_wt");
  tc_wt->Divide(3,2);
  TCanvas *tc_wt_sb = new TCanvas("tc_wt_sb","tc_wt_sb");
  tc_wt_sb->Divide(3,2);
  TLine *tlv = new TLine();
  TLine *tlh = new TLine();
  TLine *tlv2 = new TLine();
  TLine *tlh2 = new TLine();
  TBox *tb[3];
  int col[np] = {1, 2, 4};
  for (int ib=0; ib < 3; ++ib) {
    tb[ib] = new TBox();
    tb[ib]->SetLineColor(col[ib]);
    tb[ib]->SetLineWidth(3);
    tb[ib]->SetFillStyle(0);
  }
  TLatex *tlat = new TLatex();
  tlat->SetLineColor(2);
  tlat->SetNDC();
  tlat->SetTextSize(0.06);
  for (int ip=0; ip < np; ++ip) {
    for (int it=0; it < nt; ++it) {
      f = TFile::Open(Form("ana_%s_ip%d_raw.root",tt[it],ip));
      h_wt[it][ip] = (TH2F*) f->Get("epem/h_mom_the_epm_all_rec");
      const char *ttt = (it==0?"#pi^{+} and #pi^{-} (DPM)":"e^{+} and e^{-}(TDA)");
      h_wt[it][ip]->SetTitle(Form("#theta vs mom of %s",ttt));
      set_style(h_wt[it][ip]);
      double min = it==0?1.5:1.2;
      double max = it==0?1.5:1.2;
      tc_wt_sb->cd(np*it+ip+1);
      h_wt[it][ip]->Draw("colz");
      tlat->DrawLatex(0.4,0.7,rxn_tex[it]);
      tlat->DrawLatex(0.4,0.6,plab_tex[ip]);
      tc_wt->cd(np*it+ip+1);
      h_wt[it][ip]->Draw("colz");
      tb[0]->DrawBox(0,0,5,TMath::Pi()/2);
      tb[1]->DrawBox(0,0,max,TMath::Pi());
      tb[2]->DrawBox(5,0,10,TMath::Pi()/6);
    }
  }
  tc_wt->Print("physics_distributions.pdf");
  tc_wt_sb->Print("physics_distributions_sb.pdf");

  //TFile *fpim = TFile::Open("effhists_pim_flat.root");
  //TFile *felec = TFile::Open("effhists_elec_flat.root");

  //TFile *fpim = TFile::Open("pim_effhists.root");
  //TFile *felec = TFile::Open("posit_effhists.root");

  //TFile *fpim = TFile::Open("x100/effhists_pip_flat_2x100.root");
  //TFile *felec = TFile::Open("x100/effhists_elec_flat_0x100.root");

  //TFile *fpim = TFile::Open("effhists_pim_flat_200.root");
  //TFile *felec = TFile::Open("effhists_elec_flat_0x10.root");

  TFile *fpim = TFile::Open("hadd_out/hadd.pi.root");
  TFile *felec = TFile::Open("hadd_out/hadd.e.root");

  bool twod = true;
  double xx[np][100] = {{0.}}, yy[np][100] = {{0.}};
  int npt = 27;
  for (int ipc=0; ipc<npt; ++ipc) {
    const char* nn = twod?Form("prob_cut_%d/eff2d_e_id_indiv",ipc+1):Form("prob_cut_%d/eff1d_mom_e_id_indiv",ipc+1);
    cout << "eff: " << nn << endl;
    TEfficiency *eff_pim = (TEfficiency*) fpim->Get(nn);
    TEfficiency *eff_elec = (TEfficiency*) felec->Get(nn);
    for (int ip=0; ip < np; ++ip) {
      double avg_pim, avg_elec;
      if (twod) {
	avg_pim = avg_eff2d(eff_pim,h_wt[0][ip],2);
	avg_elec = avg_eff2d(eff_elec,h_wt[1][ip],1);
      } else {
	avg_pim = avg_eff(eff_pim,h_wt[0][ip]->ProjectionX());
	avg_elec = avg_eff(eff_elec,h_wt[1][ip]->ProjectionX());
      }
      xx[ip][ipc] = avg_pim;
      yy[ip][ipc] = avg_elec;
    }
  }

  TLegend *tleg = new TLegend(0.45,0.2,0.89,0.5);
  tleg->SetBorderSize(0);
  TMultiGraph *tmg = new TMultiGraph("roc","ROC curve; #varepsilon(#pi^{#pm}); #varepsilon(e^{#pm})");
  TGraph *tg[3];
  for (int ip=0; ip < np; ++ip) {
    tg[ip] = new TGraph(npt,xx[ip],yy[ip]);
    tg[ip]->SetMarkerStyle(20);
    tg[ip]->SetMarkerSize(1.5);
    tg[ip]->SetMarkerColor(col[ip]);
    tleg->AddEntry(tg[ip], plab_tex[ip], "p" );
    tmg->Add(tg[ip],"p");
  }
  TCanvas *tc_roc = new TCanvas("tc_roc","tc_roc");
  tc_roc->cd();
  tmg->Draw("a");
  gPad->Update();
  set_style(tmg->GetHistogram(),0);
  tleg->Draw();
  tc_roc->Print("roc.pdf");

}

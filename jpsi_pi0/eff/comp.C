#include "TEfficiency.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TGAxis.h"

#include <iostream>

TEfficiency* rebin(TFile *f, int rebin, const char *name) {
  cout << "rebinnig " << name << endl;
  TEfficiency *tmp = (TEfficiency*) f->Get(name);
  assert(tmp->GetDimension()==1);
  if (rebin<2) return tmp;
  TH1F* tmpN = (TH1F*)tmp->GetPassedHistogram();
  TH1F* tmpD = (TH1F*)tmp->GetTotalHistogram();
  tmpN->Rebin(rebin);
  tmpD->Rebin(rebin);
  TEfficiency *retval = new TEfficiency(*tmpN,*tmpD);
  retval->SetTitle(tmp->GetTitle());
  retval->SetMarkerStyle(20);
  return retval;
}

TEfficiency* rebin2d(TFile *f, int rebin, const char *name) {
  cout << "rebinnig " << name << endl;
  TEfficiency *tmp = (TEfficiency*) f->Get(name);
  assert(tmp->GetDimension()==2);
  if (rebin<2) return tmp;
  TH2F* tmpN = (TH2F*)tmp->GetPassedHistogram();
  TH2F* tmpD = (TH2F*)tmp->GetTotalHistogram();
  tmpN->RebinX(rebin);
  tmpN->RebinY(rebin);
  tmpD->RebinX(rebin);
  tmpD->RebinY(rebin);
  TEfficiency *retval = new TEfficiency(*tmpN,*tmpD);
  retval->SetTitle(tmp->GetTitle());
  return retval;
}

TH2F* smooth_hist2d(TH2F* h, int nbins= 1000) {
  int nbinx = h->GetXaxis()->GetNbins();
  double xmax = h->GetXaxis()->GetBinUpEdge(nbinx);
  int nbiny = h->GetYaxis()->GetNbins();
  double ymax = h->GetYaxis()->GetBinUpEdge(nbiny);
  TH2F *smooth = new TH2F(Form("%s_smooth",h->GetName()), Form("%s (smoothed)",h->GetTitle()), nbins, 0, xmax, nbins, 0, ymax);
  cout << "xmax= " << xmax << " ymax= " << ymax << endl;
  smooth->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  smooth->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  TH1F *smoothx = (TH1F*) smooth->ProjectionX();
  TH1F *smoothy = (TH1F*) smooth->ProjectionY();
  for (int i = 1; i <= nbins; i++) {
    for (int j = 1; j <= nbins; j++) {
      double xx = smoothx->GetBinCenter(i);
      double yy = smoothy->GetBinCenter(j);
      double eff= h->Interpolate(xx,yy);
      smooth->SetBinContent(i,j,eff);
    }
  }
  return smooth;
}

TEfficiency* smooth_eff2d(TEfficiency *eff)  {
  assert(eff->GetDimension()==2);
  TH2F* num = smooth_hist2d((TH2F*)eff->GetPassedHistogram());
  TH2F* den = smooth_hist2d((TH2F*)eff->GetTotalHistogram());
  TEfficiency *smooth = new TEfficiency(*num,*den);
  smooth->SetTitle(eff->GetTitle());
  return smooth;
}

TEfficiency* slice_eff2d(TEfficiency *eff, int min, int max)  {
  assert(eff->GetDimension()==2);
  TH2F* num = (TH2F*)eff->GetPassedHistogram();
  TH2F* den = (TH2F*)eff->GetTotalHistogram();
  double ymin= num->GetYaxis()->GetBinCenter(min);
  double ymax= num->GetYaxis()->GetBinCenter(max);
  TH1F* num1d = (TH1F*) num->ProjectionX("numx",min,max);
  TH1F* den1d = (TH1F*) den->ProjectionX("denx",min,max);
  TEfficiency *slice = new TEfficiency(*num1d,*den1d);
  slice->SetTitle(Form("%s (%3.0f <#theta < %3.0f)", eff->GetTitle(), ymin, ymax));
  return slice;
}

void comp(int ipc, bool indiv = false){

  //TFile *fpi = TFile::Open("effhists_pim_flat_200.root");
  //TFile *fe = TFile::Open("effhists_elec_flat_0x10.root");
  TFile *fpi = TFile::Open("hadd_out/hadd.pi.root");
  TFile *fe = TFile::Open("hadd_out/hadd.e.root");

  TEfficiency *pi = rebin(fpi, 5, Form("prob_cut_%d/eff1d_the_e_id",ipc));
  TEfficiency *pi_sd = rebin(fpi, 5, Form("prob_cut_%d/eff1d_the_e_id_sd",ipc));
  pi_sd->SetLineColor(2);
  pi_sd->SetMarkerColor(2);
  TEfficiency *pi_indiv = rebin(fpi, 5, Form("prob_cut_%d/eff1d_the_e_id_indiv",ipc));
  TEfficiency *pi_sd_indiv = rebin(fpi, 5, Form("prob_cut_%d/eff1d_the_e_id_sd_indiv",ipc));
  pi_sd_indiv->SetLineColor(2);
  pi_sd_indiv->SetMarkerColor(2);

  TEfficiency *pi_2d = smooth_eff2d(rebin2d(fpi, 5, Form("prob_cut_%d/eff2d_e_id",ipc)));
  //TEfficiency *pi_2d = rebin2d(fpi, 2, Form("prob_cut_%d/eff2d_e_id_sd",ipc));
  TEfficiency *pi_sd_2d = smooth_eff2d(rebin2d(fpi, 1, Form("prob_cut_%d/eff2d_e_id_sd",ipc)));
  TEfficiency *pi_indiv_2d = smooth_eff2d(rebin2d(fpi, 1, Form("prob_cut_%d/eff2d_e_id_indiv",ipc)));
  TEfficiency *pi_sd_indiv_2d = smooth_eff2d(rebin2d(fpi, 1, Form("prob_cut_%d/eff2d_e_id_sd_indiv",ipc)));

  TEfficiency *pi_2d_slice = slice_eff2d(rebin2d(fpi, 1, Form("prob_cut_%d/eff2d_e_id",ipc)), 25, 30);
  TEfficiency *pi_sd_2d_slice = slice_eff2d(rebin2d(fpi, 1, Form("prob_cut_%d/eff2d_e_id_sd",ipc)), 25, 30);
  pi_sd_2d_slice->SetLineColor(2);
  pi_sd_2d_slice->SetMarkerColor(2);
  TEfficiency *pi_indiv_2d_slice = slice_eff2d(rebin2d(fpi, 1, Form("prob_cut_%d/eff2d_e_id_indiv",ipc)), 25, 30);
  TEfficiency *pi_sd_indiv_2d_slice = slice_eff2d(rebin2d(fpi, 1, Form("prob_cut_%d/eff2d_e_id_sd_indiv",ipc)), 25, 30);
  pi_sd_indiv_2d_slice->SetLineColor(2);
  pi_sd_indiv_2d_slice->SetMarkerColor(2);

  TEfficiency *e = rebin(fe, 2, Form("prob_cut_%d/eff1d_the_e_id",ipc));
  TEfficiency *e_sd = rebin(fe, 2, Form("prob_cut_%d/eff1d_the_e_id_sd",ipc));
  e_sd->SetLineColor(2);
  e_sd->SetMarkerColor(2);
  TEfficiency *e_indiv = rebin(fe, 2, Form("prob_cut_%d/eff1d_the_e_id_indiv",ipc));
  TEfficiency *e_sd_indiv = rebin(fe, 2, Form("prob_cut_%d/eff1d_the_e_id_sd_indiv",ipc));
  e_sd_indiv->SetLineColor(2);
  e_sd_indiv->SetMarkerColor(2);

  TEfficiency *e_2d = smooth_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id",ipc)));
  TEfficiency *e_sd_2d = smooth_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id_sd",ipc)));
  TEfficiency *e_indiv_2d = smooth_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id_indiv",ipc)));
  TEfficiency *e_sd_indiv_2d = smooth_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id_sd_indiv",ipc)));

  TEfficiency *e_2d_slice = slice_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id",ipc)), 25, 30);
  TEfficiency *e_sd_2d_slice = slice_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id_sd",ipc)), 25, 30);
  e_sd_2d_slice->SetLineColor(2);
  e_sd_2d_slice->SetMarkerColor(2);
  TEfficiency *e_indiv_2d_slice = slice_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id_indiv",ipc)), 25, 30);
  TEfficiency *e_sd_indiv_2d_slice = slice_eff2d(rebin2d(fe, 1, Form("prob_cut_%d/eff2d_e_id_sd_indiv",ipc)), 25, 30);
  e_sd_indiv_2d_slice->SetLineColor(2);
  e_sd_indiv_2d_slice->SetMarkerColor(2);

  TCanvas *tc2d = new TCanvas("tc2d","tc2d",1600,900);
  tc2d->Divide(2,1);
  tc2d->cd(1);
  if (!indiv)
    pi_2d->Draw("colz");
  else
    pi_indiv_2d->Draw("colz");
  //gPad->SetLogz();
  //gPad->Update();
  //pi_2d->GetPaintedHistogram()->SetMinimum(1e-6);
  //pi_2d->GetPaintedHistogram()->SetMaximum(1e-1);
  tc2d->cd(2);
  if (!indiv)
    e_2d->Draw("colz");
  else
    e_indiv_2d->Draw("colz");

  TCanvas *tc2d_sd = new TCanvas("tc2d_sd","tc2d_sd",1600,900);
  tc2d_sd->Divide(2,1);
  tc2d_sd->cd(1);
  if (!indiv)
    pi_sd_2d->Draw("colz");
  else
    pi_sd_indiv_2d->Draw("colz");
  //gPad->SetLogz();
  //gPad->Update();
  //pi_sd_2d->GetPaintedHistogram()->SetMinimum(1e-6);
  //pi_sd_2d->GetPaintedHistogram()->SetMaximum(1e-1);
  tc2d_sd->cd(2);
  if (!indiv)
    e_sd_2d->Draw("colz");
  else
    e_sd_indiv_2d->Draw("colz");

  TLegend *tl = new TLegend(0.4,0.2,0.99,0.4);
  tl->AddEntry(e,"All dets","pl");
  tl->AddEntry(e_sd,"No dE/dx (STT,MVD)","pl");
  tl->SetHeader(indiv?"   With p_{indiv}>5%":"   No p_{indiv} cut");
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  tl->SetTextSize(0.06);

  TCanvas *tc = new TCanvas("tc","tc",1600,900);
  tc->Divide(2,1);
  tc->cd(1);
  if (!indiv) {
    pi->Draw();
    pi_sd->Draw("same");
    gPad->Update();
    pi->GetPaintedGraph()->GetHistogram()->SetMinimum(1e-5);
    pi->GetPaintedGraph()->GetHistogram()->SetMaximum(1e-1);
  } else {
    pi_indiv->Draw();
    pi_sd_indiv->Draw("same");
    gPad->Update();
    pi_indiv->GetPaintedGraph()->GetHistogram()->SetMinimum(1e-5);
    pi_indiv->GetPaintedGraph()->GetHistogram()->SetMaximum(1e-1);
  }
  gPad->SetLogy();

  tc->cd(2);
  if (!indiv) {
    e->Draw();
    e_sd->Draw("same");
    gPad->Update();
    e->GetPaintedGraph()->GetHistogram()->SetMinimum(0);
    e->GetPaintedGraph()->GetHistogram()->SetMaximum(1);
  } else {
    e_indiv->Draw();
    e_sd_indiv->Draw("same");
    gPad->Update();
    e_indiv->GetPaintedGraph()->GetHistogram()->SetMinimum(0);
    e_indiv->GetPaintedGraph()->GetHistogram()->SetMaximum(1);
  }
  //gPad->SetLogy();
  tl->Draw();

  TCanvas *tc_slice = new TCanvas("tc_slice","tc_slice",1600,900);
  tc_slice->Divide(2,1);
  tc_slice->cd(1);
  if (!indiv) {
    pi_2d_slice->Draw();
    pi_sd_2d_slice->Draw("same");
    gPad->Update();
    pi_2d_slice->GetPaintedGraph()->GetHistogram()->SetMinimum(1e-5);
    pi_2d_slice->GetPaintedGraph()->GetHistogram()->SetMaximum(1e-1);
  } else {
    pi_indiv_2d_slice->Draw();
    pi_sd_indiv_2d_slice->Draw("same");
    gPad->Update();
    pi_indiv_2d_slice->GetPaintedGraph()->GetHistogram()->SetMinimum(1e-5);
    pi_indiv_2d_slice->GetPaintedGraph()->GetHistogram()->SetMaximum(1e-1);
  }
  gPad->SetLogy();

  tc_slice->cd(2);
  if (!indiv) {
    e_2d_slice->Draw();
    e_sd_2d_slice->Draw("same");
    gPad->Update();
    e_2d_slice->GetPaintedGraph()->GetHistogram()->SetMinimum(0);
    e_2d_slice->GetPaintedGraph()->GetHistogram()->SetMaximum(1);
  } else {
    e_indiv_2d_slice->Draw();
    e_sd_indiv_2d_slice->Draw("same");
    gPad->Update();
    e_indiv_2d_slice->GetPaintedGraph()->GetHistogram()->SetMinimum(0);
    e_indiv_2d_slice->GetPaintedGraph()->GetHistogram()->SetMaximum(1);
  }
  //gPad->SetLogy();
  tl->Draw();

}

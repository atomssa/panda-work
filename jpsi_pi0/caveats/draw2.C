#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>

using namespace std;

void draw(int i) {

  gStyle->SetOptStat(0);

  TString path="/vol0/panda/work/jpsi_pi0/grid.out/esim_oct14_binsong_configs/runall";
  TString inFile = path+Form(".%d/bremcorr.root",i);

  Double_t config[22][3] = {
    {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
    {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
    {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
    {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
    {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
    {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

  double nbin = 100;
  double xmin = -0.3;
  double xmax = 0.4.;

  TFile *f = TFile::Open(inFile);
  TTree *t = (TTree*) f->Get("t");

  const int nChMax = 10;
  const int nNeutMax = 100;
  int nChCand;
  int nNeutCand;
  int ch_charge[nChMax];
  float ch_mom_mc[nChMax];
  float ch_mom_rec[nChMax];
  float ch_mom_cor[nChMax];
  float ch_mom_sep[nChMax];
  float ch_mom_mrg[nChMax];
  float ch_mom_stored[nChMax];
  float ch_phi[nChMax];
  float ch_the[nChMax];
  int ch_nphot_sep[nChMax];
  int ch_nphot_mrg[nChMax];
  int ch_is_prim[nChMax];

  t->SetBranchAddress("nch",&nChCand);
  t->SetBranchAddress("neut",&nNeutCand);
  t->SetBranchAddress("ch_charge",&ch_charge);
  t->SetBranchAddress("ch_mom_mc",&ch_mom_mc);
  t->SetBranchAddress("ch_mom_rec",&ch_mom_rec);
  t->SetBranchAddress("ch_mom_cor",&ch_mom_cor);
  t->SetBranchAddress("ch_mom_sep",&ch_mom_sep);
  t->SetBranchAddress("ch_mom_mrg",&ch_mom_mrg);
  t->SetBranchAddress("ch_mom_stored",&ch_mom_stored);
  t->SetBranchAddress("ch_phi",&ch_phi);
  t->SetBranchAddress("ch_the",&ch_the);
  t->SetBranchAddress("ch_nphot_sep",&ch_nphot_sep);
  t->SetBranchAddress("ch_nphot_mrg",&ch_nphot_mrg);
  t->SetBranchAddress("ch_is_prim",&ch_is_prim);

  TH1F* hrec = new TH1F("hrec","hrec",nbin,xmin,xmax);
  hrec->SetLineWidth(2);
  hrec->SetLineColor(1);

  TH1F* hsep = new TH1F("hsep","hsep",nbin,xmin,xmax);
  hsep->SetLineWidth(2);
  hsep->SetLineColor(3);

  TH1F* hcor = new TH1F("hcor","hcor",nbin,xmin,xmax);
  hcor->SetLineWidth(2);
  hcor->SetLineColor(2);

  TH1F* hmrg = new TH1F("hmrg","hmrg",nbin,xmin,xmax);
  hmrg->SetLineWidth(2);
  hmrg->SetLineColor(4);

  TH1F* hsto = new TH1F("hsto","hsto",nbin,xmin,xmax);
  hsto->SetLineWidth(2);
  hsto->SetLineColor(7);

  //TH2F* cor_vs_rec = new TH2F("cor_vs_rec","cor_vs_rec",nbin,0,config[i][0]*2.5,nbin,0,config[i][0]*2.5);
  TH2F* cor_vs_rec = new TH2F("cor_vs_rec","cor_vs_rec",nbin,xmin,xmax,nbin,xmin,xmax);

  TH2F* pmc_vs_cor = new TH2F("pmc_vs_cor","pmc_vs_cor",nbin,xmin,xmax,nbin,0,config[i][0]*2.5);
  TH2F* prec_vs_cor = new TH2F("prec_vs_cor","prec_vs_cor",nbin,xmin,xmax,nbin,0,config[i][0]*2.5);

  TH2F* pmc_vs_rec = new TH2F("pmc_vs_rec","pmc_vs_rec",nbin,xmin,xmax,nbin,0,config[i][0]*2.5);
  TH2F* prec_vs_rec = new TH2F("prec_vs_rec","prec_vs_rec",nbin,xmin,xmax,nbin,0,config[i][0]*2.5);

  int nent = t->GetEntries();
  cout << "nent = " << nent << endl;
  for (int ient = 0; ient < nent; ++ient) {
    if (ient%1000==0) cout << "ient= " << ient << "/" << nent << " nChCand= " << nChCand << endl;
    t->GetEntry(ient);
    for (int ich = 0; ich < nChCand; ++ich) {
      if (!ch_is_prim[ich]) continue;
      float cor = (ch_mom_mc[ich]-ch_mom_cor[ich])/ch_mom_mc[ich];
      float rec = (ch_mom_mc[ich]-ch_mom_rec[ich])/ch_mom_mc[ich];
      hcor->Fill(cor);
      hsep->Fill((ch_mom_mc[ich]-ch_mom_sep[ich])/ch_mom_mc[ich]);
      hmrg->Fill((ch_mom_mc[ich]-ch_mom_mrg[ich])/ch_mom_mc[ich]);
      hrec->Fill(rec);
      hsto->Fill((ch_mom_mc[ich]-ch_mom_stored[ich])/ch_mom_mc[ich]);
      cor_vs_rec->Fill(rec,cor);
      pmc_vs_cor->Fill(cor,ch_mom_mc[ich]);
      prec_vs_cor->Fill(cor,ch_mom_rec[ich]);
      pmc_vs_rec->Fill(rec,ch_mom_mc[ich]);
      prec_vs_rec->Fill(rec,ch_mom_rec[ich]);
    }
  }

  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(2,2);
  tc->cd(1);
  hcor->Draw();
  hsep->Draw("same");
  hrec->Draw("same");
  gPad->Update();
  TLegend *tl1 = new TLegend(0.6,0.4,0.9,0.6);
  tl1->SetFillStyle(0);
  tl1->SetBorderSize(0);
  tl1->SetTextSize(0.07);
  tl1->AddEntry(hcor,"Full corr.", "pl");
  tl1->AddEntry(hsep,"Sep. only", "pl");
  tl1->AddEntry(hrec,"Reco", "pl");
  tl1->Draw();
  TLatex *tt = new TLatex();
  tt->SetTextSize(0.07);
  tt->SetNDC(kTRUE);
  tt->DrawLatex(0.6,0.75,Form("p_{T}= %4.1f",config[i][0]));
  tt->DrawLatex(0.6,0.65,Form("%3.0f < #theta < %3.0f",config[i][1],config[i][2]));

  tc->cd(2);
  hcor->Draw();
  gPad->Update();
  TLegend *tl2 = new TLegend(0.6,0.4,0.9,0.6);
  tl2->SetFillStyle(0);
  tl2->SetBorderSize(0);
  tl2->SetTextSize(0.07);
  tl2->AddEntry(hcor,"Full corr.", "pl");
  tl2->Draw();
  gPad->SetGridx();
  tc->cd(3);
  hcor->Draw();
  hmrg->Draw("same");
  hsep->Draw("same");
  TLegend *tl3 = new TLegend(0.6,0.4,0.9,0.6);
  tl3->SetFillStyle(0);
  tl3->SetBorderSize(0);
  tl3->SetTextSize(0.07);
  tl3->AddEntry(hcor,"Full corr.", "pl");
  tl3->AddEntry(hsep,"Sep. only", "pl");
  tl3->AddEntry(hmrg,"Mrg. only", "pl");
  tl3->Draw();

  return;
  tc->cd(4);
  cor_vs_rec->Draw("colz");

  TCanvas *tc1 = new TCanvas("tc1","tc1");
  tc1->Divide(2,2);
  tc1->cd(1);
  pmc_vs_cor->Draw("colz");
  tc1->cd(2);
  prec_vs_cor->Draw("colz");
  tc1->cd(3);
  pmc_vs_rec->Draw("colz");
  tc1->cd(4);
  prec_vs_rec->Draw("colz");

}

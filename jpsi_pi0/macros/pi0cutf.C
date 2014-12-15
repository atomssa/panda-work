#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

void pi0cutf(int iplab) {
  double lw[3][3] = {
    {0.11,-0.05,+0.00}, /*iplab= 0*/
    {0.12,-0.03,-0.03}, /*iplab= 1*/
    {0.15,-0.06,-0.07}  /*iplab= 2*/
  };

  double up[3][3] = {
    {0.14,0.21,0.07}, /*iplab= 0*/
    {0.12,0.15,0.15}, /*iplab= 1*/
    {0.11,0.10,0.20}  /*iplab= 2*/
  };
  double plab[3] = {5.513, 8., 12.};

  TFile *file = new TFile(Form("test/ana_jpsi_brem_plab%3.1f.root",plab[iplab]));
  TH2F* h2 = (TH2F*)file->Get("gg/h_oa_gg_avg_e_g_truepi0_rec")->Clone("h2");

  h2->GetXaxis()->SetRangeUser(0,1.0);
  TCanvas*tc0 = new TCanvas("tc0","tc0");
  h2->Draw("colz");

  TF1 *flow = new TF1("flow","[2]+([0]/(x-[1]))",lw[iplab][1],3.0);
  flow->SetLineColor(2);
  flow->SetParameter(0,lw[iplab][0]);
  flow->SetParameter(1,lw[iplab][1]);
  flow->SetParameter(2,lw[iplab][2]);
  flow->Draw("same");

  TF1 *fup = new TF1("fup","[2]+([0]/(x-[1]))",up[iplab][1],3.0);
  fup->SetLineColor(1);
  fup->SetParameter(0,up[iplab][0]);
  fup->SetParameter(1,up[iplab][1]);
  fup->SetParameter(2,up[iplab][2]);
  fup->Draw("same");

  //TCanvas *tc = new TCanvas("tc","tc");
  //tc->cd();
  //fup->Draw();
  //flow->Draw("same");

}

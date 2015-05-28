#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>

using namespace std;

static const double zdet[7] = {0,6,8.5,12,17,25,200};

void padSetup(TCanvas *tc, TPad *pads[], int ii) {

  double a= 0.15;
  double b= 0.05;
  double _a = a/(1-a);
  double _b = b/(1-b);
  double w = 1/(3+_a+_b);
  double yl = ii<3?0.53:0.0;
  double yh = ii<3?1.0:0.53;
  double xl[3] = {0.0, (_a+1)*w, (_a+2)*w-0.003};
  double xh[3] = {(_a+1)*w, (_a+2)*w, (_a+_b+3)*w};

  cout << "p " << ii << "(" << xl[ii%3] << ", " << yl << ", " << xh[ii%3] << ", " << yh << ")" << " w= " << xh[ii%3] - xl[ii%3] << endl;
  pads[ii] = new TPad(Form("pads%d",ii),Form("pads%d",ii), xl[ii%3], yl,  xh[ii%3], yh);
  pads[ii]->SetBorderSize(0);
  tc->cd(0);

  pads[ii]->Draw();
  double epsilon=1e-9;
  if ((ii%3)==0) {pads[ii]->SetRightMargin(epsilon); pads[ii]->SetLeftMargin(0.2);}
  if ((ii%3)==1) {pads[ii]->SetLeftMargin(epsilon); pads[ii]->SetRightMargin(epsilon); }
  if ((ii%3)==2) {pads[ii]->SetLeftMargin(epsilon); pads[ii]->SetRightMargin(0.1);}

  if (ii<3) {
    pads[ii]->SetTopMargin(0.05);
    pads[ii]->SetBottomMargin(epsilon);
  } else {
    pads[ii]->SetTopMargin(epsilon);
    pads[ii]->SetBottomMargin(0.15);
  }
  pads[ii]->SetTicks(0,1);
  pads[ii]->cd();
}

void raddep(int iconfig = 2){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  //int iconfig = 8; double max = 0.2;
  double max = 0.2;

  static const int nconfig = 25;
  Double_t config[nconfig][3] = {
    /*0*/{0.15, 30.0, 45.0}, /*1*/{0.25, 30.0, 45.0}, /*2*/{0.5, 30.0, 45.0}, /*3*/{0.5, 45.0, 60.0}, /*4*/{0.5, 60.0, 75.0},
    /*5*/{0.5, 75.0, 90.0}, /*6*/{0.5, 90.0, 105.0}, /*7*/{0.5, 105.0, 120.0}, /*8*/{1.0, 30.0, 45.0}, /*9*/{1.0, 45.0, 60.0},
    /*10*/{1.0, 60.0, 75.0}, /*11*/{1.0, 75.0, 90.0}, /*12*/{1.5, 30.0, 45.0}, /*13*/{1.5, 45.0, 60.0}, /*14*/{2.0, 30.0, 45.0},
    /*15*/{0.15, 10.0, 20.0}, /*16*/{0.25, 10.0, 20.0}, /*17*/{0.5, 10.0, 20.0}, /*18*/{1.0, 10.0, 20.0}, /*19*/{1.5, 10.0, 20.0},
    /*20*/{2.0, 10.0, 20.0}, /*21*/{2.5, 10.0, 20.0}, /*22*/{0.2, 30.0, 45.0}, /*23*/{0.2, 90.0, 105.0}, /*24*/{0.2, 10.0, 20.0}};

  TFile *f = TFile::Open(Form("../grid.out/esim_oct14_binsong_configs/all.ibs.xr-0.1_0.3/bremcorr.all.ibs.cfg.%d_hists.root",iconfig));

  TH1F* h_rec_1brem[6];
  TH1F* h_out_1brem[6];

  TCanvas *tc[6];
  TLatex *tt = new TLatex();
  tt->SetTextSize(0.06);
  tt->SetNDC(kTRUE);

  TCanvas *tcall = new TCanvas("tcall","tcall");
  //tcall->cd();
  tcall->Divide(3,2);
  TPad *pads[6];

  for (int i=0; i<6; ++i) {
    //tc[i] = new TCanvas(Form("tc_rad_dep_%d",i),Form("tc_rad_dep_%d",i),1000,1000);
    //padSetup(tcall, pads, i);
    tcall->cd(i+1);

    h_out_1brem[i] = (TH1F*) f->Get(Form("h_out_1brem_%d",i));
    //h_out_1brem[i]->Scale(1./h_out_1brem[i]->GetEntries());
    //h_out_1brem[i]->SetTitle(Form(";%s",h_out_1brem[i]->GetXaxis()->GetTitle()));
    h_out_1brem[i]->SetTitle(";(p - p_{KF})/p");
    h_out_1brem[i]->SetLineColor(2);

    h_rec_1brem[i] = (TH1F*) f->Get(Form("h_rec_1brem_%d",i));
    //h_rec_1brem[i]->Scale(1./h_rec_1brem[i]->GetEntries());
    //h_rec_1brem[i]->SetTitle(Form(";%s",h_rec_1brem[i]->GetXaxis()->GetTitle()));
    h_rec_1brem[i]->SetTitle(";(p - p_{KF})/p");
    h_rec_1brem[i]->SetLineColor(4);

    h_rec_1brem[i]->Rebin(4);
    h_out_1brem[i]->Rebin(4);

    //h_rec_1brem[i]->SetMaximum(max);
    //h_rec_1brem[i]->Draw(i==5?"":"same");

    if (h_rec_1brem[i]->GetMaximum()>h_out_1brem[i]->GetMaximum()) {
      h_rec_1brem[i]->Draw();
      h_out_1brem[i]->Draw("same");
    } else {
      h_out_1brem[i]->Draw();
      h_rec_1brem[i]->Draw("same");
    }

    gPad->SetGridx();

    tt->SetTextSize(0.06);
    tt->SetTextColor(1);
    bool fwd = config[iconfig][1] < 15.0 && 15.0<config[iconfig][2];
    int s= fwd ? zdet[i] : i*7;
    int e= fwd ? zdet[i+1] : (i+1)*7;
    const char* var = (fwd?"Z_{True}":"R_{True}");
    //double off=(i==2||i==5)?-0.1:0.0;
    double off=-0.11;
    tt->DrawLatex(off+0.46,0.85,Form("%d < %s(cm) < %d", s, var, e));
    tt->DrawLatex(off+0.53,0.75,Form("p_{T}= %4.1f GeV/c",config[iconfig][0]));
    tt->DrawLatex(off+0.53,0.65,Form("%3.0f#circ < #theta < %3.0f#circ",config[iconfig][1],config[iconfig][2]));
    tt->SetTextColor(4);
    tt->SetTextSize(0.08);
    tt->DrawLatex(off+0.63,0.55,Form("p = p_{MC}"));
    tt->SetTextColor(2);
    tt->DrawLatex(off+0.58,0.43,Form("p = p_{MC} - E_{#gamma}"));
    //tc[i]->Print(Form("rad_dep_config%d_r%d.pdf",iconfig,i));
  }

}

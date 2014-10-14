void set_style(TH1* h, int col=4) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);

  if (col>0) {
    h->SetLineWidth(3);
    h->SetLineColor(col);
  }
}
void figs(int page) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);

  // 0 - bg, 1 - Sig
  TFile *f_brem[2];
  f_brem[0] = TFile::Open("hists/ana_bg_brem.root");
  f_brem[1] = TFile::Open("hists/ana_jpsi_brem.root");

  TCanvas *tc;
  TLegend *tl[2];
  tl[0] = new TLegend(0.35,0.6,0.9,0.8);
  tl[0]->SetFillStyle(0);
  tl[0]->SetBorderSize(0);
  tl[1] = new TLegend(0.35,0.6,0.9,0.8);
  tl[1]->SetFillStyle(0);
  tl[1]->SetBorderSize(0);

  switch(page) {

  case 1:
    {
      tc = new TCanvas("tc_oa_vs_m","tc_oa_vs_m",1400,1000);
      tc->Divide(2,1);
      TH1F *h_oa_vs_m[2];
      for (int i=0; i<2; ++i) {
	h_oa_vs_m[i] = (TH1F*) f_brem[i]->Get("h_oa_vs_mass_gg_epair");
	tc->cd(i+1);
	h_oa_vs_m[i]->Draw("colz");
      }
      break;
    }

  case 2:
    {
      tc = new TCanvas("tc_oa_vs_mgg","tc_oa_vs_mgg",1400,1000);
      tc->Divide(2,1);
      TH1F *h_oa_vs_mgg[2];
      for (int i=0; i<2; ++i) {
	h_oa_vs_mgg[i] = (TH1F*) f_brem[i]->Get("h_oa_gg_epair_vs_mass_gg");
	tc->cd(i+1);
	h_oa_vs_mgg[i]->Draw("colz");
      }
      break;
    }

  case 3:
    {
      tc = new TCanvas("tc_m_vs_mgg","tc_m_vs_mgg",1400,1000);
      tc->Divide(2,1);
      TH1F *h_m_vs_mgg[2];
      for (int i=0; i<2; ++i) {
	h_m_vs_mgg[i] = (TH1F*) f_brem[i]->Get("h_mass_gg_epair_vs_mass_gg");
	tc->cd(i+1);
	h_m_vs_mgg[i]->Draw("colz");
      }
      break;
    }

  case 4:
    {
      tc = new TCanvas("tc_m_vs_mgg","tc_m_vs_mgg",1400,1000);
      tc->Divide(2,1);
      TH1F *h_m_vs_mgg_btb[2];
      for (int i=0; i<2; ++i) {
	h_m_vs_mgg_btb[i] = (TH1F*) f_brem[i]->Get("h_oa_vs_mass_gg_epair_btb");
	tc->cd(i+1);
	h_m_vs_mgg_btb[i]->Draw("colz");
      }
      break;
    }

  case 5:
    {
      tc = new TCanvas("tc_m_vs_mgg","tc_m_vs_mgg",1400,1000);
      tc->Divide(2,1);
      TH1F *h_mgg_btb[2];
      for (int i=0; i<2; ++i) {
	h_mgg_btb[i] = (TH1F*) f_brem[i]->Get("h_m_gg_btb");
	tc->cd(i+1);
	h_mgg_btb[i]->Draw("colz");
      }
      break;
    }

  case 6:
    {
      tc = new TCanvas("tc_m_vs_mgg","tc_m_vs_mgg",1400,1000);
      tc->Divide(2,1);
      TH1F *h_m_vs_mgg_cts[2];
      for (int i=0; i<2; ++i) {
	h_m_vs_mgg_cts[i] = (TH1F*) f_brem[i]->Get("h_oa_vs_mass_gg_epair_cts");
	tc->cd(i+1);
	h_m_vs_mgg_cts[i]->Draw("colz");
      }
      break;
    }

  case 7:
    {
      tc = new TCanvas("tc_m_vs_mgg","tc_m_vs_mgg",1400,1000);
      tc->Divide(2,1);
      TH1F *h_mgg_cts[2];
      for (int i=0; i<2; ++i) {
	h_mgg_cts[i] = (TH1F*) f_brem[i]->Get("h_m_gg_cts");
	tc->cd(i+1);
	h_mgg_cts[i]->Draw("colz");
      }
      break;
    }

  case 8:
    {
      tc = new TCanvas("tc_prob_4c_epem","tc_prob_4c_epem",1400,1000);
      tc->Divide(2,1);
      TH1F *h_prob_4c_cts[2];
      for (int i=0; i<2; ++i) {
	h_prob_4c_cts[i] = (TH1F*) f_brem[i]->Get("h_4c_prob_vs_m_cts_epem_4c");
	tc->cd(i+1);
	h_prob_4c_cts[i]->Draw("colz");
	//gPad->SetLogy();
      }
      break;
    }

  case 9:
    {
      tc = new TCanvas("tc_prob_4c_epem","tc_prob_4c_epem",1400,1000);
      tc->Divide(2,1);
      TH1F *h_prob_4c_cts[2];
      for (int i=0; i<2; ++i) {
	h_prob_4c_cts[i] = (TH1F*) f_brem[i]->Get("h_4c_prob_cts_epempi0");
	tc->cd(i+1);
	h_prob_4c_cts[i]->Draw("colz");
	gPad->SetLogy();
      }
      break;
    }


  }


  TText *tt = new TText();
  tc->cd(1);
  tt->SetTextSize(0.08);
  tt->DrawTextNDC(0.44,0.9,"Background");
  tc->cd(2);
  tt->DrawTextNDC(0.65,0.9,"Signal");


  tc->Print(Form("figs/2014.10.13/fig_p%d.pdf",page));

}

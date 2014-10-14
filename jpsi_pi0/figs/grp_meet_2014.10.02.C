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
  TFile *f_raw[2], *f_brem[2];
  f_raw[0] = TFile::Open("hists/grp_meet_2014.10.02/ana_bg_raw.root");
  f_brem[0] = TFile::Open("hists/grp_meet_2014.10.02/ana_bg_brem.root");
  f_raw[1] = TFile::Open("hists/grp_meet_2014.10.02/ana_jpsi_raw.root");
  f_brem[1] = TFile::Open("hists/grp_meet_2014.10.02/ana_jpsi_brem.root");

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
      TH1F *ng[2];
      TH1F *ng_tr[2];
      for (int i=0; i<2; ++i) {
	ng[i] = (TH1F*) f_raw[i]->Get("h_num_g");
	set_style(ng[i]);
	ng_tr[i] = (TH1F*) f_raw[i]->Get("h_num_g_tr");
	set_style(ng_tr[i],2);
	tl[i]->AddEntry(ng[i],"All","l");
	tl[i]->AddEntry(ng_tr[i],"Truth Match","l");
      }
      tc = new TCanvas("tc_ng","tc_ng",1400,1000);
      tc->Divide(2,1);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	ng_tr[i]->Draw();
	ng[i]->Draw("same");
	tl[i]->Draw();
      }
      break;
    }

  case 2:
    {
      TH1F *eg[2];
      TH1F *eg_tr[2];
      for (int i=0; i<2; ++i) {
	eg[i] = (TH1F*) f_raw[i]->Get("h_e_g");
	set_style(eg[i]);
	eg_tr[i] = (TH1F*) f_raw[i]->Get("h_e_g_tr");
	set_style(eg_tr[i],2);
	tl[i]->AddEntry(eg[i],"All","l");
	tl[i]->AddEntry(eg_tr[i],"Truth Match","l");
      }
      tc = new TCanvas("tc_eg","tc_eg",1400,1000);
      tc->Divide(2,1);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	eg[i]->Draw();
	eg_tr[i]->Draw("same");
	gPad->SetLogy();
	tl[i]->Draw();
      }
      break;
    }

  case 3:
    {
      TH1F *nch[2];
      TH1F *nch_tr[2];
      for (int i=0; i<2; ++i) {
	nch[i] = (TH1F*) f_raw[i]->Get(i==0?"h_num_pipm":"h_num_epm");
	set_style(nch[i]);
	nch_tr[i] = (TH1F*) f_raw[i]->Get(i==0?"h_num_pipm_tr":"h_num_epm_tr");
	set_style(nch_tr[i],2);
	tl[i]->AddEntry(nch[i],"All","l");
	tl[i]->AddEntry(nch_tr[i],"Truth Match","l");
      }
      tc = new TCanvas("tc_nch","tc_nch",1400,1000);
      tc->Divide(2,1);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	nch_tr[i]->Draw();
	nch[i]->Draw("same");
	tl[i]->Draw();
      }
      break;
    }

  case 4:
    {
      TH2F *p_th_ch[2];
      TH2F *p_th_ch_tr[2];
      for (int i=0; i<2; ++i) {
	p_th_ch[i] = (TH2F*) f_raw[i]->Get(i==0?"h_mom_the_pipm":"h_mom_the_epm");
	set_style(p_th_ch[i],-1);
	p_th_ch_tr[i] = (TH2F*) f_raw[i]->Get(i==0?"h_mom_the_pipm_tr":"h_mom_the_epm_tr");
	set_style(p_th_ch_tr[i],-1);
	//tl[i]->AddEntry(p_th_ch[i],"All","l");
	//tl[i]->AddEntry(p_th_ch_tr[i],"Truth Match","l");
      }
      tc = new TCanvas("tc_p_th_ch","tc_p_th_ch",1400,1000);
      tc->Divide(2,1);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	p_th_ch_tr[i]->Draw("colz");
	//p_th_ch[i]->Draw("colz,same");
	//tl[i]->Draw();
      }
      break;
    }

  case 5:
    {
      TH1F *all[2];
      TH1F *ftm[2];
      TH1F *nst[2];
      for (int i=0; i<2; ++i) {
	all[i] = (TH1F*) f_raw[i]->Get("hpi0m_all");
	all[i]->GetXaxis()->SetTitle("M_{#gamma-#gamma}[GeV/c^{2}]");
	set_style(all[i],1);
	tl[i]->AddEntry(all[i],"All #gamma-#gamma pairs", "l");
	ftm[i] = (TH1F*) f_raw[i]->Get("hpi0m_ftm");
	set_style(ftm[i],2);
	tl[i]->AddEntry(ftm[i],"Truth Matched #gamma-#gamma pairs", "l");
	nst[i] = (TH1F*) f_raw[i]->Get("h_m_pi0n");
	set_style(nst[i]);
	tl[i]->AddEntry(all[i],"#gamma-#gamma nearest to M_{#pi^{0}}^{PDG}", "l");
      }
      tc = new TCanvas("tc_pi0m","tc_pi0m",1400,1000);
      tc->Divide(2);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	all[i]->Draw();
	gPad->SetLogy();
	all[i]->SetMinimum(5);
	ftm[i]->Draw("same");
	nst[i]->Draw("same");
	tl[i]->Draw();
      }
      break;
    }

  case 6:
    {
      tl[1] = new TLegend(0.15,0.6,0.8,0.8);
      tl[1]->SetFillStyle(0);
      tl[1]->SetBorderSize(0);
      TH1F *all[2];
      TH1F *ftm[2];
      TH1F *ftm_jpsi[2];
      cout << "fetching" << endl;
      for (int i=0; i<2; ++i) {
	all[i] = (TH1F*) f_raw[i]->Get(i==0?"h_m_pippim":"h_m_epem");
	set_style(all[i],1);
	tl[i]->AddEntry(all[i],"All pairs", "l");
	ftm[i] = (TH1F*) f_raw[i]->Get(i==0?"h_m_pippim_tr":"h_m_epem_tr");
	set_style(ftm[i],2);
	tl[i]->AddEntry(ftm[i],"Truth Matched pairs", "l");
	ftm_jpsi[i] = (TH1F*) f_raw[i]->Get("hjpsim_ftm");
	set_style(ftm_jpsi[i],4);
	if (i==1) tl[i]->AddEntry(ftm_jpsi[i],"Truth Matched J/#psi", "l");
      }
      cout << "done fetching" << endl;
      tc = new TCanvas("tc_ch_pair_m","tc_ch_pair_m",1400,1000);
      tc->Divide(2);

      for (int i=0; i<2; ++i) {
	cout << "pad " << i << endl;
	tc->cd(1+i);
	all[i]->Draw();
	ftm[i]->Draw("same");
	if (i==1) ftm_jpsi[i]->Draw("same");
	tl[i]->Draw();
	cout << "pad " << i << "done " << endl;
      }

      TLatex *tt0 = new TLatex();
      tt0->SetNDC();
      tt0->SetTextSize(0.04);
      tc->cd(1);
      tt0->DrawLatex(0.45,0.8,"#pi^{#pm} mass hypothesis");
      tc->cd(2);
      tt0->DrawLatex(0.6,0.8,"e^{#pm} mass hypothesis");

      break;
    }

  case 7:
    {
      tl[1] = new TLegend(0.15,0.6,0.8,0.8);
      tl[1]->SetFillStyle(0);
      tl[1]->SetBorderSize(0);
      TH1F *ftm[2];
      TH1F *ftm_brem[2];
      for (int i=0; i<2; ++i) {
	ftm[i] = (TH1F*) f_raw[i]->Get(i==0?"h_m_epem":"hjpsim_ftm");
	set_style(ftm[i],1);
	tl[i]->AddEntry(ftm[i],"No Brem Correction", "l");
	ftm_brem[i] = (TH1F*) f_brem[i]->Get(i==0?"h_m_epem":"hjpsim_ftm");
	set_style(ftm_brem[i],4);
	tl[i]->AddEntry(ftm_brem[i],"Brem Corrected", "l");
      }
      tc = new TCanvas("tc_brem_comp_m","tc_brem_comp_m",1400,1000);
      tc->Divide(2);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	if (i==0) {
	  ftm[i]->Draw();
	  ftm_brem[i]->Draw("same");
	} else {
	  ftm_brem[i]->Draw();
	  ftm[i]->Draw("same");
	}
	tl[i]->Draw();
      }

      TLatex *tt0 = new TLatex();
      tt0->SetNDC();
      tt0->SetTextSize(0.04);
      tc->cd(1);
      tt0->DrawLatex(0.45,0.8,"e^{#pm} mass hypothesis");
      tc->cd(2);
      tt0->DrawLatex(0.6,0.8,"e^{#pm} mass hypothesis");

      break;
    }


  case 8:
    {
      TH1F *chi2[2];
      TH1F *prob[2];
      for (int i=0; i<2; ++i) {
	chi2[i] = (TH1F*) f_raw[i]->Get("h_4c_chi2_epempi0");
	set_style(chi2[i],i==0?4:2);
	prob[i] = (TH1F*) f_raw[i]->Get("h_4c_prob_epempi0");
	set_style(prob[i],i==0?4:2);
      }

      tl[0]->AddEntry(chi2[0],"Background");
      tl[0]->AddEntry(chi2[1],"Signal");
      tl[1]->AddEntry(prob[0],"Background");
      tl[1]->AddEntry(prob[1],"Signal");

      tc = new TCanvas("tc_4c_chi2","tc_4c_chi2",1400,1000);
      tc->Divide(2);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	if (i==0) {
	  chi2[1]->DrawNormalized();
	  chi2[0]->DrawNormalized("same");
	} else {
	  gPad->SetLogy();
	  prob[0]->DrawNormalized();
	  prob[1]->DrawNormalized("same");
	}
	tl[i]->Draw();
      }

      TLatex *tt0 = new TLatex();
      tt0->SetNDC();
      tt0->SetTextSize(0.04);
      tc->cd(1);
      tt0->DrawLatex(0.45,0.85,"No Brem. Correction");
      tt0->DrawLatex(0.45,0.8,"e^{#pm} mass hypothesis");
      tc->cd(2);
      tt0->DrawLatex(0.45,0.85,"No Brem. Correction");
      tt0->DrawLatex(0.45,0.8,"e^{#pm} mass hypothesis");
      break;
    }

  case 9:
    {
      TH1F *chi2[2];
      TH1F *prob[2];
      for (int i=0; i<2; ++i) {
	chi2[i] = (TH1F*) f_brem[i]->Get("h_4c_chi2_epempi0");
	set_style(chi2[i],i==0?4:2);
	prob[i] = (TH1F*) f_brem[i]->Get("h_4c_prob_epempi0");
	set_style(prob[i],i==0?4:2);
      }

      tl[0]->AddEntry(chi2[0],"Background");
      tl[0]->AddEntry(chi2[1],"Signal");
      tl[1]->AddEntry(prob[0],"Background");
      tl[1]->AddEntry(prob[1],"Signal");

      tc = new TCanvas("tc_4c_chi2","tc_4c_chi2",1400,1000);
      tc->Divide(2);
      for (int i=0; i<2; ++i) {
	tc->cd(1+i);
	if (i==0) {
	  chi2[1]->DrawNormalized();
	  chi2[0]->DrawNormalized("same");
	} else {
	  gPad->SetLogy();
	  prob[0]->DrawNormalized();
	  prob[1]->DrawNormalized("same");
	}
	tl[i]->Draw();
      }

      TLatex *tt0 = new TLatex();
      tt0->SetNDC();
      tt0->SetTextSize(0.04);
      tc->cd(1);
      tt0->DrawLatex(0.45,0.85,"With Brem. Correction");
      tt0->DrawLatex(0.45,0.8,"e^{#pm} mass hypothesis");
      tc->cd(2);
      tt0->DrawLatex(0.45,0.85,"With Brem. Correction");
      tt0->DrawLatex(0.45,0.8,"e^{#pm} mass hypothesis");

      break;
    }

  case 10:
    {

      TH1F* sig = (TH1F*) f_brem[1]->Get("h_m_epem")->Clone("sig");
      TH1F* bg = (TH1F*) f_brem[0]->Get("h_m_epem")->Clone("bg");
      TH1F* tot = (TH1F*) f_brem[0]->Get("h_m_epem")->Clone("tot");
      tot->Add((TH1F*) f_brem[1]->Get("h_m_epem"));

      set_style(sig,2);
      set_style(bg,4);
      set_style(tot,1);

      tc = new TCanvas("tc_m","tc_m", 1400,1000);
      tc->cd();
      tot->Draw();
      sig->Draw("same");
      bg->Draw("same");

      break;
    }


  }

  if (page!=8&&page!=9&&page!=10) {
    TText *tt = new TText();
    tc->cd(1);
    tt->SetTextSize(0.08);
    tt->DrawTextNDC(0.44,0.83,"Background");
    tc->cd(2);
    tt->DrawTextNDC(0.65,0.83,"Signal");
  }

  tc->Print(Form("grp_meet/fig_p%d.pdf",page));


}

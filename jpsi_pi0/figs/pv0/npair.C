void npair() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  static const int nstep = 5;
  TH1F* sig_nep[nplab][nstep], *sig_ngg[nplab][nstep];
  TH1F* bg_nep[nplab][nstep], *bg_ngg[nplab][nstep];
  TFile *fsig[nplab], *fbg[nplab];
  for (int iplab = 0; iplab < nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("%s/hists/pcm.jul.2015/anav2_jpsi_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"),plab[iplab]));
    fbg[iplab] = TFile::Open(Form("%s/hists/pcm.jul.2015/anav2_pip_pim_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"),plab[iplab]));
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
  sig_s0_s1->Print(Form("%s/figs/2015.09.15/num_ep_pair_vs_cut_sig.pdf",bdir));
  sig_s0_s1_noarrow->Print(Form("%s/figs/2015.09.15/num_ep_pair_vs_cut_sig_noarrow.pdf",bdir));

}

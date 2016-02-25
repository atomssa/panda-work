void ana() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  const int nstep = 5;
  const int ntbin = 12;

  double mmin = 2.96, mmax = 3.22;
  int immin = 0, immax = 0;
  double x[nstep] = {1, 2 ,3, 4, 5};

  double yield_sig[nplab][nstep]= {{0.0}}, yield_bg[nplab][nstep]= {{0.0}};
  double eff_sig[nplab][nstep]= {{0.0}}, eff_bg[nplab][nstep]= {{0.0}};
  double stob[nplab][nstep]= {{0.0}};

  TGraph *tg_sig_eff[nplab], *tg_bg_eff[nplab];
  TGraph *tg_sig_yield[nplab], *tg_bg_yield[nplab];
  TGraph *tg_stob[nplab];

  TMultiGraph *tmg_yield[nplab], *tmg_eff[nplab];

  TCanvas *tc_mep_all_steps[nplab];
  TCanvas *tc_mep_steps[nplab][nstep];

  TFile *fsig[nplab], *fbg[nplab];
  TH1F* hmsg[nplab][nstep], *hmbg[nplab][nstep], *hmfg[nplab][nstep];
  TH1F* hmsg_tb[nplab][ntbin], *hmbg_tb[nplab][ntbin], *hmfg_tb[nplab][ntbin];

  TLatex *tl[6][nplab];
  for (int ii = 0; ii < 6; ++ii) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tl[ii][iplab] = new TLatex();
      tl[ii][iplab]->SetNDC(true);
      //tl[ii][iplab]->SetTextColor(ii==0?2:1);
      if (ii==0) tl[ii][iplab]->SetTextSize((iplab==0?1.5:2.0)*tl[ii][iplab]->GetTextSize());
      if (ii==1) tl[ii][iplab]->SetTextSize(0.8*(iplab==0?1.5:2.0)*tl[ii][iplab]->GetTextSize());
      if (ii==2) tl[ii][iplab]->SetTextSize(1.8*tl[ii][iplab]->GetTextSize());
      if (ii==3) tl[ii][iplab]->SetTextSize(0.8*1.8*tl[ii][iplab]->GetTextSize());
    }
  }

  TLegend *leg[3];
  leg[0] = new TLegend(0.2,0.5,0.8,0.7);
  leg[1] = new TLegend(0.15,0.65,0.7,0.75);
  leg[2] = new TLegend(0.2,0.25,0.8,0.45);
  for (int jj = 0; jj < 3; ++jj) {
    leg[jj]->SetFillStyle(0);
    leg[jj]->SetBorderSize(0);
    leg[jj]->SetTextSize(0.07);
  }

  for (int iplab=0; iplab<nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("%s/hists/pcm.jul.2015/anav2_jpsi_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));
    fbg[iplab] = TFile::Open(Form("%s/hists/pcm.jul.2015/anav2_pip_pim_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));

    for (int istep=0; istep<nstep; ++istep) {

      hmfg[iplab][istep] = (TH1F*) fsig[iplab]->Get(Form("hmep_%d",istep))->Clone(Form("hmep_fg_p%d_s%d",iplab,istep));
      hmsg[iplab][istep] = (TH1F*) fsig[iplab]->Get(Form("hmep_%d",istep))->Clone(Form("hmep_sig_p%d_s%d",iplab,istep));
      hmbg[iplab][istep] = (TH1F*) fbg[iplab]->Get(Form("hmep_%d",istep))->Clone(Form("hmep_bg_p%d_s%d",iplab,istep));
      hmfg[iplab][istep]->Add(hmbg[iplab][istep]);
      hmfg[iplab][istep]->SetTitle(";M_{inv}[GeV/c^{2}]");
      if (immin==0||immax==0) {
	immin = hmsg[iplab][istep]->GetXaxis()->FindBin(mmin);
	immax = hmsg[iplab][istep]->GetXaxis()->FindBin(mmax);
      }

      // Hackery to combine the pid steps for SIG and BG
      if (istep==0) {
	yield_sig[iplab][istep] = hmsg[iplab][istep]->Integral(immin, immax);
	yield_bg[iplab][istep] = hmbg[iplab][istep]->Integral(immin, immax);
	eff_sig[iplab][istep] = istep==0?100.:100*yield_sig[iplab][istep]/yield_sig[iplab][0];
	eff_bg[iplab][istep] = istep==0?100.:100*yield_bg[iplab][istep]/yield_bg[iplab][0];
	stob[iplab][istep] = yield_sig[iplab][istep]/yield_bg[iplab][istep];
      } else if (istep>1) {
	yield_sig[iplab][istep-1] = hmsg[iplab][istep]->Integral(immin, immax);
	yield_bg[iplab][istep-1] = hmbg[iplab][istep]->Integral(immin, immax);
	eff_sig[iplab][istep-1] = istep==0?100.:100*yield_sig[iplab][istep-1]/yield_sig[iplab][0];
	eff_bg[iplab][istep-1] = istep==0?100.:100*yield_bg[iplab][istep-1]/yield_bg[iplab][0];
	stob[iplab][istep-1] = yield_sig[iplab][istep-1]/yield_bg[iplab][istep-1];
      }

      set_style(hmfg[iplab][istep], 1, 2, true);
      set_style(hmsg[iplab][istep], 2, 2, true);
      set_style(hmbg[iplab][istep], 4, 2, true);

      hmfg[iplab][istep]->GetXaxis()->SetRangeUser(0,4.0);
      hmsg[iplab][istep]->GetXaxis()->SetRangeUser(0,4.0);
      hmbg[iplab][istep]->GetXaxis()->SetRangeUser(0,4.0);

    }

    tc_mep_all_steps[iplab] = new TCanvas(Form("tc_mep_all_steps_p%d",iplab),Form("tc_mep_all_steps_p%d",iplab));
    tc_mep_all_steps[iplab]->Divide(3,2);
    for (int istep=0; istep<nstep; ++istep) {
      tc_mep_steps[iplab][istep] = new TCanvas(Form("tc_mep_steps_p%d_s%d",iplab,istep),Form("tc_mep_steps_p%d_s%d",iplab,istep));
      draw_fg_sig_bg(tc_mep_all_steps[iplab], istep+1, hmfg[iplab][istep], hmsg[iplab][istep], hmbg[iplab][istep], true);
      draw_fg_sig_bg(tc_mep_steps[iplab][istep], 0, hmfg[iplab][istep], hmsg[iplab][istep], hmbg[iplab][istep], true);
      //tc_mep_steps[iplab][istep]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_mep_steps[iplab][istep]->GetName()));
    }
    //tc_mep_all_steps[iplab]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_mep_all_steps[iplab]->GetName()));

    tmg_yield[iplab] = new TMultiGraph(Form("tmg_yield_p%d",iplab),Form("tmg_yield_p%d",iplab));
    tg_sig_yield[iplab] = new TGraph(nstep-1,x,yield_sig[iplab]);
    set_style(tg_sig_yield[iplab], 1, 20, 3);
    //set_style(tg_sig_yield[iplab], iplab==0?1:(iplab==1?2:4), 20, 1);
    tg_sig_yield[iplab]->SetMinimum(0);
    tmg_yield[iplab]->Add(tg_sig_yield[iplab],"p");
    tg_bg_yield[iplab] = new TGraph(nstep-1,x,yield_bg[iplab]);
    //set_style(tg_bg_yield[iplab], iplab==0?1:(iplab==1?2:4), 24, 1);
    set_style(tg_bg_yield[iplab], 1, 24, 3);
    tg_bg_yield[iplab]->SetMinimum(0);
    tmg_yield[iplab]->Add(tg_bg_yield[iplab],"p");
    for (int jj = 0; jj < 3; ++jj) {
      if (iplab==0) {
	leg[jj]->SetHeader(jj==0?"Raw yields":(jj==1?"Yields relative to no cuts":"S/B"));
	leg[jj]->AddEntry(tg_sig_yield[iplab],"Signal","p");
	leg[jj]->AddEntry(tg_bg_yield[iplab],"Background","p");
      }
    }

    tmg_eff[iplab] = new TMultiGraph(Form("tmg_eff_p%d",iplab),Form("tmg_eff_p%d",iplab));
    tg_sig_eff[iplab] = new TGraph(nstep-1,x,eff_sig[iplab]);
    //set_style(tg_sig_eff[iplab], iplab==0?1:(iplab==1?2:4), 20, 1);
    set_style(tg_sig_eff[iplab], 1, 20, 2);
    tg_sig_eff[iplab]->SetMinimum(0);
    tmg_eff[iplab]->Add(tg_sig_eff[iplab],"p");
    tg_bg_eff[iplab] = new TGraph(nstep-1,x,eff_bg[iplab]);
    //set_style(tg_bg_eff[iplab], iplab==0?1:(iplab==1?2:4), 24, 1);
    set_style(tg_bg_eff[iplab], 1, 24, 2);
    tg_bg_eff[iplab]->SetMinimum(0);
    tmg_eff[iplab]->Add(tg_bg_eff[iplab],"p");

    tg_stob[iplab] = new TGraph(nstep-1,x,stob[iplab]);
    //set_style(tg_stob[iplab], iplab==0?1:(iplab==1?2:4), 20, 1);
    set_style(tg_stob[iplab], 1, 20, 3);
    tg_stob[iplab]->SetMinimum(0);

  }

  const char* bin_label[5] = {"EID(BG)","EID(SIG)","& N_{e^{+}e^{-}}=1  ","N_{#pi^{0}}>0","BtoB"};

  TCanvas *yield_step = new TCanvas("yield_step","yield_step",1800,1000);
  TPad *pad_yield[nplab];
  TH1F *hdummy_yield[nplab];
  double max_yield = 0;
  for (int iplab = 0; iplab < nplab; ++iplab) { if( yield_sig[iplab][0] > max_yield) max_yield = yield_sig[iplab][0];  }
  for (int iplab=0; iplab<nplab; ++iplab) {
    double xl = (iplab==0?0:0.1)+iplab*0.3;
    double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
    cout << "xl= " << xl << " xh= " << xh << endl;
    pad_yield[iplab] = new TPad(Form("pad_yield%d",iplab),Form("pad_yield%d",iplab),xl,0.0,xh,1.0);
    yield_step->cd(0);
    pad_yield[iplab]->Draw();
    double epsilon=1e-9;
    if (iplab==0) {pad_yield[iplab]->SetRightMargin(epsilon); pad_yield[iplab]->SetLeftMargin(0.2);}
    if (iplab==1) {pad_yield[iplab]->SetLeftMargin(epsilon); pad_yield[iplab]->SetRightMargin(epsilon); }
    if (iplab==2) {pad_yield[iplab]->SetLeftMargin(epsilon); pad_yield[iplab]->SetRightMargin(0.1);}
    pad_yield[iplab]->SetTopMargin(0.05);
    pad_yield[iplab]->SetBottomMargin(0.15);
    pad_yield[iplab]->SetTicks(0,1);
    pad_yield[iplab]->cd();
    tmg_yield[iplab]->Draw("a");
    hdummy_yield[iplab] = tmg_yield[iplab]->GetHistogram();
    for (int ibin=0; ibin<5; ++ibin){
      int iibin = 0;
      if (ibin==0) iibin = 5;
      else if (ibin==1) iibin = 32;
      else if (ibin==2) iibin = 39;
      else iibin = 36+(ibin-2)*30;
      hdummy_yield[iplab]->GetXaxis()->SetBinLabel(iibin, bin_label[ibin]);
    }
    hdummy_yield[iplab]->SetLabelSize(0.05,"Y");
    hdummy_yield[iplab]->SetLabelSize(iplab==0?0.065:0.08,"X");
    hdummy_yield[iplab]->SetLabelOffset(0.005,"Y");
    hdummy_yield[iplab]->SetTitleSize(0.06,"Y");
    hdummy_yield[iplab]->SetTitleOffset(1.7,"Y");
    hdummy_yield[iplab]->SetTitle(";;counts in 2.96<M_{e^{+}e^{-}}[GeV/c^{2}]<3.22");
    hdummy_yield[iplab]->SetMinimum(0);
    hdummy_yield[iplab]->SetMaximum(max_yield*1.1);
    tl[0][iplab]->DrawLatex(iplab==0?0.35:(iplab==1?0.25:0.15),0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[1][iplab]->DrawLatex(iplab==0?0.40:(iplab==1?0.35:0.25),0.75,Form("(%3.2f<|u|<%3.2f)",tmin[iplab],tmax[iplab]));
    if (iplab==2) leg[0]->Draw();
  }
  //yield_step->Print(Form("%s/figs/2015.09.15/yield_vs_steps.pdf",bdir));
  yield_step->Print(Form("%s/yield_vs_steps.pdf",bdir));

  TCanvas *eff_step = new TCanvas("eff_step","eff_step");
  TPad *pad_eff[nplab];
  TH1F *hdummy_eff[nplab];
  for (int iplab=0; iplab<nplab; ++iplab) {
    double xl = (iplab==0?0:0.1)+iplab*0.3;
    double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
    cout << "xl= " << xl << " xh= " << xh << endl;
    pad_eff[iplab] = new TPad(Form("pad_eff%d",iplab),Form("pad_eff%d",iplab),xl,0.0,xh,1.0);
    eff_step->cd(0);
    pad_eff[iplab]->Draw();
    double epsilon=1e-9;
    if (iplab==0) {pad_eff[iplab]->SetRightMargin(epsilon); pad_eff[iplab]->SetLeftMargin(0.2);}
    if (iplab==1) {pad_eff[iplab]->SetLeftMargin(epsilon); pad_eff[iplab]->SetRightMargin(epsilon); }
    if (iplab==2) {pad_eff[iplab]->SetLeftMargin(epsilon); pad_eff[iplab]->SetRightMargin(0.1);}
    pad_eff[iplab]->SetTicks(0,1);
    pad_eff[iplab]->cd();
    tmg_eff[iplab]->Draw("a");
    hdummy_eff[iplab] = tmg_eff[iplab]->GetHistogram();
    hdummy_eff[iplab]->SetLabelSize(0.05,"Y");
    hdummy_eff[iplab]->SetLabelSize(0.05,"X");
    hdummy_eff[iplab]->SetLabelOffset(0.005,"Y");
    hdummy_eff[iplab]->SetTitleSize(0.06,"Y");
    hdummy_eff[iplab]->SetTitleOffset(1.1,"Y");
    hdummy_eff[iplab]->SetTitle(";step number;efficiency in 2.96<M_{e^{+}e^{-}}[GeV/c^{2}]<3.22");
    hdummy_eff[iplab]->SetMinimum(0);
    hdummy_eff[iplab]->SetMaximum(105);
    tl[0][iplab]->DrawLatex(iplab==0?0.35:(iplab==1?0.25:0.15),0.23,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[1][iplab]->DrawLatex(iplab==0?0.40:(iplab==1?0.35:0.25),0.18,Form("(%3.2f<|u|<%3.2f)",tmin[iplab],tmax[iplab]));
    if (iplab==2) leg[1]->Draw();
  }
  //eff_step->Print(Form("%s/figs/2015.09.15/eff_vs_steps.pdf",bdir));

  TCanvas *stob_step = new TCanvas("stob_step","stob_step",2200,1000);
  stob_step->Divide(3,1);
  TH1F *hdummy_stob[nplab];
  for (int iplab=0; iplab<nplab; ++iplab) {
    stob_step->cd(1+iplab);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    tg_stob[iplab]->Draw("ap");
    hdummy_stob[iplab] = tg_stob[iplab]->GetHistogram();
    for (int ibin=0; ibin<5; ++ibin){
      int iibin = 0;
      if (ibin==0) iibin = 5;
      else if (ibin==1) iibin = 32;
      else if (ibin==2) iibin = 39;
      else iibin = 36+(ibin-2)*30;
      hdummy_stob[iplab]->GetXaxis()->SetBinLabel(iibin, bin_label[ibin]);
    }
    hdummy_stob[iplab]->SetLabelSize(0.05,"Y");
    hdummy_stob[iplab]->SetLabelSize(0.07,"X");
    hdummy_stob[iplab]->SetLabelOffset(0.005,"Y");
    hdummy_stob[iplab]->SetTitleSize(0.06,"Y");
    hdummy_stob[iplab]->SetTitleOffset(1.1,"Y");
    hdummy_stob[iplab]->SetTitle(";;");
    hdummy_stob[iplab]->SetMinimum(0);
    //hdummy_stob[iplab]->SetMaximum(105);
    tl[2][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    tl[3][iplab]->DrawLatex(0.35,0.75,Form("(%3.2f<|u|<%3.2f)",tmin[iplab],tmax[iplab]));
    tl[3][iplab]->DrawLatex(0.25,0.69,"Signal/Background");
    //if (iplab==0) leg[2]->Draw();
  }
  //stob_step->Print(Form("%s/figs/2015.09.15/stob_vs_steps.pdf",bdir));
  //stob_step->Print(Form("%s/stob_vs_steps.pdf",bdir));

}

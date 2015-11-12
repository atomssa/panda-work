void ana2() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  // Prepare to fit...
  const int ntu = 2;
  const char* toru[2] = {"t","u"};

  bool pol3 = true;

  //double mmin = 2.0, mmax = 5.0;
  double mmin = 2.9, mmax = 3.3; //-> 3sigma
  //double mmin = 2.8, mmax = 3.3;
  //double mmin = 2.0, mmax = 4.0;
  int immin = 0, immax = 0;

  const int ntbin_max = 20;

  TGraph *tg_sig_eff[nplab], *tg_bg_eff[nplab];
  TGraph *tg_sig_yield[nplab], *tg_bg_yield[nplab];
  TMultiGraph *tmg_yield[nplab], *tmg_eff[nplab];;

  TCanvas *tc_mep_tbins[ntu][nplab];

  TH1F* hmsg_tu[ntu][nplab][ntbin_max], *hmbg_tu[ntu][nplab][ntbin_max], *hmfg_tu[ntu][nplab][ntbin_max];
  //TH1F* hmsg_tu_nofit[ntu][nplab][ntbin_max], *hmbg_tu_nofit[ntu][nplab][ntbin_max], *hmfg_tu_nofit[ntu][nplab][ntbin_max]; // DIRTY
  TF1* fmsg_tu[ntu][nplab][ntbin_max], *fmbg_tu[ntu][nplab][ntbin_max], *fmfg_tu[ntu][nplab][ntbin_max];

  TLatex *tl[10][nplab];
  for (int ii = 0; ii < 10; ++ii) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tl[ii][iplab] = new TLatex();
      tl[ii][iplab]->SetNDC(true);
      tl[ii][iplab]->SetTextColor(ii==0?2:1);
      tl[ii][iplab]->SetTextSize((ii==0?2.5:(iplab==0?1.5:2.0))*tl[ii][iplab]->GetTextSize());
    }
  }

  int mar[2] = {20,21};
  int mar_cnt[2] = {24,25};
  int col[3] = {1,2,4};
  double tvalid[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt_er[ntu][nplab][ntbin_max]={{{0.0}}};

  double yield_bg[ntu][nplab][ntbin_max]=={{{0.0}}};
  double stob[ntu][nplab][ntbin_max]=={{{0.0}}};
  double stob_er[ntu][nplab][ntbin_max]=={{{0.0}}};

  int nptok[ntu][nplab] = {{0}};
  double chi2ndf[ntu][nplab][ntbin_max]={{{0.0}}};
  double prob[ntu][nplab][ntbin_max]={{{0.0}}};

  TGraphErrors *tg_yield[ntu][nplab];
  TGraphErrors *tg_yield_cnt[ntu][nplab];
  TGraphErrors *tg_stob[ntu][nplab];
  TGraph *tg_chi2ndf[ntu][nplab];
  TGraph *tg_prob[ntu][nplab];

  TMultiGraph *tmg_yield[ntu];
  TMultiGraph *tmg_yield_pbp[ntu][nplab];
  TMultiGraph *tmg_yield_cnt[ntu];
  TMultiGraph *tmg_yield_cnt_pbp[ntu][nplab];
  TMultiGraph *tmg_stob_pbp[ntu][nplab];
  TMultiGraph *tmg_chi2ndf[ntu];
  TMultiGraph *tmg_prob[ntu];

  for (int itu = 0; itu < ntu; ++itu) {
    tmg_yield[itu] = new TMultiGraph(Form("tmg_yield_%s",toru[itu]),Form("tmg_yield_%s",toru[itu]));
    tmg_yield_cnt[itu] = new TMultiGraph(Form("tmg_yield_cnt_%s",toru[itu]),Form("tmg_yield_cnt_%s",toru[itu]));
    tmg_chi2ndf[itu] = new TMultiGraph(Form("tmg_chi2ndf_%s",toru[itu]),Form("tmg_chi2ndf_%s",toru[itu]));
    tmg_prob[itu] = new TMultiGraph(Form("tmg_prob_%s",toru[itu]),Form("tmg_prob_%s",toru[itu]));
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tmg_yield_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_pbp%s_t%d",toru[itu],iplab),Form("tmg_yield_pbp%s_p%d",toru[itu],iplab));
      tmg_yield_cnt_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_cnt_pbp%s_p%d",toru[itu],iplab),Form("tmg_yield_cnt_pbp%s_p%d",toru[itu],plab));
      tmg_stob_pbp[itu][iplab] = new TMultiGraph(Form("tmg_stob_pbp%s_t%d",toru[itu],iplab),Form("tmg_stob_pbp%s_p%d",toru[itu],iplab));
    }
  }

  TCanvas *tctmp = new TCanvas("tctmp","tctmp");

  TLine *tline = new TLine();
  tline->SetLineColor(4);
  tline->SetLineWidth(2);
  int _ntbin[nplab] = {0};
  double t[nplab][ntbin_max] = {{0.0}}, t_er[ntbin_max] = {0.};
  double t_cnt[nplab][ntbin_max]= {{0.0}};

  TFile *fsig[nplab], *fbg[nplab];

  TLegend *legend[nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    legend[iplab] = new TLegend(0.2,0.6,0.8,0.78);
    legend[iplab]->SetFillStyle(0);
    legend[iplab]->SetBorderSize(0);
    legend[iplab]->SetTextSize(0.08);
  }

  for (int iplab=0; iplab<nplab; ++iplab) {

    fsig[iplab] = TFile::Open(Form("%s/hists/note.aug.2015/anav2_jpsi_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));
    fbg[iplab] = TFile::Open(Form("%s/hists/note.aug.2015/anav2_pip_pim_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));

    vector<double> tu_bins;
    TVectorD *vv =  (TVectorT<double>*)fsig[iplab]->Get("tu_binning");
    get_vector(tu_bins,vv);
    cout << "tu_binningsize= " << tu_bins.size() << endl;
    const int ntbin = tu_bins.size()>0?tu_bins.size()-1:12;
    _ntbin[iplab] = ntbin;

    if (tu_bins.size()==0) {
      for (int itbin = 0; itbin < ntbin; ++itbin) {
	t[iplab][itbin] = -0.45+0.1*itbin;
	t_cnt[iplab][itbin] = t[iplab][itbin]+0.01;
      }
    } else {
      for (int itbin=0; itbin < ntbin; ++itbin) {
	t[iplab][itbin] = (tu_bins[itbin+1]+tu_bins[itbin])/2.0;
	t_cnt[iplab][itbin] = t[iplab][itbin]+0.01;
      }
    }

    for (int itu = 0; itu < ntu; ++itu) {
      tc_mep_tbins[itu][iplab] = new TCanvas(Form("fitted_mep_%sbins_p%d",toru[itu],iplab),Form("fitted_mep_%sbins_p%d",toru[itu],iplab));
      tc_mep_tbins[itu][iplab]->Divide(3,3);
      for (int itbin=0; itbin<ntbin; ++itbin) {

	double tbin_width = tu_bins[itbin+1] - tu_bins[itbin];

	hmfg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_fg_p%d_%s%db",iplab,toru[itu],itbin));
	hmsg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_sig_p%d_%s%d",iplab,toru[itu],itbin));
	hmbg_tu[itu][iplab][itbin] = (TH1F*) fbg[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_bg_p%d_%s%d",iplab,toru[itu],itbin));
	if (iplab==0) hmbg_tu[itu][iplab][itbin]->Scale(81874.0/816807.0);
	if (iplab==1) hmbg_tu[itu][iplab][itbin]->Scale(224120.0/888292.0);
	if (iplab==2) hmbg_tu[itu][iplab][itbin]->Scale(10*189015.0/889395.0);
	hmfg_tu[itu][iplab][itbin]->Add(hmbg_tu[itu][iplab][itbin]);
	hmfg_tu[itu][iplab][itbin]->SetTitle(	Form("%4.2f < t < %4.2f;M_{inv}", tu_bins[itbin], tu_bins[itbin+1]));

	// Do counting before rebin, for more precise control
	immin = hmsg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmin);
	immax = hmsg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmax);
	double integral = hmsg_tu[itu][iplab][itbin]->Integral(immin, immax);
	yield_cnt[itu][iplab][itbin] = integral;
	yield_cnt_er[itu][iplab][itbin] = TMath::Sqrt(integral);
	double integral_bg = hmbg_tu[itu][iplab][itbin]->Integral(immin, immax);
	yield_bg[itu][iplab][itbin] = integral_bg;

	set_style_ana(hmfg_tu[itu][iplab][itbin], 1, 4, true);
	set_style_ana(hmsg_tu[itu][iplab][itbin], 2, 4, true);
	set_style_ana(hmbg_tu[itu][iplab][itbin], 4, 4, true);

	hmfg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.5);
	hmsg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.5);
	hmbg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.5);

	binw = hmfg_tu[itu][iplab][itbin]->GetBinWidth(3);

	tctmp->cd();
	// Define and setup fit function
	double fitmin = 2.3;
	double fitmax = 3.8;
	if (pol3) {
	  fmfg_tu[itu][iplab][itbin] = new TF1(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),fitFunctionPol3,fitmin,fitmax,7);
	  fmfg_tu[itu][iplab][itbin]->SetParameters(1,1,1,1,100,3.1,0.1);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(5,3.0,3.2);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(6,0.05,0.3);
	  fmbg_tu[itu][iplab][itbin] = new TF1(Form("fmep_bg_p%d_%s%db",iplab,toru[itu],itbin),background3,fitmin,fitmax,4);
	  fmsg_tu[itu][iplab][itbin] = new TF1(Form("fmep_sg_p%d_%s%db",iplab,toru[itu],itbin),gaussianPeak,fitmin,fitmax,3);
	} else {
	  fmfg_tu[itu][iplab][itbin] = new TF1(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),fitFunctionPol2,fitmin,fitmax,6);
	  fmfg_tu[itu][iplab][itbin]->SetParameters(1,1,1,100,3.1,0.1);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(4,3.0,3.2);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(5,0.05,0.1);
	  fmbg_tu[itu][iplab][itbin] = new TF1(Form("fmep_bg_p%d_%s%db",iplab,toru[itu],itbin),background2,fitmin,fitmax,3);
	  fmsg_tu[itu][iplab][itbin] = new TF1(Form("fmep_sg_p%d_%s%db",iplab,toru[itu],itbin),gaussianPeak,fitmin,fitmax,3);
	}
	fmfg_tu[itu][iplab][itbin]->SetNpx(500);
	fmfg_tu[itu][iplab][itbin]->SetLineWidth(2);
	fmfg_tu[itu][iplab][itbin]->SetLineColor(kCyan);
	fmbg_tu[itu][iplab][itbin]->SetLineColor(kBlack);
	fmsg_tu[itu][iplab][itbin]->SetLineColor(kBlue);

	if (integral>15) {
	  hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),"0Q");
	  fmfg_tu[itu][iplab][itbin]->SetParameter(pol3?4:3, fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3));
	  hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin), "Q+R", "ep");
	  cout << "par4 = " << fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3) << endl;
	  yield[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3);
	  yield_er[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParError(pol3?4:3);
	  tvalid[itu][iplab][nptok[itu][iplab]] = t[iplab][itbin];

	  for (int ipar=0; ipar<(pol3?4:3); ++ipar) fmbg_tu[itu][iplab][itbin]->SetParameter(ipar, fmfg_tu[itu][iplab][itbin]->GetParameter(ipar));
	  for (int ipar=(pol3?4:3); ipar<(pol3?7:6); ++ipar) fmsg_tu[itu][iplab][itbin]->SetParameter(ipar-(pol3?4:3), fmfg_tu[itu][iplab][itbin]->GetParameter(ipar));

	  stob[itu][iplab][nptok[itu][iplab]] = yield[itu][iplab][nptok[itu][iplab]]/yield_bg[itu][iplab][itbin];
	  double rel_er_num = yield_er[itu][iplab][nptok[itu][iplab]]/yield[itu][iplab][nptok[itu][iplab]];
	  double rel_er_den = sqrt(yield_bg[itu][iplab][nptok[itu][iplab]])/yield_bg[itu][iplab][nptok[itu][iplab]];
	  double rel_er = TMath::Hypot(rel_er_num, rel_er_den);
	  //stob_er[itu][iplab][nptok[itu][iplab]] = rel_er * stob[itu][iplab][nptok[itu][iplab]];
	  stob_er[itu][iplab][nptok[itu][iplab]] = 0.0; //rel_er * stob[itu][iplab][nptok[itu][iplab]];

	  yield[itu][iplab][nptok[itu][iplab]] /= tbin_width;
	  yield_er[itu][iplab][nptok[itu][iplab]] /= tbin_width;

	  nptok[itu][iplab]++;
	}

	yield_cnt[itu][iplab][itbin] /= tbin_width;
	yield_cnt_er[itu][iplab][itbin] /= tbin_width;

	if (itbin<9) { // Drawing only here. we have 3x3
	  tc_mep_tbins[itu][iplab]->cd(itbin+1);
	  hmfg_tu[itu][iplab][itbin]->Draw();
	  hmfg_tu[itu][iplab][itbin]->SetMinimum(0);
	  hmfg_tu[itu][iplab][itbin]->SetMaximum(fmfg_tu[itu][iplab][itbin]->GetMaximum()*1.1);
	  hmsg_tu[itu][iplab][itbin]->SetLineStyle(9);
	  hmbg_tu[itu][iplab][itbin]->SetLineStyle(9);
	  //hmsg_tu[itu][iplab][itbin]->Draw("hist,same");
	  hmbg_tu[itu][iplab][itbin]->Draw("hist,same");

	  fmbg_tu[itu][iplab][itbin]->Draw("same");
	  fmsg_tu[itu][iplab][itbin]->Draw("same");

	  chi2ndf[itu][iplab][itbin] = fmfg_tu[itu][iplab][itbin]->GetChisquare(); //fmfg_tu[itu][iplab][itbin]->GetNDF();
	  prob[itu][iplab][itbin] = fmfg_tu[itu][iplab][itbin]->GetProb();
	  tl[0][iplab]->DrawLatex(0.15,0.7,Form("%s",hmsg_tu[itu][iplab][itbin]->GetTitle()));

	  //tline->DrawLine(mmin,0,mmin,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());
	  //tline->DrawLine(mmax,0,mmax,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());

	}
      }

      tc_mep_tbins[itu][iplab]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_mep_tbins[itu][iplab]->GetName()));

      tg_yield[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],yield[itu][iplab],t_er,yield_er[itu][iplab]);
      tg_yield[itu][iplab]->SetMarkerStyle(mar[itu]);
      tg_yield[itu][iplab]->SetMarkerSize(1);
      tg_yield[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield[itu][iplab]->SetLineWidth(2);
      tmg_yield[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_pbp[itu][iplab]->Add(tg_yield[itu][iplab],"p");

      tg_yield_cnt[itu][iplab] = new TGraphErrors(ntbin-1,t_cnt[iplab],yield_cnt[itu][iplab],t_er,yield_cnt_er[itu][iplab]);
      tg_yield_cnt[itu][iplab]->SetMarkerStyle(mar_cnt[itu]);
      tg_yield_cnt[itu][iplab]->SetMarkerSize(1);
      tg_yield_cnt[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield_cnt[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield_cnt[itu][iplab]->SetLineWidth(2);
      tmg_yield_cnt[itu]->Add(tg_yield_cnt[itu][iplab],"p");
      tmg_yield_cnt[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_cnt_pbp[itu][iplab]->Add(tg_yield_cnt[itu][iplab],"p");
      tmg_yield_cnt_pbp[itu][iplab]->Add(tg_yield[itu][iplab],"p");
      if (itu==0) {
	legend[iplab]->AddEntry(tg_yield[itu][iplab],"Fit (fg. histo)","pl");
	legend[iplab]->AddEntry(tg_yield_cnt[itu][iplab],"Count (sig. histo)","pl");
      }

      tg_stob[itu][iplab] = new TGraphErrors(ntbin-1,tvalid[itu][iplab],stob[itu][iplab],t_er,stob_er[itu][iplab]);
      tg_stob[itu][iplab]->SetMarkerStyle(21);
      tg_stob[itu][iplab]->SetMarkerSize(1);
      tg_stob[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_stob[itu][iplab]->SetLineColor(col[iplab]);
      tg_stob[itu][iplab]->SetLineWidth(2);
      tmg_stob_pbp[itu][iplab]->Add(tg_stob[itu][iplab],"p");

      tg_chi2ndf[itu][iplab] = new TGraph(ntbin-1,t[iplab],chi2ndf[itu][iplab]);
      tg_chi2ndf[itu][iplab]->SetMarkerStyle(mar[itu]);
      tg_chi2ndf[itu][iplab]->SetMarkerSize(1);
      tg_chi2ndf[itu][iplab]->SetMarkerColor(col[iplab]);
      tmg_chi2ndf[itu]->Add(tg_chi2ndf[itu][iplab],"p");

      tg_prob[itu][iplab] = new TGraph(ntbin-1,t[iplab],prob[itu][iplab]);
      tg_prob[itu][iplab]->SetMarkerStyle(mar[itu]);
      tg_prob[itu][iplab]->SetMarkerSize(1);
      tg_prob[itu][iplab]->SetMarkerColor(col[iplab]);
      tmg_prob[itu]->Add(tg_prob[itu][iplab],"p");
    }
  }

  tctmp->Close();

  TCanvas *tc_yield_pbp[ntu];
  TCanvas *tc_yield_cnt_pbp[ntu];
  TCanvas *tc_stob_pbp[ntu];
  TCanvas *tc_chi2ndf[ntu];
  TCanvas *tc_prob[ntu];

  TPad *pad_stob[ntu][nplab];
  TH1F *hdummy_stob[ntu][nplab];

  TPad *pad_yield[ntu][nplab];
  TH1F *hdummy_yield[ntu][nplab];

  TPad *pad_yield_cnt[ntu][nplab];
  TH1F *hdummy_yield_cnt[ntu][nplab];

  double max_yield[ntu] = {0.0};
  for (int itu = 0; itu < ntu; ++itu) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double gr_max = 0;
      for (int itbin = 0; itbin < _ntbin[iplab]; ++itbin) {
	double top = yield_cnt[itu][iplab][itbin]+yield_cnt_er[itu][iplab][itbin];
	if (top > gr_max) { gr_max = top; }
      }
      if(gr_max  > max_yield[itu]) max_yield[itu] = gr_max;
    }
  }

  double max_stob[ntu] = {0.0};
  for (int itu = 0; itu < ntu; ++itu) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double gr_max = 0;
      for (int itbin = 0; itbin < _ntbin[iplab]; ++itbin) {
	//double top = stob[itu][iplab][itbin]+stob_er[itu][iplab][itbin];
	double top = stob[itu][iplab][itbin];
	if (top > gr_max) { gr_max = top; }
      }
      if(gr_max  > max_stob[itu]) max_stob[itu] = gr_max;
    }
  }


  for (int itu = 0; itu < ntu; ++itu) {

    tc_stob_pbp[itu] = new TCanvas(Form("stob_pbp%s",toru[itu]),Form("stob_pbp%s",toru[itu]));
    tc_stob_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {

      //double xl = (iplab==0?0:0.1)+iplab*0.3;
      //double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
      //pad_stob[itu][iplab] = new TPad(Form("pad_stob_%s_%d",toru[itu],iplab),Form("pad_stob_%s_%d",toru[itu],iplab),xl,0.0,xh,1.0);
      //tc_stob_pbp[itu]->cd(0);
      //pad_stob[itu][iplab]->Draw();
      //double epsilon=1e-9;
      //if (iplab==0) {pad_stob[itu][iplab]->SetRightMargin(epsilon); pad_stob[itu][iplab]->SetLeftMargin(0.2);}
      //if (iplab==1) {pad_stob[itu][iplab]->SetLeftMargin(epsilon);  pad_stob[itu][iplab]->SetRightMargin(epsilon); }
      //if (iplab==2) {pad_stob[itu][iplab]->SetLeftMargin(epsilon);  pad_stob[itu][iplab]->SetRightMargin(0.1);}
      //pad_stob[itu][iplab]->SetTicks(0,1);
      //pad_stob[itu][iplab]->cd();
      //tmg_stob_pbp[itu][iplab]->Draw("a");
      //hdummy_stob[itu][iplab] = (TH1F*)tmg_stob_pbp[itu][iplab]->GetHistogram();
      //set_style(hdummy_stob[itu][iplab]);
      //hdummy_stob[itu][iplab]->SetLabelSize(0.065,"Y");
      //hdummy_stob[itu][iplab]->SetLabelSize(0.065,"X");
      //hdummy_stob[itu][iplab]->SetLabelOffset(0.005,"Y");
      //hdummy_stob[itu][iplab]->SetTitleSize(0.06,"X");
      //hdummy_stob[itu][iplab]->SetTitleSize(0.06,"Y");
      //hdummy_stob[itu][iplab]->SetTitleOffset(1.5,"Y");
      //hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background",(itu==0?"t":"u")));
      //hdummy_stob[itu][iplab]->SetMinimum(0);
      //hdummy_stob[itu][iplab]->SetMaximum(max_stob[itu]*1.2);

      tc_stob_pbp[itu]->cd(iplab+1);
      tmg_stob_pbp[itu][iplab]->Draw("a");
      hdummy_stob[itu][iplab] = (TH1F*)tmg_stob_pbp[itu][iplab]->GetHistogram();
      set_style(hdummy_stob[itu][iplab]);
      tmg_stob_pbp[itu][iplab]->SetMinimum(0.0);
      //tmg_stob_pbp[itu][iplab]->SetMaximum(iplab==0?18:(iplab==1?200:1200));
      tmg_stob_pbp[itu][iplab]->SetMaximum(iplab==0?6.5:(iplab==1?60:350));
      hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background",(itu==0?"t":"u")));
      hdummy_stob[itu][iplab]->SetNdivisions(605);
      double epsilon = 1e-9;
      TPad *_pad = (TPad*) tc_stob_pbp[itu]->GetPad(iplab);
      _pad->SetRightMargin(epsilon);
      tl[1][iplab]->DrawLatex((iplab==0?0.2:0.15),0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

    }
    tc_stob_pbp[itu]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_stob_pbp[itu]->GetName()));

    tc_yield_pbp[itu] = new TCanvas(Form("fitted_yield_pbp%s",toru[itu]),Form("fitted_yield_pbp%s",toru[itu]));
    tc_yield_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double xl = (iplab==0?0:0.1)+iplab*0.3;
      double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
      pad_yield[itu][iplab] = new TPad(Form("pad_yield_%s_%d",toru[itu],iplab),Form("pad_yield_%s_%d",toru[itu],iplab),xl,0.0,xh,1.0);
      tc_yield_pbp[itu]->cd(0);
      pad_yield[itu][iplab]->Draw();
      double epsilon=1e-9;
      if (iplab==0) {pad_yield[itu][iplab]->SetRightMargin(epsilon); pad_yield[itu][iplab]->SetLeftMargin(0.2);}
      if (iplab==1) {pad_yield[itu][iplab]->SetLeftMargin(epsilon);  pad_yield[itu][iplab]->SetRightMargin(epsilon); }
      if (iplab==2) {pad_yield[itu][iplab]->SetLeftMargin(epsilon);  pad_yield[itu][iplab]->SetRightMargin(0.1);}
      pad_yield[itu][iplab]->SetTicks(0,1);
      pad_yield[itu][iplab]->cd();
      tmg_yield_pbp[itu][iplab]->Draw("a");
      hdummy_yield[itu][iplab] = (TH1F*)tmg_yield_pbp[itu][iplab]->GetHistogram();
      hdummy_yield[itu][iplab]->SetLabelSize(0.065,"Y");
      hdummy_yield[itu][iplab]->SetLabelSize(0.065,"X");
      hdummy_yield[itu][iplab]->SetLabelOffset(0.005,"Y");
      hdummy_yield[itu][iplab]->SetTitleSize(0.06,"X");
      hdummy_yield[itu][iplab]->SetTitleSize(0.06,"Y");
      hdummy_yield[itu][iplab]->SetTitleOffset(1.5,"Y");
      hdummy_yield[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];dN_{J/#psi}/d%s[GeV^{-2}]",(itu==0?"t":"u"),(itu==0?"t":"u")));
      hdummy_yield[itu][iplab]->SetMinimum(0);
      hdummy_yield[itu][iplab]->SetMaximum((max_yield[0]>max_yield[1]?max_yield[0]:max_yield[1])*1.4);
      tmg_yield_pbp[itu][iplab]->SetMinimum(0.0);
      tl[1][iplab]->DrawLatex((iplab==0?0.2:0.15),0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    }
    tc_yield_pbp[itu]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_yield_pbp[itu]->GetName()));

    tc_yield_cnt_pbp[itu] = new TCanvas(Form("fitted_yield_cnt_pbp%s",toru[itu]),Form("fitted_yield_cnt_pbp%s",toru[itu]));
    tc_yield_cnt_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double xl = (iplab==0?0:0.1)+iplab*0.3;
      double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
      pad_yield_cnt[itu][iplab] = new TPad(Form("pad_yield_%s_%d",toru[itu],iplab),Form("pad_yield_%s_%d",toru[itu],iplab),xl,0.0,xh,1.0);
      tc_yield_cnt_pbp[itu]->cd(0);
      pad_yield_cnt[itu][iplab]->Draw();
      double epsilon=1e-9;
      if (iplab==0) {pad_yield_cnt[itu][iplab]->SetRightMargin(epsilon); pad_yield_cnt[itu][iplab]->SetLeftMargin(0.2);}
      if (iplab==1) {pad_yield_cnt[itu][iplab]->SetLeftMargin(epsilon);  pad_yield_cnt[itu][iplab]->SetRightMargin(epsilon); }
      if (iplab==2) {pad_yield_cnt[itu][iplab]->SetLeftMargin(epsilon);  pad_yield_cnt[itu][iplab]->SetRightMargin(0.1);}
      pad_yield_cnt[itu][iplab]->SetTicks(0,1);
      pad_yield_cnt[itu][iplab]->cd();
      tmg_yield_cnt_pbp[itu][iplab]->Draw("a");
      hdummy_yield_cnt[itu][iplab] = (TH1F*)tmg_yield_cnt_pbp[itu][iplab]->GetHistogram();

      hdummy_yield_cnt[itu][iplab]->SetLabelSize(0.065,"Y");
      hdummy_yield_cnt[itu][iplab]->SetLabelSize(0.065,"X");
      hdummy_yield_cnt[itu][iplab]->SetLabelOffset(0.005,"Y");
      hdummy_yield_cnt[itu][iplab]->SetTitleSize(0.06,"X");
      hdummy_yield_cnt[itu][iplab]->SetTitleSize(0.06,"Y");
      hdummy_yield_cnt[itu][iplab]->SetTitleOffset(1.5,"Y");
      hdummy_yield_cnt[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];dN_{J/#psi}/d%s[GeV^{-2}]",(itu==0?"t":"u"),(itu==0?"t":"u")));
      hdummy_yield_cnt[itu][iplab]->SetMinimum(0);
      hdummy_yield_cnt[itu][iplab]->SetMaximum((max_yield[0]>max_yield[1]?max_yield[0]:max_yield[1])*1.4);

      tmg_yield_cnt_pbp[itu][iplab]->SetMinimum(0.0);
      tl[1][iplab]->DrawLatex(iplab==0?0.2:0.15,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

      legend[iplab]->Draw();

    }
    tc_yield_cnt_pbp[itu]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_yield_cnt_pbp[itu]->GetName()));
    //tc_yield_cnt_pbp[itu] = new TCanvas(Form("fitted_yield_cnt_pbp%s",toru[itu]),Form("fitted_yield_cnt_pbp%s",toru[itu]));
    //tc_yield_cnt_pbp[itu]->Divide(3,1);
    //for (int iplab = 0; iplab < nplab; ++iplab) {
    //  tc_yield_cnt_pbp[itu]->cd(1+iplab);
    //  tmg_yield_cnt_pbp[itu][iplab]->Draw("a");
    //  tmg_yield_cnt_pbp[itu][iplab]->SetMinimum(0.0);
    //}

    continue;
    tc_chi2ndf[itu] = new TCanvas(Form("fit_chi2ndf_%s",toru[itu]),Form("fit_chi2ndf_%s",toru[itu]));
    tc_chi2ndf[itu]->cd();
    tmg_chi2ndf[itu]->Draw("a");
    tc_prob[itu] = new TCanvas(Form("fit_prob_%s",toru[itu]),Form("fit_prob_%s",toru[itu]));
    tc_prob[itu]->cd();
    tmg_prob[itu]->Draw("a");
    tc_yield[itu]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_yield[itu]->GetName()));
  }

}

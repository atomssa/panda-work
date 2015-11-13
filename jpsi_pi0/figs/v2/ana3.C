void ana3() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  gROOT->LoadMacro(Form("%s/tda/tda.C",bdir));
  gROOT->LoadMacro(Form("%s/figs/v2/ananote.C",bdir));

  // Prepare to fit...
  const int ntu = 2;
  const char* toru[2] = {"t","u"};

  bool pol3 = false;

  //double mmin = 2.0, mmax = 5.0;
  //double mmin = 2.9, mmax = 3.3; //-> 3sigma
  double mmin = 2.8, mmax = 3.3;
  //double mmin = 2.0, mmax = 4.0;
  //double mmin = 0, mmax = 10;
  int immin = 0, immax = 0;

  const int ntbin_max = 30;

  bool msv = true;
  const double Lumi[3] = msv? {95.5,103.0,108.0} : {2e3,2e3,2e3};
  const double Br = 5.94e-2;

  //12.25 & 5.513 & -0.092 & 0.59
  //16.87 & 8.0  & -1.3 & 0.43
  //24.35 & 12.0 & -2.85 & 0.3
  //const double tvalidmin[nplab] = {-0.092, -1.3, -2.85};
  const double tvalidmin[nplab] = {-0.092, -1.0, -1.0};
  const double tvalidmax[nplab] = {0.59, 0.43, 0.3};

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
      //tl[ii][iplab]->SetTextColor(ii==0?2:1);
      if (ii==0) tl[ii][iplab]->SetTextSize(1.7 * tl[ii][iplab]->GetTextSize() );
      if (ii==1) tl[ii][iplab]->SetTextSize( (iplab==0?1.5:2.0)*tl[ii][iplab]->GetTextSize() ) ;
      if (ii==2) tl[ii][iplab]->SetTextSize(1.5 * tl[ii][iplab]->GetTextSize() );
      if (ii==3) tl[ii][iplab]->SetTextSize(1.2 * tl[ii][iplab]->GetTextSize() );
    }
  }

  int mar[2] = {20,21};
  int mar_cnt[2] = {24,25};
  int mar_cor[2] = {24,25};
  int col[3] = {1,2,4};
  double tvalid[ntu][nplab][ntbin_max]={{{0.0}}};
  double _tvalid_min[ntu][nplab][ntbin_max]={{{0.0}}};
  double _tvalid_max[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double eff_cor[ntu][nplab][ntbin_max]={{{0.0}}};
  double eff_cor_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cor[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cor_er[ntu][nplab][ntbin_max]={{{0.0}}};

  double yield_bg[ntu][nplab][ntbin_max]=={{{0.0}}};
  double yield_bg_er[ntu][nplab][ntbin_max]=={{{0.0}}};
  double stob[ntu][nplab][ntbin_max]=={{{0.0}}};
  double stob_er[ntu][nplab][ntbin_max]=={{{0.0}}};

  int nptok[ntu][nplab] = {{0}};

  TGraphErrors *tg_yield[ntu][nplab];
  TGraphErrors *tg_yield_cnt[ntu][nplab];
  TGraphErrors *tg_yield_cor[ntu][nplab];
  TGraphErrors *tg_stob[ntu][nplab];

  TMultiGraph *tmg_yield[ntu];
  TMultiGraph *tmg_yield_pbp[ntu][nplab];
  TMultiGraph *tmg_yield_cnt[ntu];
  TMultiGraph *tmg_yield_cnt_pbp[ntu][nplab];
  TMultiGraph *tmg_yield_cor[ntu];
  TMultiGraph *tmg_yield_cor_pbp[ntu][nplab];

  TMultiGraph *tmg_stob_pbp[ntu][nplab];

  for (int itu = 0; itu < ntu; ++itu) {
    tmg_yield[itu] = new TMultiGraph(Form("tmg_yield_%s",toru[itu]),Form("tmg_yield_%s",toru[itu]));
    tmg_yield_cnt[itu] = new TMultiGraph(Form("tmg_yield_cnt_%s",toru[itu]),Form("tmg_yield_cnt_%s",toru[itu]));
    tmg_yield_cor[itu] = new TMultiGraph(Form("tmg_yield_cor_%s",toru[itu]),Form("tmg_yield_cor_%s",toru[itu]));
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tmg_yield_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_pbp%s_t%d",toru[itu],iplab),Form("tmg_yield_pbp%s_p%d",toru[itu],iplab));
      tmg_yield_cnt_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_cnt_pbp%s_p%d",toru[itu],iplab),Form("tmg_yield_cnt_pbp%s_p%d",toru[itu],iplab));
      tmg_yield_cor_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_cor_pbp%s_p%d",toru[itu],iplab),Form("tmg_yield_cor_pbp%s_p%d",toru[itu],iplab));
      tmg_stob_pbp[itu][iplab] = new TMultiGraph(Form("tmg_stob_pbp%s_t%d",toru[itu],iplab),Form("tmg_stob_pbp%s_p%d",toru[itu],iplab));
    }
  }

  TCanvas *tctmp = new TCanvas("tctmp","tctmp");

  TLine *tline = new TLine();
  tline->SetLineColor(4);
  tline->SetLineWidth(2);
  int _ntbin[nplab] = {0};
  double _t[nplab][ntbin_max] = {{0.0}}, t_er[ntbin_max] = {0.};
  double _t_min[nplab][ntbin_max] = {{0.0}}, _t_max[nplab][ntbin_max] = {{0.0}};
  double t[nplab][ntbin_max] = {{0.0}}, u[nplab][ntbin_max] = {{0.0}};
  double t_cnt[nplab][ntbin_max]= {{0.0}};
  double u_cnt[nplab][ntbin_max]= {{0.0}};

  TFile *fsig[nplab], *fbg[nbg][nplab], *feff[nplab];
  TH1F *h_teff_num[nplab], *h_ueff_num[nplab];
  TH1F *h_teff_den[nplab], *h_ueff_den[nplab];

  TLegend *legend[ntu][nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    for (int itu=0; itu < ntu; ++itu) {
      legend[itu][iplab] = new TLegend(0.2,0.7,0.8,0.89);
      legend[itu][iplab]->SetFillStyle(0);
      legend[itu][iplab]->SetBorderSize(0);
      legend[itu][iplab]->SetTextSize(0.08);
    }
  }

  TLegend *legend2[ntu][nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    for (int itu=0; itu < ntu; ++itu) {
      legend2[itu][iplab] = new TLegend(0.2,0.7,0.8,0.89);
      legend2[itu][iplab]->SetFillStyle(0);
      legend2[itu][iplab]->SetBorderSize(0);
      legend2[itu][iplab]->SetTextSize(0.06);
    }
  }

  TLegend *legend3[ntu][nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    for (int itu=0; itu < ntu; ++itu) {
      legend3[itu][iplab] = new TLegend(0.4,0.25,0.99,0.35);
      legend3[itu][iplab]->SetFillStyle(0);
      legend3[itu][iplab]->SetBorderSize(0);
      legend3[itu][iplab]->SetTextSize(0.08);
    }
  }

  //const double scale[nplab][nbg] = {{1.,1.,1.,0.64}, {1.,1.,1.,0.64}, {1.,1.,1.,0.64}};
  // temporary, until fixed in effing macro
  double nevt_sim[5][3] = {{814794.0,888292.0,898721.0}, {32780.0,50142.0,51860.0}, {214780.0,174864.0,160099.0}, {570751.0,609044.0,527506.0}, {200000.0,200000.0,200000.0}};
  double nevt_xsect[5][3] = {{4.0e11, 1e11, 2e10}, {32780.0,50142.0,51860.0}, {1.15e12, 3.15e11, 6.84e10}, {3.19e12, 1.14e12, 2.92e11}, {94243.0, 157947.3, 177361.2}};
  double pi0pi0jpsi_scale[3] = {1.0};
  for (int iplab=0; iplab < nplab; ++iplab) {
    nevt_xsect[4][iplab] = nevt_xsect[1][iplab]*(nevt_xsect[2][iplab]/nevt_xsect[0][iplab]);
    pi0pi0jpsi_scale[iplab] = nevt_xsect[4][iplab]/nevt_sim[4][iplab];
    cout << "nxsect ip=" << iplab << " = " << nevt_xsect[4][iplab] << endl;
    cout << "nsim ip=" << iplab << " = " << nevt_sim[4][iplab] << endl;
    cout << "scale ip=" << iplab << " = " << pi0pi0jpsi_scale[iplab] << endl;
  }

  for (int iplab=0; iplab<nplab; ++iplab) {
    //for (int iplab=0; iplab<1; ++iplab) {

    int pass = 15;

    if (msv)
      fsig[iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_%s_msv_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    else
      fsig[iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    feff[iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0jpsi_%s4eff_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    //fbg[0][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pipm_%s_p%d_pass13.root",bdir,(ibrem==0?"raw":"brem"), iplab));
    fbg[0][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[1][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi02pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[2][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pipm2_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[3][iplab] = TFile::Open(Form("%s/hists/note.v2.oct.2015/anav2_pi0pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));

    h_teff_den[iplab] = (TH1F*) feff[iplab]->Get("tu/httrumc");
    h_ueff_den[iplab] = (TH1F*) feff[iplab]->Get("tu/hutrumc");
    h_teff_num[iplab] = (TH1F*) feff[iplab]->Get(pass<15?"tu/htrecgg_mc":"tu/htrecgg_mcut");
    h_ueff_num[iplab] = (TH1F*) feff[iplab]->Get(pass<15?"tu/hurecgg_mc":"tu/hurecgg_mcut");

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
	u_cnt[iplab][itbin] = t[iplab][itbin]+0.01;
      }
    } else {
      for (int itbin=0; itbin < ntbin; ++itbin) {
	_t[iplab][itbin] = (tu_bins[itbin+1]+tu_bins[itbin])/2.0;
	_t_min[iplab][itbin] = tu_bins[itbin];
	_t_max[iplab][itbin] = tu_bins[itbin+1];
	t[iplab][itbin] = get_mean(h_teff_den[iplab], tu_bins[itbin], tu_bins[itbin+1]);
	u[iplab][itbin] = get_mean(h_ueff_den[iplab], tu_bins[itbin], tu_bins[itbin+1]);
	if (t[iplab][itbin]<tu_bins[itbin]||t[iplab][itbin]>tu_bins[itbin+1]) {cout<<"Binning Error"<<endl; return;};
	if (u[iplab][itbin]<tu_bins[itbin]||u[iplab][itbin]>tu_bins[itbin+1]) {cout<<"Binning Error"<<endl; return;};
	t_cnt[iplab][itbin] = t[iplab][itbin]+0.01;
	u_cnt[iplab][itbin] = u[iplab][itbin]+0.01;

	//cout << "h_teff_num[iplab]= " << h_teff_num[iplab] << " h_teff_den[iplab]= " << h_teff_den[iplab] << endl;
	eff_cor[0][iplab][itbin] = integrate_content(h_teff_num[iplab], tu_bins[itbin], tu_bins[itbin+1]);
	double d0 = integrate_content(h_teff_den[iplab], tu_bins[itbin], tu_bins[itbin+1]);
	eff_cor[0][iplab][itbin] /= d0;
	double tmp0 = eff_cor[0][iplab][itbin];
	eff_cor_er[0][iplab][itbin] = TMath::Sqrt(tmp0*(1-tmp0))/TMath::Sqrt(d0);

	eff_cor[1][iplab][itbin] = integrate_content(h_ueff_num[iplab], tu_bins[itbin], tu_bins[itbin+1]);
	double d1 = integrate_content(h_ueff_den[iplab], tu_bins[itbin], tu_bins[itbin+1]);
	eff_cor[1][iplab][itbin] /= d1;
	double tmp1 = eff_cor[1][iplab][itbin];
	eff_cor_er[1][iplab][itbin] = TMath::Sqrt(tmp1*(1-tmp1))/TMath::Sqrt(d1);

      }
    }

    for (int itu = 0; itu < ntu; ++itu) {
      tc_mep_tbins[itu][iplab] = new TCanvas(Form("fitted_mep_%sbins_p%d",toru[itu],iplab),Form("fitted_mep_%sbins_p%d",toru[itu],iplab));
      tc_mep_tbins[itu][iplab]->Divide(3,2);

      cout << "\\begin{table}[hbpt]" << endl;
      cout << "  \\begin{center}" << endl;
      cout << "    \\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" << endl;

      cout << " \\hline " << endl;
      cout << "$" << toru[itu] << "_{min}$ & $" << toru[itu] << "_{max}$ & $" << toru[itu]
	   << "$ & Sig. cnt & Sig. fit & Bg. cnt & S/B & eff & dN corr & rel. err\\\\ " << endl;
      //<< " & Signal count (GeV$^-1$) & Signal+Bg fit (GeV$^-1$) & Bg count & S/B & eff & dN corr \\\\ " << endl;
      cout << " \\hline " << endl;

      for (int itbin=0; itbin<ntbin; ++itbin) {

	if (itu==0&&(tvalidmin[iplab]>t[iplab][itbin]||t[iplab][itbin]>tvalidmax[iplab])) {continue;}
	if (itu==1&&(tvalidmin[iplab]>u[iplab][itbin]||u[iplab][itbin]>tvalidmax[iplab])) {continue;}

	double tbin_width = tu_bins[itbin+1] - tu_bins[itbin];

	hmfg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_fg_p%d_%s%db",iplab,toru[itu],itbin));
	hmsg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_sig_p%d_%s%d",iplab,toru[itu],itbin));
	for (int ibg=0; ibg<nbg; ++ibg) {
	  // old
	  // hmbg_tu[itu][iplab][itbin] = (TH1F*) fbg[ibg][iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("hmep_bg_p%d_%s%d",iplab,toru[itu],itbin));
	  TH1F* htmp = (TH1F*) fbg[ibg][iplab]->Get(Form("tu_bins/hmep%s%d",toru[itu],itbin))->Clone(Form("tmp_hmep_bg_p%d_%s%d",iplab,toru[itu],itbin));
	  if (msv) htmp->Scale(0.0478);
	  if (ibg==nbg-1) htmp->Scale(pi0pi0jpsi_scale[iplab]);
	  if (ibg==0)
	    hmbg_tu[itu][iplab][itbin] = (TH1F*) htmp->Clone(Form("hmep_bg_p%d_%s%d",iplab,toru[itu],itbin));
	  else
	    hmbg_tu[itu][iplab][itbin]->Add(htmp);
	}
	//double yield_bg_scale = iplab==0?81874.0/816807.0:(iplab==1?224120.0/888292.0:10*189015.0/889395.0);
	//hmbg_tu[itu][iplab][itbin]->Scale(yield_bg_scale);

	double nsim_bg = iplab==0?816807.0:(iplab==1?888292.0:8893950.0);

	hmfg_tu[itu][iplab][itbin]->Add(hmbg_tu[itu][iplab][itbin]);
	//hmfg_tu[itu][iplab][itbin]->SetTitle(	Form("%4.2f < t < %4.2f;M_{inv}", tu_bins[itbin], tu_bins[itbin+1]));
	hmfg_tu[itu][iplab][itbin]->SetTitle(";M_{inv}");

	// Do counting before rebin, for more precise control
	immin = hmsg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmin);
	immax = hmsg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmax);
	double integral = hmsg_tu[itu][iplab][itbin]->Integral(immin, immax);
	yield_cnt[itu][iplab][nptok[itu][iplab]] = integral;
	yield_cnt_er[itu][iplab][nptok[itu][iplab]] = TMath::Sqrt(integral);
	yield_cnt[itu][iplab][nptok[itu][iplab]] /= tbin_width;
	yield_cnt_er[itu][iplab][nptok[itu][iplab]] /= tbin_width;

	double integral_bg = hmbg_tu[itu][iplab][itbin]->Integral(immin, immax);
	double integral_bg_full = hmbg_tu[itu][iplab][itbin]->Integral();
	yield_bg[itu][iplab][nptok[itu][iplab]] = integral_bg;
	double nsim_in_bin = integral_bg*nsim_bg/integral_bg_full;
	//yield_bg_er[itu][iplab][nptok[itu][iplab]] = 0.0; //TMath::Sqrt(integral_bg);
	yield_bg_er[itu][iplab][nptok[itu][iplab]] = integral_bg*sqrt(nsim_in_bin)/nsim_in_bin;
	yield_bg[itu][iplab][nptok[itu][iplab]] /= tbin_width;
	yield_bg_er[itu][iplab][nptok[itu][iplab]] /= tbin_width;
	//cout << "yieldbg= " << 	yield_bg[itu][iplab][nptok[itu][iplab]] << "yieldbger = " << 	yield_bg_er[itu][iplab][nptok[itu][iplab]] << endl;

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

	//if (integral>15) {
	//if (tvalidmin[iplab]<t[iplab][itbin]&&t[iplab][itbin]<tvalidmax[iplab]) {
	if (1) {
	  hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin),"0Q");
	  fmfg_tu[itu][iplab][itbin]->SetParameter(pol3?4:3, fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3));
	  hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%db",iplab,toru[itu],itbin), "Q+R", "ep");

	  //cout << "par4 = " << fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3) << endl;
	  yield[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3);
	  yield_er[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParError(pol3?4:3);

	  tvalid[itu][iplab][nptok[itu][iplab]] = itu==0?t[iplab][itbin]:u[iplab][itbin];
	  _tvalid_min[itu][iplab][nptok[itu][iplab]] = _t_min[iplab][itbin];
	  _tvalid_max[itu][iplab][nptok[itu][iplab]] = _t_max[iplab][itbin];
	  yield[itu][iplab][nptok[itu][iplab]] /= tbin_width;
	  yield_er[itu][iplab][nptok[itu][iplab]] /= tbin_width;

	  for (int ipar=0; ipar<(pol3?4:3); ++ipar)
	    fmbg_tu[itu][iplab][itbin]->SetParameter(ipar, fmfg_tu[itu][iplab][itbin]->GetParameter(ipar));
	  for (int ipar=(pol3?4:3); ipar<(pol3?7:6); ++ipar)
	    fmsg_tu[itu][iplab][itbin]->SetParameter(ipar-(pol3?4:3), fmfg_tu[itu][iplab][itbin]->GetParameter(ipar));

	  stob[itu][iplab][nptok[itu][iplab]] = yield[itu][iplab][nptok[itu][iplab]]/yield_bg[itu][iplab][nptok[itu][iplab]];
	  stob_er[itu][iplab][nptok[itu][iplab]] =
	    calc_err_r(yield[itu][iplab][nptok[itu][iplab]], yield_bg[itu][iplab][nptok[itu][iplab]],
		       yield_er[itu][iplab][nptok[itu][iplab]], yield_bg_er[itu][iplab][nptok[itu][iplab]]);

	  yield_cor[itu][iplab][nptok[itu][iplab]] = yield_cnt[itu][iplab][nptok[itu][iplab]]/eff_cor[itu][iplab][itbin];
	  yield_cor_er[itu][iplab][nptok[itu][iplab]] =
	    calc_err_r(yield_cnt[itu][iplab][nptok[itu][iplab]], eff_cor[itu][iplab][itbin],
		       yield_cnt_er[itu][iplab][nptok[itu][iplab]], eff_cor_er[itu][iplab][itbin]) ;
	  //cout << "t=  " << t[iplab][nptok[itu][iplab]] << " eff= " << eff_cor[itu][iplab][itbin]
	  // << " pm " << eff_cor_er[itu][iplab][itbin] << endl;

	  yield_cor[itu][iplab][nptok[itu][iplab]] /= Lumi[iplab]*Br;
	  yield_cor_er[itu][iplab][nptok[itu][iplab]] /= Lumi[iplab]*Br;

	  cout << " " << Form("%.2g", _tvalid_min[itu][iplab][nptok[itu][iplab]])
	       << " & " << Form("%.2g", _tvalid_max[itu][iplab][nptok[itu][iplab]])
	       << " & " << Form("%4.2f", tvalid[itu][iplab][nptok[itu][iplab]])
	       << " & " << Form("%.1f", yield_cnt[itu][iplab][nptok[itu][iplab]])
	       << "$\\pm$" << Form("%.1f",yield_cnt_er[itu][iplab][nptok[itu][iplab]])
	       << " &  " << Form("%.1f",yield[itu][iplab][nptok[itu][iplab]])
	       << "$\\pm$" << Form("%.1f",yield_er[itu][iplab][nptok[itu][iplab]])
	       << " &  " << Form((iplab<2?"%.1f":"%.2f"),yield_bg[itu][iplab][nptok[itu][iplab]])
	       << "$\\pm$" << Form((iplab<2?"%.1f":"%.1f"),yield_bg_er[itu][iplab][nptok[itu][iplab]])
	       << " &  " << Form("%.1f",stob[itu][iplab][nptok[itu][iplab]])
	       << "$\\pm$" << Form("%.1f",stob_er[itu][iplab][nptok[itu][iplab]])
	       << " & " << Form("%.1f",100*eff_cor[itu][iplab][nptok[itu][iplab]])
	       << "$\\pm$" << Form("%.1f",100*eff_cor_er[itu][iplab][nptok[itu][iplab]])
	       << " &  " << Form("%.1f",yield_cor[itu][iplab][nptok[itu][iplab]])
	       << "$\\pm$" << Form("%.1f",yield_cor_er[itu][iplab][nptok[itu][iplab]])
	       << " & " << Form("%.1f",100*yield_cor_er[itu][iplab][nptok[itu][iplab]]/yield_cor[itu][iplab][nptok[itu][iplab]])
	       << " \\\\ " << endl;
	  cout << " \\hline " << endl;
	}

	if (nptok[itu][iplab]<6) { // Drawing only here. we have 3x3
	  tc_mep_tbins[itu][iplab]->cd(nptok[itu][iplab]+1);
	  tc_mep_tbins[itu][iplab]->GetPad(nptok[itu][iplab]+1)->SetBottomMargin(0.15);
	  hmfg_tu[itu][iplab][itbin]->Draw();
	  hmfg_tu[itu][iplab][itbin]->GetYaxis()->SetNdivisions(505);
	  hmfg_tu[itu][iplab][itbin]->SetMinimum(0);
	  hmfg_tu[itu][iplab][itbin]->SetMaximum(fmfg_tu[itu][iplab][itbin]->GetMaximum()*1.1);
	  hmsg_tu[itu][iplab][itbin]->SetLineStyle(9);
	  hmbg_tu[itu][iplab][itbin]->SetLineStyle(9);
	  //hmsg_tu[itu][iplab][itbin]->Draw("hist,same");
	  hmbg_tu[itu][iplab][itbin]->Draw("hist,same");

	  fmbg_tu[itu][iplab][itbin]->Draw("same");
	  fmsg_tu[itu][iplab][itbin]->Draw("same");

	  //tl[0][iplab]->DrawLatex(0.15,0.7,Form("%s",hmsg_tu[itu][iplab][itbin]->GetTitle()));
	  tl[0][iplab]->DrawLatex(0.15,0.92,Form("%4.2f<t[GeV^{2}]<%4.2f", tu_bins[itbin], tu_bins[itbin+1]));

	  //tline->DrawLine(mmin,0,mmin,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());
	  //tline->DrawLine(mmax,0,mmax,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());

	}

	nptok[itu][iplab]++;

      }

      cout << "    \\end{tabular}" << endl;
      cout << "  \\end{center}" << endl;
      cout << "\\end{table}" << endl;

      //tc_mep_tbins[itu][iplab]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_mep_tbins[itu][iplab]->GetName()));

      tg_yield[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],yield[itu][iplab],t_er,yield_er[itu][iplab]);
      tg_yield[itu][iplab]->SetMarkerStyle(mar[itu]);
      tg_yield[itu][iplab]->SetMarkerSize(1);
      tg_yield[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield[itu][iplab]->SetLineWidth(2);
      tmg_yield[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_pbp[itu][iplab]->Add(tg_yield[itu][iplab],"p");

      tg_yield_cnt[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],yield_cnt[itu][iplab],t_er,yield_cnt_er[itu][iplab]);
      tg_yield_cnt[itu][iplab]->SetMarkerStyle(mar_cnt[itu]);
      tg_yield_cnt[itu][iplab]->SetMarkerSize(1);
      tg_yield_cnt[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield_cnt[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield_cnt[itu][iplab]->SetLineWidth(2);
      tmg_yield_cnt[itu]->Add(tg_yield_cnt[itu][iplab],"p");
      tmg_yield_cnt[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_cnt_pbp[itu][iplab]->Add(tg_yield_cnt[itu][iplab],"p");
      tmg_yield_cnt_pbp[itu][iplab]->Add(tg_yield[itu][iplab],"p");
      legend[itu][iplab]->AddEntry(tg_yield_cnt[itu][iplab],"Count (sig. histo)","pl");
      legend[itu][iplab]->AddEntry(tg_yield[itu][iplab],"Fit (fg. histo)","pl");

      //tg_yield_cor[itu][iplab] = new TGraphErrors(ntbin-1, itu==0?t_cnt[iplab]:u_cnt[iplab], yield_cor[itu][iplab], t_er, yield_cor_er[itu][iplab]);
      tg_yield_cor[itu][iplab] = new TGraphErrors(nptok[itu][iplab], tvalid[itu][iplab], yield_cor[itu][iplab], t_er, yield_cor_er[itu][iplab]);
      tg_yield_cor[itu][iplab]->SetMarkerStyle(mar_cor[itu]);
      tg_yield_cor[itu][iplab]->SetMarkerSize(1.5);
      tg_yield_cor[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield_cor[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield_cor[itu][iplab]->SetLineWidth(2);
      tmg_yield_cor[itu]->Add(tg_yield_cor[itu][iplab],"p");
      //tmg_yield_cor[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_cor_pbp[itu][iplab]->Add(tg_yield_cor[itu][iplab],"p");
      legend2[itu][iplab]->AddEntry(tg_yield_cor[itu][iplab],"Eff. corrected MC","pl");

      tg_stob[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],stob[itu][iplab],t_er,stob_er[itu][iplab]);
      tg_stob[itu][iplab]->SetMarkerStyle(21);
      tg_stob[itu][iplab]->SetMarkerSize(1);
      tg_stob[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_stob[itu][iplab]->SetLineColor(col[iplab]);
      tg_stob[itu][iplab]->SetLineWidth(2);
      tmg_stob_pbp[itu][iplab]->Add(tg_stob[itu][iplab],"p");
      legend3[itu][iplab]->AddEntry(tg_stob[itu][iplab],"S/B ratio","pl");

    }
  }

  tctmp->Close();

  TCanvas *tc_yield_pbp[ntu];
  TCanvas *tc_yield_cnt_pbp[ntu];
  TCanvas *tc_yield_cor_pbp[ntu];
  TCanvas *tc_stob_pbp[ntu];

  TPad *pad_stob[ntu][nplab];
  TH1F *hdummy_stob[ntu][nplab];

  TPad *pad_yield[ntu][nplab];
  TH1F *hdummy_yield[ntu][nplab];

  TPad *pad_yield_cnt[ntu][nplab];
  TH1F *hdummy_yield_cnt[ntu][nplab];

  TPad *pad_yield_cor[ntu][nplab];
  TH1F *hdummy_yield_cor[ntu][nplab];

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

  cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;

  for (int itu = 0; itu < ntu; ++itu) {

    tc_stob_pbp[itu] = new TCanvas(Form("stob_pbp%s",toru[itu]),Form("stob_pbp%s",toru[itu]),1200,500);
    //tc_stob_pbp[itu] = new TCanvas(Form("stob_pbp%s",toru[itu]),Form("stob_pbp%s",toru[itu]));
    tc_stob_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tc_stob_pbp[itu]->cd(iplab+1);
      tmg_stob_pbp[itu][iplab]->Draw("a");
      hdummy_stob[itu][iplab] = (TH1F*)tmg_stob_pbp[itu][iplab]->GetHistogram();

      hdummy_stob[itu][iplab]->SetLabelSize(iplab==2?0.045:0.065,"Y");
      hdummy_stob[itu][iplab]->SetLabelSize(0.065,"X");
      hdummy_stob[itu][iplab]->SetLabelOffset(0.005,"Y");
      hdummy_stob[itu][iplab]->SetTitleSize(0.06,"X");
      hdummy_stob[itu][iplab]->SetTitleSize(iplab==2?0.045:(iplab==1?0.058:0.06),"Y");
      hdummy_stob[itu][iplab]->SetTitleOffset(iplab==2?1.55:(iplab==1?1.1:1.0),"Y");
      if (iplab==0)
	hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background          ",(itu==0?"t":"u")));
      else if (iplab==1)
	hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background          ",(itu==0?"t":"u")));
      else
	hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background          ",(itu==0?"t":"u")));
      hdummy_stob[itu][iplab]->SetNdivisions(605);
      //hdummy_stob[itu][iplab]->GetXaxis()->SetRangeUser(-3,3);
      //set_style(hdummy_stob[itu][iplab]);

      tmg_stob_pbp[itu][iplab]->SetMinimum(0.0);
      //tmg_stob_pbp[itu][iplab]->SetMaximum(iplab==0?20:(iplab==1?90:1300));
      //tmg_stob_pbp[itu][iplab]->SetMaximum(iplab==0?18:(iplab==1?200:1200));
      //tmg_stob_pbp[itu][iplab]->SetMaximum(iplab==0?6.5:(iplab==1?60:350));

      double epsilon = 1e-9;
      TPad *_pad = (TPad*) tc_stob_pbp[itu]->GetPad(iplab);
      _pad->SetRightMargin(0.03);
      double _x_st = (iplab==0?0.2:(iplab==1?0.2:0.2));
      tl[2][iplab]->DrawLatex(_x_st,0.95,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
      tl[3][iplab]->DrawLatex(_x_st+0.34,0.23,"S: #bar{p}p#rightarrow#pi^{0}J/#psi");
      tl[3][iplab]->DrawLatex(_x_st+0.34,0.17,"B: #bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-}");
      legend3[itu][iplab]->Draw();
    }
    //tc_stob_pbp[itu]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_stob_pbp[itu]->GetName()));

    /*
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
    */

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
      hdummy_yield_cnt[itu][iplab]->SetMaximum((max_yield[0]>max_yield[1]?max_yield[0]:max_yield[1])*1.2);

      tmg_yield_cnt_pbp[itu][iplab]->SetMinimum(0.0);
      tl[1][iplab]->DrawLatex(iplab==0?0.33:0.15,0.93,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

      legend[itu][iplab]->Draw();

    }
    //tc_yield_cnt_pbp[itu]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_yield_cnt_pbp[itu]->GetName()));

    tc_yield_cor_pbp[itu] = new TCanvas(Form("fitted_yield_cor_pbp%s",toru[itu]),Form("fitted_yield_cor_pbp%s",toru[itu]), 1400, 550);
    tc_yield_cor_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tc_yield_cor_pbp[itu]->cd(iplab+1);

      double epsilon=1e-9;
      pad_yield_cor[itu][iplab] = (TPad*) tc_yield_cor_pbp[itu]->GetPad(iplab+1);
      pad_yield_cor[itu][iplab]->SetRightMargin(0.03);
      pad_yield_cor[itu][iplab]->SetLeftMargin(0.2);

      tmg_yield_cor_pbp[itu][iplab]->Draw("a");

      hdummy_yield_cor[itu][iplab] = (TH1F*)tmg_yield_cor_pbp[itu][iplab]->GetHistogram();
      hdummy_yield_cor[itu][iplab]->SetLabelSize(0.065,"Y");
      hdummy_yield_cor[itu][iplab]->SetLabelSize(0.065,"X");
      hdummy_yield_cor[itu][iplab]->SetLabelOffset(0.005,"Y");
      hdummy_yield_cor[itu][iplab]->SetTitleSize(0.06,"X");
      hdummy_yield_cor[itu][iplab]->SetTitleSize(0.06,"Y");
      hdummy_yield_cor[itu][iplab]->SetTitleOffset(1.5,"Y");
      hdummy_yield_cor[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];d#sigma_{J/#psi-#pi^{0}}/d%s[pb/GeV^{2}]",(itu==0?"t":"u"),(itu==0?"t":"u")));
      if (iplab==1) {
	hdummy_yield_cor[itu][iplab]->GetXaxis()->SetNdivisions(508);
      }
      //hdummy_yield_cor[itu][iplab]->SetMinimum(0);
      //hdummy_yield_cor[itu][iplab]->SetMaximum((max_yield[0]>max_yield[1]?max_yield[0]:max_yield[1])*1.4);

      tmg_yield_cor_pbp[itu][iplab]->SetMinimum(0.0);

      TF1 *funkyfunk = get_func(iplab, tvalid[itu][iplab][0], tvalid[itu][iplab][nptok[itu][iplab]-1] );
      funkyfunk->SetLineColor(col[iplab]);
      funkyfunk->Draw("same");
      legend2[itu][iplab]->AddEntry(funkyfunk,"TDA model","pl");
      legend2[itu][iplab]->Draw();

      tl[2][iplab]->DrawLatex(0.33,0.93,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

    }
    //tc_yield_cor_pbp[itu]->Print(Form("%s/figs/2015.09.15/%s.pdf",bdir,tc_yield_cor_pbp[itu]->GetName()));

  }

}

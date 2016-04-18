#include "ananote.C"

void ana4(int icth = 8) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  bool verbose = false;
  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  //gROOT->LoadMacro(Form("%s/figs/pv0/ananote.C",bdir));
  //const int nplab =3;

  // Prepare to fit...
  const int ntu = 4;// until indexing bug with u is fixed
  const char* toru[2] = {"t","u"};

  bool pol3 = true;

  //double mmin = 2.0, mmax = 5.0;
  //double mmin = 2.9, mmax = 3.3; //-> 3sigma
  double mmin = 2.8, mmax = 3.3;
  //double mmin = 2.0, mmax = 4.0;
  //double mmin = 0, mmax = 10;
  int immin = 0, immax = 0;

  const int ntbin_max = 8;

  bool msv = false;
  const double Lumi_msv[3] = {95.5,103.0,108.0};
  const double Lumi_full[3] = {2.e3,2.e3,2.e3};
  //const double Lumi[3] = {0.};

  const double Br = 5.94e-2;

  //12.25 & 5.513 & -0.092 & 0.59
  //16.87 & 8.0  & -1.3 & 0.43
  //24.35 & 12.0 & -2.85 & 0.3
  //const double tvalidmin[nplab] = {-0.092, -1.3, -2.85};
  const double tvalidmin[nplab] = {-0.092, -1.0, -1.0};
  const double tvalidmax[nplab] = {0.59, 0.43, 0.3};

  TGraph *tg_sig_eff[nplab], *tg_bg_eff[nplab];
  TGraph *tg_sig_yield[nplab], *tg_bg_yield[nplab];

  TCanvas *tc_mep_tbins[ntu][nplab];

  TH1F* hmsg_tu[ntu][nplab][ntbin_max], *hmbg_tu[ntu][nplab][ntbin_max], *hmfg_tu[ntu][nplab][ntbin_max];
  //TH1F* hmsg_tu_nofit[ntu][nplab][ntbin_max], *hmbg_tu_nofit[ntu][nplab][ntbin_max], *hmfg_tu_nofit[ntu][nplab][ntbin_max]; // DIRTY
  TF1* fmsg_tu[ntu][nplab][ntbin_max], *fmbg_tu[ntu][nplab][ntbin_max], *fmfg_tu[ntu][nplab][ntbin_max];

  TF1* fcosth_tu[ntu][nplab];

  TLatex *tl[10][nplab];
  for (int ii = 0; ii < 10; ++ii) {
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tl[ii][iplab] = new TLatex();
      tl[ii][iplab]->SetNDC(true);
      //tl[ii][iplab]->SetTextColor(ii==0?2:1);
      if (ii==0) tl[ii][iplab]->SetTextSize(1.7 * tl[ii][iplab]->GetTextSize() );
      if (ii==1) tl[ii][iplab]->SetTextSize( (iplab==0?1.5:2.0)*tl[ii][iplab]->GetTextSize() ) ;
      if (ii==2) tl[ii][iplab]->SetTextSize(1.3 * tl[ii][iplab]->GetTextSize() );
      if (ii==3) tl[ii][iplab]->SetTextSize(1.3 * tl[ii][iplab]->GetTextSize() );
      if (ii==5) tl[ii][iplab]->SetTextSize(1.2 * tl[ii][iplab]->GetTextSize() );
    }
  }

  int mar[ntu] = {20,20,21,21};
  int mar_cnt[ntu] = {24,24,25,25};
  int mar_cor[ntu] = {24,24,25,25};

  int col[3] = {1,2,4};

  double tvalid[ntu][nplab][ntbin_max]={{{0.0}}};
  double _tvalid_min[ntu][nplab][ntbin_max]={{{0.0}}};
  double _tvalid_max[ntu][nplab][ntbin_max]={{{0.0}}};

  double yield[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt_fg[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt_fg_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cnt_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double eff_cor[ntu][nplab][ntbin_max]={{{0.0}}};
  double eff_cor_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cor[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_cor_er[ntu][nplab][ntbin_max]={{{0.0}}};

  double s_yield_bg[ntu][nplab][ntbin_max]={{{0.0}}};
  double s_yield_bg_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg_cor[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg_cor_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg0[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg0_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg0_cor[ntu][nplab][ntbin_max]={{{0.0}}};
  double yield_bg0_cor_er[ntu][nplab][ntbin_max]={{{0.0}}};
  double d_yield_bg_min[ntu][nplab]={{0.}},d_yield_bg_max[ntu][nplab]={{0.}};

  double stob[ntu][nplab][ntbin_max]={{{0.0}}};
  double stob_er[ntu][nplab][ntbin_max]={{{0.0}}};

  int nptok[ntu][nplab] = {{0}};

  TGraphErrors *tg_yield[ntu][nplab];
  TGraphErrors *tg_yield_cnt[ntu][nplab];
  TGraphErrors *tg_yield_cor[ntu][nplab];
  TGraphErrors *tg_stob[ntu][nplab];
  TGraphErrors *tg_yield_bg[ntu][nplab];
  TGraphErrors *tg_yield_bg_cor[ntu][nplab];
  TGraphErrors *tg_yield_bg0_cor[ntu][nplab];
  TH1F *h_yield_bg[ntu][nplab];
  TH1F *h_yield_bg_cor[ntu][nplab];
  TH1F *h_yield_bg0_cor[ntu][nplab];
  TH1F* hcostfit[ntu][nplab];

  TMultiGraph *tmg_yield[ntu];
  TMultiGraph *tmg_yield_pbp[ntu][nplab];
  TMultiGraph *tmg_yield_cnt[ntu];
  TMultiGraph *tmg_yield_cnt_pbp[ntu][nplab];
  TMultiGraph *tmg_yield_cor[ntu];
  TMultiGraph *tmg_yield_cor_pbp[ntu][nplab];

  TMultiGraph *tmg_stob_pbp[ntu][nplab];

  for (int itu = 0; itu < ntu; ++itu) {
    tmg_yield[itu] = new TMultiGraph(Form("tmg_yield_%s%d",toru[itu/2],itu%2),Form("tmg_yield_%s%d",toru[itu/2],itu%2));
    tmg_yield_cnt[itu] = new TMultiGraph(Form("tmg_yield_cnt_%s%d",toru[itu/2],itu%2),Form("tmg_yield_cnt_%s%d",toru[itu/2],itu%2));
    tmg_yield_cor[itu] = new TMultiGraph(Form("tmg_yield_cor_%s%d",toru[itu/2],itu%2),Form("tmg_yield_cor_%s%d",toru[itu/2],itu%2));
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tmg_yield_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_pbp_%s%d_p%d",toru[itu/2],itu%2,iplab),Form("tmg_yield_pbp_%s%d_p%d",toru[itu/2],itu%2,iplab));
      tmg_yield_cnt_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_cnt_pbp_%s%d_p%d",toru[itu/2],itu%2,iplab),Form("tmg_yield_cnt_pbp_%s%d_p%d",toru[itu/2],itu%2,iplab));
      tmg_yield_cor_pbp[itu][iplab] = new TMultiGraph(Form("tmg_yield_cor_pbp%s%d_p%d",toru[itu/2],itu%2,iplab),Form("tmg_yield_cor_pbp_%s%d_p%d",toru[itu/2],itu%2,iplab));
      tmg_stob_pbp[itu][iplab] = new TMultiGraph(Form("tmg_stob_pbp_%s%d_t%d",toru[itu/2],itu%2,iplab),"");//,Form("tmg_stob_pbp_%s%d_p%d",toru[itu/2],itu%2,iplab));
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
  TH1F *h_eff_den[ntu][nplab], *h_eff_num[ntu][nplab];

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
      legend2[itu][iplab] = new TLegend(0.325,0.23,0.8,0.37);
      legend2[itu][iplab]->SetFillStyle(0);
      legend2[itu][iplab]->SetBorderSize(0);
      legend2[itu][iplab]->SetTextSize(0.06);
    }
  }

  TLegend *legend3[ntu][nplab];
  for (int iplab=0; iplab < nplab; ++iplab) {
    for (int itu=0; itu < ntu; ++itu) {
      //legend3[itu][iplab] = new TLegend(0.4,0.25,0.99,0.35);
      legend3[itu][iplab] = new TLegend(0.13,0.15,0.69,0.25);
      legend3[itu][iplab]->SetFillStyle(0);
      legend3[itu][iplab]->SetBorderSize(0);
      legend3[itu][iplab]->SetTextSize(0.08);
    }
  }

  double nevt_sim_sg[3] = {32780.0,50142.0,51860.0};
  double nevt_xsect_sg[3] = {32780.0,50142.0,51860.0};
  double nevt_sim_bg[4][3] = {{814794.0,888292.0,898721.0}, {214780.0,174864.0,160099.0}, {570751.0,609044.0,527506.0}, {1.98e6,2.0e6,1.97e7}};
  double nevt_xsect_bg[4][3] = {{4.0e11, 1e11, 2e10}, {1.15e12, 3.15e11, 6.84e10}, {3.19e12, 1.14e12, 2.92e11}, {94243.0, 157947.3, 177361.2}};

  // Undo event scaling applied in anav2 module, to do the scaling anew with proper number of events
  //double pi0pi0jpsi_unscale_anav2[3] = {200000.0/94243.0,200000.0/157947.3,200000.0/177361.2};
  double pi0pi0jpsi_unscale_anav2[3] = {1.0,1.0,1.0};
  double pi0pi0jpsi_re_scale[3] = {1.0};
  for (int iplab=0; iplab < nplab; ++iplab) {
    //nevt_xsect_bg[3][iplab] = nevt_xsect_sg[iplab]*(nevt_xsect_bg[1][iplab]/nevt_xsect_bg[0][iplab]);
    pi0pi0jpsi_re_scale[iplab] = pi0pi0jpsi_unscale_anav2[iplab]*nevt_xsect_bg[3][iplab]/nevt_sim_bg[3][iplab];
    //cout << "nxsect ip=" << iplab << " = " << nevt_xsect_bg[3][iplab] << endl;
    //cout << "nsim ip=" << iplab << " = " << nevt_sim_bg[3][iplab] << endl;
    //cout << "scale ip=" << iplab << " = " << pi0pi0jpsi_re_scale[iplab] << endl;
  }

  ofstream dat_out;
  bool write_dat = false;
  int iwt = 0;
  bool kin_fitted = true;

  for (int iplab=0; iplab<nplab; ++iplab) {
  //for (int iplab=0; iplab<1; ++iplab) {

    int pass = 18;

    int ncth = (iplab==0)?198:(iplab==1?183:187);
    //int ncth = (iplab==0)?99:(iplab==1?92:94);

    if (icth>=ncth) continue;

    cout << "icth=" << icth << " iplab= " << iplab << endl;

    if (write_dat)
      dat_out.open(Form("%s/%s%s/data_table_p%d_cth%d.dat",
			(msv?"msv":"full"),
			(iwt==0||iwt==1?Form("wt%d",iwt):"nowt"),
			(kin_fitted?"_kinfit":"_nokinfit"),
			iplab,icth));

    if (msv)
      fsig[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_%s_msv_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    else
      //fsig[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/cth.2x/anav2_pi0jpsi_%s4cth.2x.%d_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), icth, iplab, pass));
      fsig[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/cth/anav2_pi0jpsi_%s4cth.%d_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), icth, iplab, pass));

    feff[iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_%s4eff_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[0][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[1][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pi0jpsi_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[2][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi02pipm_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));
    fbg[3][iplab] = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pipm2_%s_p%d_pass%d.root",bdir,(ibrem==0?"raw":"brem"), iplab, pass));

    //cout << "fsig[iplab] = " << fsig[iplab] << endl;

    // this has to be fixed (differentiate for each t bin)
    for (int itu=0; itu < ntu; ++itu) {
      //cout << Form("epeff/hepcosth_jpsi_mc_itu%d",itu) << endl;
      //cout << Form("epeff/%shepcosth_jpsi_rec_itu%d",(kin_fitted?"f_":""),itu) << endl;
      h_eff_den[itu][iplab] = (TH1F*) feff[iplab]->Get(Form("epeff/hepcosth_jpsi_mc%s_itu%d",(iwt==0||iwt==1?Form("_wt%d",iwt):""),itu));
      h_eff_num[itu][iplab] = (TH1F*) feff[iplab]->Get(Form("epeff/%shepcosth_jpsi_rec%s_itu%d",(kin_fitted?"f_":""),(iwt==0||iwt==1?Form("_wt%d",iwt):""),itu));
    }

    //cout << "Getting efficincies" << endl;
    vector<double> tu_bins;
    const int ntbin = ntbin_max;
    for (int itbin = 0; itbin < ntbin; ++itbin) {
      //double d_tmin = -1.0 + (2.0*itbin/8.0);
      //double d_tmax = -1.0 + (2.0*(itbin+1.0)/8.0);
      // temporary!!
      //double d_tmin = (1.0*itbin/6.0);
      //double d_tmax = (1.0*(itbin+1.0)/6.0);
      double d_tmin = -0.8 + (2.0*itbin/10.0);
      double d_tmax = -0.8 + (2.0*(itbin+1.0)/10.0);

      if (itbin==0) {
	tu_bins.push_back(d_tmin);
      }
      tu_bins.push_back(d_tmax);

      t[iplab][itbin] = (d_tmin+d_tmax)/2.0;
      u[iplab][itbin] = (d_tmin+d_tmax)/2.0;

      t_cnt[iplab][itbin] = t[iplab][itbin]+0.01;
      u_cnt[iplab][itbin] = t[iplab][itbin]+0.01;

      for (int itu=0; itu < ntu; ++itu) {
	eff_cor[itu][iplab][itbin] = integrate_content(h_eff_num[itu][iplab], d_tmin, d_tmax);
	double d0 = integrate_content(h_eff_den[itu][iplab], d_tmin, d_tmax);
	if (verbose)
	  cout << "tmin= " << d_tmin << " t_max= " << d_tmax << " EffNum = " << eff_cor[itu][iplab][itbin]
	       << " EffDen = " << d0 << " Eff= " << 	eff_cor[itu][iplab][itbin] / d0 << endl;
	eff_cor[itu][iplab][itbin] /= d0;
	double tmp0 = eff_cor[itu][iplab][itbin];
	eff_cor_er[itu][iplab][itbin] = TMath::Sqrt(tmp0*(1-tmp0))/TMath::Sqrt(d0);
      }
    }
    _ntbin[iplab] = ntbin;

    double bbb[30] = {0.0};
    for (unsigned int ibin=0; ibin < tu_bins.size(); ++ibin) { bbb[ibin] = tu_bins[ibin]; }
    h_yield_bg[1][iplab] = new TH1F(Form("h_yield_bg_p%d_1",iplab),Form("h_yield_bg_p%d_1",iplab),ntbin-1,bbb);
    h_yield_bg[3][iplab] = new TH1F(Form("h_yield_bg_p%d_3",iplab),Form("h_yield_bg_p%d_3",iplab),ntbin-1,bbb);
    h_yield_bg_cor[1][iplab] = new TH1F(Form("h_yield_bg_cor_p%d_1",iplab),Form("h_yield_bg_cor_p%d_1",iplab),ntbin-1,bbb);
    h_yield_bg_cor[3][iplab] = new TH1F(Form("h_yield_bg_cor_p%d_3",iplab),Form("h_yield_bg_cor_p%d_3",iplab),ntbin-1,bbb);
    h_yield_bg0_cor[1][iplab] = new TH1F(Form("h_yield_bg_cor0_p%d_1",iplab),Form("h_yield_bg_cor0_p%d_1",iplab),ntbin-1,bbb);
    h_yield_bg0_cor[3][iplab] = new TH1F(Form("h_yield_bg_cor0_p%d_3",iplab),Form("h_yield_bg_cor0_p%d_3",iplab),ntbin-1,bbb);

    for (int itu = 0; itu < ntu; ++itu) {

      if (itu==0||itu==2) continue;

      tc_mep_tbins[itu][iplab] = new TCanvas(Form("fitted_mep_%s%d_p%d",toru[itu/2],itu%2,iplab),Form("fitted_mep_%s%d_p%d",toru[itu/2],itu%2,iplab));
      tc_mep_tbins[itu][iplab]->Divide(3,3);

      for (int itbin=0; itbin<ntbin; ++itbin) {

	if (itbin >= ntbin_max) break;

	// maybe use content != 0 or relative error minimum for condition to add here
	//if (itu==0&&(tvalidmin[iplab]>t[iplab][itbin]||t[iplab][itbin]>tvalidmax[iplab])) {continue;}
	//if (itu==1&&(tvalidmin[iplab]>u[iplab][itbin]||u[iplab][itbin]>tvalidmax[iplab])) {continue;}

	//double tbin_width = tu_bins[itbin+1] - tu_bins[itbin];
	double tbin_width = 1.0;

	const char *name_in = Form("%scosthbins/%shmep_%s%d_cth%d%s",toru[itu/2],(kin_fitted?"f_":""),toru[itu/2],itu%2,itbin,(iwt==0||iwt==1?Form("_wt%d",iwt):""));
	//cout << "name_in= " << name_in << endl;

	const char *name_out = Form("hmep_fg_p%d_%s%d_cth%d",iplab,toru[itu/2],itu%2,itbin);
	hmfg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(name_in)->Clone(name_out);
	name_out = Form("hmep_sg_p%d_%s%d_cth%d",iplab,toru[itu/2],itu%2,itbin);
	hmsg_tu[itu][iplab][itbin] = (TH1F*) fsig[iplab]->Get(name_in)->Clone(name_out);

	immin = hmsg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmin);
	immax = hmsg_tu[itu][iplab][itbin]->GetXaxis()->FindBin(mmax);

	//cout << "sig title= " << hmfg_tu[itu][iplab][itbin]->GetTitle() << endl;
	//cout << "tmean = " << t[iplab][itbin] << endl;

	//if (hmsg_tu[itu][iplab][itbin]->Integral()<2) continue;

	if (eff_cor[itu][iplab][itbin]<0.01||(eff_cor[itu][iplab][itbin]!=eff_cor[itu][iplab][itbin])) continue;

	//cout << "========= iplab= " << iplab << "itbin= " << itbin << " ===========" << endl;
	for (int ibg=0; ibg<nbg; ++ibg) {
	  const char* ftmp_in = Form("%scosthbins/hmep_%s%d_cth%d",toru[itu/2],toru[itu/2],itu%2,itbin);
	  const char* ftmp_out = Form("tmp_hmep_bg_p%d_%s%d_cth%d",iplab,toru[itu/2],itu%2,itbin);
	  TH1F* htmp = (TH1F*) fbg[ibg][iplab]->Get(ftmp_in)->Clone(ftmp_out);

	  if (msv) htmp->Scale(0.0478);
	  if (ibg==1) htmp->Scale(pi0pi0jpsi_re_scale[iplab]*50);

	  if (ibg==0)
	    hmbg_tu[itu][iplab][itbin] = (TH1F*) htmp->Clone(Form("hmep_bg_p%d_%s%d_cth%d",iplab,toru[itu/2],itu%2,itbin));
	  else
	    hmbg_tu[itu][iplab][itbin]->Add(htmp);

	  double integral_bg = htmp->Integral(immin, immax);
	  double integral_bg_full = htmp->Integral();
	  double entries_bg_full = htmp->GetEntries();
	  if (integral_bg > 0
	      //&& ibg!=1 && ibg!=2
	      //&& !(ibg==3&&itu==0)
	      ) {
	    //double nsim_in_bin = nevt_sim_bg[ibg][iplab]*(integral_bg/integral_bg_full);
	    double nrec_in_bin = entries_bg_full*(integral_bg/integral_bg_full);
	    double _scale = entries_bg_full/integral_bg_full;
	    //double yield_bg_er_indiv = _scale*sqrt(nrec_in_bin);
	    double yield_bg_er_indiv = integral_bg* sqrt(nrec_in_bin)/nrec_in_bin;
	    //double yield_bg_er_indiv = integral_bg*sqrt(nsim_in_bin)/nsim_in_bin;
	    //cout << "ibg= " << ibg << " y= " <<  integral_bg
	    //	 << " ent= " << entries_bg_full << " int= " << integral_bg_full << " scale = " << _scale
	    //	 << " nrec= " << nrec_in_bin <<  " nrec_er= " << sqrt(nrec_in_bin)
	    //	 << " yieldbger = " << yield_bg_er_indiv << endl;
	    yield_bg[itu][iplab][nptok[itu][iplab]] += integral_bg;
	    double yield_bg_er_tot = TMath::Hypot(yield_bg_er[itu][iplab][nptok[itu][iplab]], yield_bg_er_indiv);
	    yield_bg_er[itu][iplab][nptok[itu][iplab]] = yield_bg_er_tot;

	    if (ibg==0) {
	      yield_bg0[itu][iplab][nptok[itu][iplab]] += integral_bg;
	      yield_bg0_er[itu][iplab][nptok[itu][iplab]] = yield_bg_er_indiv;
	    }

	  }
	}

	hmfg_tu[itu][iplab][itbin]->Add(hmbg_tu[itu][iplab][itbin]);
	//hmfg_tu[itu][iplab][itbin]->SetTitle(	Form("%4.2f < t < %4.2f;M_{inv}", tu_bins[itbin], tu_bins[itbin+1]));
	hmfg_tu[itu][iplab][itbin]->SetTitle(";M_{inv}");

	// Do counting before rebin, for more precise control
	double integral = hmsg_tu[itu][iplab][itbin]->Integral(immin, immax);
	yield_cnt[itu][iplab][nptok[itu][iplab]] = integral;
	yield_cnt_er[itu][iplab][nptok[itu][iplab]] = TMath::Sqrt(integral);
	yield_cnt[itu][iplab][nptok[itu][iplab]] /= tbin_width;
	yield_cnt_er[itu][iplab][nptok[itu][iplab]] /= tbin_width;

	double integral_fg = hmfg_tu[itu][iplab][itbin]->Integral(immin, immax);
	yield_cnt_fg[itu][iplab][nptok[itu][iplab]] = integral_fg;
	yield_cnt_fg_er[itu][iplab][nptok[itu][iplab]] = TMath::Sqrt(integral_fg);

	set_style_ana(hmfg_tu[itu][iplab][itbin], 1, 2, false);
	set_style_ana(hmsg_tu[itu][iplab][itbin], 2, 2, false);
	set_style_ana(hmbg_tu[itu][iplab][itbin], 4, 2, false);

	hmfg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.5);
	hmsg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.5);
	hmbg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(1.3,4.5);

	binw = hmfg_tu[itu][iplab][itbin]->GetBinWidth(3);

	tctmp->cd();
	// Define and setup fit function
	double fitmin = 2.3;
	double fitmax = 3.8;
	if (pol3) {
	  fmfg_tu[itu][iplab][itbin] = new TF1(Form("fmep_fg_p%d_%s%d_%d",iplab,toru[itu/2],itu%2,itbin),fitFunctionPol3,fitmin,fitmax,7);
	  fmfg_tu[itu][iplab][itbin]->SetParameters(1,1,1,1,0,3.1,0.1);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(5,3.0,3.2);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(6,0.05,0.3);
	  fmbg_tu[itu][iplab][itbin] = new TF1(Form("fmep_bg_p%d_%s%d_%d",iplab,toru[itu/2],itu%2,itbin),background3,fitmin,fitmax,4);
	  fmsg_tu[itu][iplab][itbin] = new TF1(Form("fmep_sg_p%d_%s%d_%d",iplab,toru[itu/2],itu%2,itbin),gaussianPeak,fitmin,fitmax,3);
	} else {
	  fmfg_tu[itu][iplab][itbin] = new TF1(Form("fmep_fg_p%d_%s%d_%d",iplab,toru[itu/2],itu%2,itbin),fitFunctionPol2,fitmin,fitmax,6);
	  fmfg_tu[itu][iplab][itbin]->SetParameters(1,1,1,0,3.1,0.1);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(4,3.0,3.2);
	  fmfg_tu[itu][iplab][itbin]->SetParLimits(5,0.05,0.1);
	  fmbg_tu[itu][iplab][itbin] = new TF1(Form("fmep_bg_p%d_%s%d_%d", iplab,toru[itu/2],itu%2,itbin),background2,fitmin,fitmax,3);
	  fmsg_tu[itu][iplab][itbin] = new TF1(Form("fmep_sg_p%d_%s%d_%db",iplab,toru[itu/2],itu%2,itbin),gaussianPeak,fitmin,fitmax,3);
	}
	fmfg_tu[itu][iplab][itbin]->SetNpx(500);
	fmfg_tu[itu][iplab][itbin]->SetLineWidth(2);
	fmfg_tu[itu][iplab][itbin]->SetLineColor(kBlack);
	fmbg_tu[itu][iplab][itbin]->SetLineColor(kBlue);
	fmsg_tu[itu][iplab][itbin]->SetLineColor(kRed);

	//if (integral>15) {
	//if (tvalidmin[iplab]<t[iplab][itbin]&&t[iplab][itbin]<tvalidmax[iplab]) {
	//if (hmfg_tu[itu][iplab][itbin]->Integral(immin,immax)>10) {
	if (true) {


	  hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%d_%d",iplab,toru[itu/2],itu%2,itbin),"0Q");
	  fmfg_tu[itu][iplab][itbin]->SetParameter(pol3?4:3, fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3));
	  TFitResultPtr fit_res = hmfg_tu[itu][iplab][itbin]->Fit(Form("fmep_fg_p%d_%s%d_%d",iplab,toru[itu/2],itu%2,itbin), "Q+RS", "ep");
	  for (int ipar=0; ipar<(pol3?4:3); ++ipar)
	    fmbg_tu[itu][iplab][itbin]->SetParameter(ipar, fmfg_tu[itu][iplab][itbin]->GetParameter(ipar));
	  for (int ipar=(pol3?4:3); ipar<(pol3?7:6); ++ipar)
	    fmsg_tu[itu][iplab][itbin]->SetParameter(ipar-(pol3?4:3), fmfg_tu[itu][iplab][itbin]->GetParameter(ipar));

	  double params2[3], params3[4], covmat2[9], covmat3[16];
	  get_covmat(pol3?4:3, pol3?7:6, fmfg_tu[itu][iplab][itbin], fit_res, pol3?params3:params2, pol3?covmat3:covmat2);

	  if (verbose) {
	    if (itu==1&&itbin==0&&iplab==0) {
	      cout << "fit result full:" << endl;
	      print_pars(pol3?7:6, fmfg_tu[itu][iplab][itbin]);
	      cout << "fit result bgonly:" << endl;
	      print_pars(pol3?4:3, fmbg_tu[itu][iplab][itbin]);
	      cout << "fit result copied:" << endl;
	      print_pars(pol3?4:3, pol3?params3:params2);
	      cout << "fit covariance matrix:" << endl;
	      print_covmat(pol3?7:6, fit_res->GetCovarianceMatrix().GetMatrixArray());
	      cout << "copied matrix:" << endl;
	      print_covmat(pol3?4:3, pol3?covmat3:covmat2);
	    }
	  }

	  //cout << "par4 = " << fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3) << endl;
	  // from fit
	  //yield[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3);
	  //yield_er[itu][iplab][nptok[itu][iplab]] = fmfg_tu[itu][iplab][itbin]->GetParError(pol3?4:3);

	  yield[itu][iplab][nptok[itu][iplab]] = yield_cnt_fg[itu][iplab][nptok[itu][iplab]] - fmbg_tu[itu][iplab][itbin]->Integral(mmin, mmax);
	  double e2 = fmbg_tu[itu][iplab][itbin]->IntegralError(mmin, mmax, pol3?params3:params2, pol3?covmat3:covmat2);
	  yield_er[itu][iplab][nptok[itu][iplab]] = TMath::Hypot(yield_cnt_fg_er[itu][iplab][nptok[itu][iplab]],e2);

	  if (verbose){
	    cout << " fgInt " << yield_cnt_fg[itu][iplab][nptok[itu][iplab]]
		 << " bgFuncInt= " << fmbg_tu[itu][iplab][itbin]->Integral(mmin, mmax)
		 << " bgFuncIntErr = " << e2
		 << " yieldNew = " << yield[itu][iplab][nptok[itu][iplab]] << " yieldOld= " << fmfg_tu[itu][iplab][itbin]->GetParameter(pol3?4:3) << endl;
	  }

	  tvalid[itu][iplab][nptok[itu][iplab]] = itu==0?t[iplab][itbin]:u[iplab][itbin];
	  _tvalid_min[itu][iplab][nptok[itu][iplab]] = _t_min[iplab][itbin];
	  _tvalid_max[itu][iplab][nptok[itu][iplab]] = _t_max[iplab][itbin];
	  yield[itu][iplab][nptok[itu][iplab]] /= tbin_width;
	  yield_er[itu][iplab][nptok[itu][iplab]] /= tbin_width;

	  stob[itu][iplab][nptok[itu][iplab]] = yield_bg[itu][iplab][nptok[itu][iplab]]!=0?yield[itu][iplab][nptok[itu][iplab]]/yield_bg[itu][iplab][nptok[itu][iplab]]:0.0;
	  stob_er[itu][iplab][nptok[itu][iplab]] =
	    calc_err_r(yield[itu][iplab][nptok[itu][iplab]], yield_bg[itu][iplab][nptok[itu][iplab]],
		       yield_er[itu][iplab][nptok[itu][iplab]], yield_bg_er[itu][iplab][nptok[itu][iplab]]);

	  bool cnt = false;
	  if (cnt){
	    yield_cor[itu][iplab][nptok[itu][iplab]] = yield_cnt[itu][iplab][nptok[itu][iplab]]/eff_cor[itu][iplab][itbin];
	    yield_cor_er[itu][iplab][nptok[itu][iplab]] =
	      calc_err_r(yield_cnt[itu][iplab][nptok[itu][iplab]], eff_cor[itu][iplab][itbin],
			 yield_cnt_er[itu][iplab][nptok[itu][iplab]], eff_cor_er[itu][iplab][itbin]) ;
	  } else {
	    yield_cor[itu][iplab][nptok[itu][iplab]] = yield[itu][iplab][nptok[itu][iplab]]/eff_cor[itu][iplab][itbin];
	    yield_cor_er[itu][iplab][nptok[itu][iplab]] =
	      calc_err_r(yield[itu][iplab][nptok[itu][iplab]], eff_cor[itu][iplab][itbin],
			 yield_er[itu][iplab][nptok[itu][iplab]], eff_cor_er[itu][iplab][itbin]) ;
	  }

	  yield_bg_cor[itu][iplab][nptok[itu][iplab]] = yield_bg[itu][iplab][nptok[itu][iplab]]/eff_cor[itu][iplab][itbin];
	  yield_bg_cor_er[itu][iplab][nptok[itu][iplab]] =
	    calc_err_r(yield_bg[itu][iplab][nptok[itu][iplab]], eff_cor[itu][iplab][itbin],
		       yield_bg_er[itu][iplab][nptok[itu][iplab]], eff_cor_er[itu][iplab][itbin]) ;
	  // contribution from pi0pippim
	  yield_bg0_cor[itu][iplab][nptok[itu][iplab]] = yield_bg0[itu][iplab][nptok[itu][iplab]]/eff_cor[itu][iplab][itbin];
	  yield_bg0_cor_er[itu][iplab][nptok[itu][iplab]] =
	    calc_err_r(yield_bg0[itu][iplab][nptok[itu][iplab]], eff_cor[itu][iplab][itbin],
		       yield_bg0_er[itu][iplab][nptok[itu][iplab]], eff_cor_er[itu][iplab][itbin]) ;

	  //cout << "t=  " << t[iplab][nptok[itu][iplab]] << " eff= " << eff_cor[itu][iplab][itbin]
	  // << " pm " << eff_cor_er[itu][iplab][itbin] << endl;

	  //yield_cor[itu][iplab][nptok[itu][iplab]] /= (msv?Lumi_msv[iplab]:Lumi_full[iplab])*Br;
	  //yield_cor_er[itu][iplab][nptok[itu][iplab]] /= (msv?Lumi_msv[iplab]:Lumi_full[iplab])*Br;

	}

	if (nptok[itu][iplab]<9) { // Drawing only here. we have 3x3

	  tc_mep_tbins[itu][iplab]->cd(nptok[itu][iplab]+1);
	  tc_mep_tbins[itu][iplab]->GetPad(nptok[itu][iplab]+1)->SetBottomMargin(0.15);
	  hmfg_tu[itu][iplab][itbin]->Draw();

	  hmfg_tu[itu][iplab][itbin]->GetXaxis()->SetRangeUser(2.3,3.8);
	  hmfg_tu[itu][iplab][itbin]->GetYaxis()->SetNdivisions(505);
	  //hmfg_tu[itu][iplab][itbin]->SetMinimum(0);
	  //hmfg_tu[itu][iplab][itbin]->SetMaximum(fmfg_tu[itu][iplab][itbin]->GetMaximum()*1.1);
	  hmsg_tu[itu][iplab][itbin]->SetLineStyle(9);
	  hmbg_tu[itu][iplab][itbin]->SetLineStyle(9);

	  hmsg_tu[itu][iplab][itbin]->Draw("hist,same");
	  hmbg_tu[itu][iplab][itbin]->Draw("hist,same");

	  fmfg_tu[itu][iplab][itbin]->Draw("same");
	  fmbg_tu[itu][iplab][itbin]->Draw("same");
	  fmsg_tu[itu][iplab][itbin]->Draw("same");

	  //tl[0][iplab]->DrawLatex(0.15,0.7,Form("%s",hmsg_tu[itu][iplab][itbin]->GetTitle()));
	  //tl[0][iplab]->DrawLatex(0.15,0.92,Form("%4.2f<t[GeV^{2}]<%4.2f", tu_bins[itbin], tu_bins[itbin+1]));

	  //tline->DrawLine(mmin,0,mmin,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());
	  //tline->DrawLine(mmax,0,mmax,0.75*hmfg_tu[itu][iplab][itbin]->GetMaximum());

	  //gPad->SetLogy();

	}

	nptok[itu][iplab]++;

      }

      //tc_mep_tbins[itu][iplab]->Print(Form("%s/figs/v2/%s/%s.pdf",bdir,(msv?"msv":"full"),tc_mep_tbins[itu][iplab]->GetName()));
      if (verbose) {
	for (int itbin=0; itbin < ntbin; ++itbin) {
	  cout << "itbin " << itbin << "/" << nptok[itu][iplab]
	    //<< " tvalid " << tvalid[itu][iplab][itbin] << " t_er[itbin]= " << t_er[itbin]
	       << " yield " << yield[itu][iplab][itbin]  << " yield_er " << yield_er[itu][iplab][itbin]
	       << " yield_cnt " << yield_cnt[itu][iplab][itbin] << " yield_cnt_er " << yield_cnt_er[itu][iplab][itbin]
	       << " stob " << stob[itu][iplab][itbin] << " stob_er " << stob_er[itu][iplab][itbin]
	       << " eff " << eff_cor[itu][iplab][itbin]
	       << " yield_cor " << yield_cor[itu][iplab][itbin] << " yield_cor_er " << yield_cor_er[itu][iplab][itbin] << endl;
	}
      }

      tg_yield_bg[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],yield_bg[itu][iplab],t_er,yield_bg_er[itu][iplab]);
      graph_to_hist(tg_yield_bg[itu][iplab],h_yield_bg[itu][iplab],d_yield_bg_min[itu][iplab],d_yield_bg_max[itu][iplab]);
      h_yield_bg[itu][iplab]->SetLineWidth(2);
      h_yield_bg[itu][iplab]->SetLineColor(col[iplab]);
      h_yield_bg[itu][iplab]->SetFillColor(col[iplab]);
      h_yield_bg[itu][iplab]->SetFillStyle(3004);

      tg_yield_bg_cor[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],yield_bg_cor[itu][iplab],t_er,yield_bg_cor_er[itu][iplab]);
      double ddummy1= 0, ddummy2= 0;
      graph_to_hist(tg_yield_bg_cor[itu][iplab],h_yield_bg_cor[itu][iplab],ddummy1,ddummy2);
      h_yield_bg_cor[itu][iplab]->SetLineWidth(2);
      h_yield_bg_cor[itu][iplab]->SetLineColor(col[iplab]);
      h_yield_bg_cor[itu][iplab]->SetFillColor(col[iplab]);
      h_yield_bg_cor[itu][iplab]->SetFillStyle(3004);

      tg_yield_bg0_cor[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],yield_bg0_cor[itu][iplab],t_er,yield_bg0_cor_er[itu][iplab]);
      graph_to_hist(tg_yield_bg0_cor[itu][iplab],h_yield_bg0_cor[itu][iplab],ddummy1,ddummy2);
      h_yield_bg0_cor[itu][iplab]->SetLineWidth(2);
      h_yield_bg0_cor[itu][iplab]->SetLineColor(col[iplab]);
      h_yield_bg0_cor[itu][iplab]->SetFillColor(col[iplab]);
      h_yield_bg0_cor[itu][iplab]->SetFillStyle(1001);

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
      legend[itu][iplab]->AddEntry(h_yield_bg[itu][iplab],"Background","plf");

      //tg_yield_cor[itu][iplab] = new TGraphErrors(ntbin-1, itu==0?t_cnt[iplab]:u_cnt[iplab], yield_cor[itu][iplab], t_er, yield_cor_er[itu][iplab]);
      tg_yield_cor[itu][iplab] = new TGraphErrors(nptok[itu][iplab], tvalid[itu][iplab], yield_cor[itu][iplab], t_er, yield_cor_er[itu][iplab]);
      tg_yield_cor[itu][iplab]->SetMarkerStyle(mar_cor[itu]);
      tg_yield_cor[itu][iplab]->SetMarkerSize(1.5);
      tg_yield_cor[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_yield_cor[itu][iplab]->SetLineColor(col[iplab]);
      tg_yield_cor[itu][iplab]->SetLineWidth(2);

      const char *hfit_name = Form("hcosthfit_p%d_%s%d",iplab,toru[itu/2],itu%2);
      hcostfit[itu][iplab] = new TH1F(hfit_name,hfit_name,8,-0.8,0.8);
      double d_fitmin, d_fitmax;
      for (int ipt=0; ipt < nptok[itu][iplab]; ++ipt) {
	double x,val,err;
	tg_yield_cor[itu][iplab]->GetPoint(ipt,x,val);
	err = tg_yield_cor[itu][iplab]->GetErrorY(ipt);
	int ibin = hcostfit[itu][iplab]->GetXaxis()->FindBin(x);
	hcostfit[itu][iplab]->SetBinContent(ibin,val);
	hcostfit[itu][iplab]->SetBinError(ibin,err);
	if (ipt==0) d_fitmin = hcostfit[itu][iplab]->GetBinLowEdge(ibin);
	if (ipt==nptok[itu][iplab]-1) d_fitmax = hcostfit[itu][iplab]->GetBinLowEdge(ibin) + hcostfit[itu][iplab]->GetBinWidth(ibin);
      }

      const char *func_name = Form("fcosth_p%d_%s%d",iplab,toru[itu/2],itu%2);
      fcosth_tu[itu][iplab] = new TF1(func_name, "[0]*(1+[1]*x*x)", d_fitmin, d_fitmax);
      fcosth_tu[itu][iplab]->SetLineColor(col[iplab]);
      fcosth_tu[itu][iplab]->SetParameter(1,1.0);

      hcostfit[itu][iplab]->Fit(func_name,"0Q");
      //tg_yield_cor[itu][iplab]->Fit(func_name, "0Q");
      fcosth_tu[itu][iplab]->SetParameter(0,fcosth_tu[itu][iplab]->GetParameter(0));
      fcosth_tu[itu][iplab]->SetParameter(1,fcosth_tu[itu][iplab]->GetParameter(1));
      hcostfit[itu][iplab]->Fit(func_name, "RINOEQ+");

      //tg_yield_cor[itu][iplab]->Fit(func_name, "RNOEQ+");
      //tg_yield_cor[itu][iplab]->Fit(func_name, "RNO");
      tg_yield_cor[itu][iplab]->GetListOfFunctions()->Add(fcosth_tu[itu][iplab]);

      if (write_dat)
	dat_out << Form("icth= %d itu= %d iplab= %d A= %5.2f #pm %5.2f",
			icth, itu, iplab, fcosth_tu[itu][iplab]->GetParameter(1),
			fcosth_tu[itu][iplab]->GetParError(1)) << endl;
      //cout << Form("icth= %d itu= %d iplab= %d A= %5.2f #pm %5.2f", icth, itu, iplab, fcosth_tu[itu][iplab]->GetParameter(1), fcosth_tu[itu][iplab]->GetParError(1)) << endl;

      tmg_yield_cor[itu]->Add(tg_yield_cor[itu][iplab],"p");
      //tmg_yield_cor[itu]->Add(tg_yield[itu][iplab],"p");
      tmg_yield_cor_pbp[itu][iplab]->Add(tg_yield_cor[itu][iplab],"p");
      legend2[itu][iplab]->AddEntry(tg_yield_cor[itu][iplab],"Eff. corrected yield","ep");
      legend2[itu][iplab]->AddEntry(h_yield_bg_cor[itu][iplab],"Background","f");

      tg_stob[itu][iplab] = new TGraphErrors(nptok[itu][iplab],tvalid[itu][iplab],stob[itu][iplab],t_er,stob_er[itu][iplab]);
      tg_stob[itu][iplab]->SetMarkerStyle(21);
      tg_stob[itu][iplab]->SetMarkerSize(1);
      tg_stob[itu][iplab]->SetMarkerColor(col[iplab]);
      tg_stob[itu][iplab]->SetLineColor(col[iplab]);
      tg_stob[itu][iplab]->SetLineWidth(2);
      tmg_stob_pbp[itu][iplab]->Add(tg_stob[itu][iplab],"p");
      legend3[itu][iplab]->AddEntry(tg_stob[itu][iplab],"S/B ratio","pl");

      if (write_dat) {
	tc_mep_tbins[itu][iplab]->Close();
      }

    }

    if (write_dat) {
      fsig[iplab]->Close();
      feff[iplab]->Close();
      fbg[0][iplab]->Close();
      fbg[1][iplab]->Close();
      fbg[2][iplab]->Close();
      fbg[3][iplab]->Close();
      dat_out.close();
    }
  }

  if (write_dat) {
    tctmp->Close();
  }

  if (write_dat){
    return;
  }

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
    if (itu==0||itu==2) continue;
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
    if (itu==0||itu==2) continue;
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

  //cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;
  //cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;
  //cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;
  //cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;

  bool fancy = true;
  for (int itu = 0; itu < ntu; ++itu) {

    if (itu ==0 || itu==2) continue;

    /*
    tc_stob_pbp[itu] = new TCanvas(Form("stob_pbp_%s%d",toru[itu/2],itu%2),Form("stob_pbp_%s%d",toru[itu/2],itu%2),1200,500);
    //tc_stob_pbp[itu] = new TCanvas(Form("stob_pbp%s",toru[itu]),Form("stob_pbp%s",toru[itu]));
    tc_stob_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      tc_stob_pbp[itu]->cd(iplab+1);
      tmg_stob_pbp[itu][iplab]->Draw("a");

      if (fancy) {
    	hdummy_stob[itu][iplab] = (TH1F*)tmg_stob_pbp[itu][iplab]->GetHistogram();

    	hdummy_stob[itu][iplab]->SetLabelSize(iplab==2?0.045:0.065,"Y");
    	hdummy_stob[itu][iplab]->SetLabelSize(0.065,"X");
    	hdummy_stob[itu][iplab]->SetLabelOffset(0.005,"Y");
    	hdummy_stob[itu][iplab]->SetTitleSize(0.06,"X");
    	hdummy_stob[itu][iplab]->SetTitleSize(iplab==2?0.045:(iplab==1?0.058:0.06),"Y");
    	hdummy_stob[itu][iplab]->SetTitleOffset(iplab==2?1.55:(iplab==1?1.1:1.0),"Y");

    	hdummy_stob[itu][iplab]->SetTitle(Form(";cos(#theta);Signal/Background          "));
    	//if (iplab==0)
    	//	hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background          ",(itu==0?"t":"u")));
    	//else if (iplab==1)
    	//	hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background          ",(itu==0?"t":"u")));
    	//else
    	//hdummy_stob[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];Signal/Background          ",(itu==0?"t":"u")));
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
    	//tl[3][iplab]->DrawLatex(_x_st+0.34,0.23,"S: #bar{p}p#rightarrow#pi^{0}J/#psi");
    	//tl[3][iplab]->DrawLatex(_x_st+0.34,0.17,"B: #bar{p}p#rightarrow#pi^{0}#pi^{+}#pi^{-}");
    	legend3[itu][iplab]->Draw();

      }

    }

    //tc_stob_pbp[itu]->Print(Form("%s/figs/v2/%s/%s.pdf",bdir,(msv?"msv":"full"),tc_stob_pbp[itu]->GetName()));
    //tc_stob_pbp[itu]->Print(Form("cth_%s_%s.pdf",tc_stob_pbp[itu]->GetName(),(msv?"msv":"full")));

    //tc_yield_cnt_pbp[itu] = new TCanvas(Form("fitted_yield_cnt_pbp_%s%d",toru[itu/2],itu%2),Form("fitted_yield_cnt_pbp_%s%d",toru[itu/2],itu%2));
    //tc_yield_cnt_pbp[itu]->Divide(3,1);
    //for (int iplab = 0; iplab < nplab; ++iplab) {
    //  tc_yield_cnt_pbp[itu]->cd(1+iplab);
    //  tmg_yield_cnt_pbp[itu][iplab]->Draw("a");
    //}
    //
    //return;

    */

    fancy = true;
    tc_yield_cnt_pbp[itu] = new TCanvas(Form("fitted_yield_cnt_pbp_%s%d",toru[itu/2],itu%2),Form("fitted_yield_cnt_pbp_%s%d",toru[itu/2],itu%2));
    //tc_yield_cnt_pbp[itu]->Divide(3,1);
    for (int iplab = 0; iplab < nplab; ++iplab) {
      double xl = (iplab==0?0:0.1)+iplab*0.3;
      double xh = 0.1+(iplab+1)*0.3+(iplab==1?0.001:0.0);
      pad_yield_cnt[itu][iplab] = new TPad(Form("pad_yield_%s%d_%d",toru[itu/2],itu%2,iplab),Form("pad_yield_%s%d_%d",toru[itu/2],itu%2,iplab),xl,0.0,xh,1.0);
      tc_yield_cnt_pbp[itu]->cd(0);
      pad_yield_cnt[itu][iplab]->Draw();

      if (fancy) {
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
	//hdummy_yield_cnt[itu][iplab]->SetTitle(Form(";%s[GeV^{2}];dN_{J/#psi}/d%s[GeV^{-2}]",(itu==0?"t":"u"),(itu==0?"t":"u")));
	hdummy_yield_cnt[itu][iplab]->SetTitle(Form(";cos(#theta);count"));
	hdummy_yield_cnt[itu][iplab]->SetMinimum(0);
	hdummy_yield_cnt[itu][iplab]->SetMaximum((max_yield[0]>max_yield[1]?max_yield[0]:max_yield[1])*1.2);

	tmg_yield_cnt_pbp[itu][iplab]->SetMinimum(0.0);
	tl[1][iplab]->DrawLatex(iplab==0?0.33:0.15,0.93,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

	h_yield_bg[itu][iplab]->GetXaxis()->SetRangeUser(d_yield_bg_min[itu][iplab]-0.005,d_yield_bg_max[itu][iplab]+0.005);
	h_yield_bg[itu][iplab]->Draw("same,hist");

	legend[itu][iplab]->Draw();
      }

    }
    //tc_yield_cnt_pbp[itu]->Print(Form("%s/figs/v2/%s/%s.pdf",bdir,(msv?"msv":"full"),tc_yield_cnt_pbp[itu]->GetName()));
    tc_yield_cnt_pbp[itu]->Print(Form("cth_%s_%s.pdf",tc_yield_cnt_pbp[itu]->GetName(),(msv?"msv":"full")));

    tc_yield_cor_pbp[itu] = new TCanvas(Form("fitted_yield_cor_pbp_%s%d",toru[itu/2],itu%2),Form("fitted_yield_cor_pbp_%s%d",toru[itu/2],itu%2), 10, itu/2==0?10:560, 1400, 550);
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

      hdummy_yield_cor[itu][iplab]->GetXaxis()->SetNdivisions(505);
      hdummy_yield_cor[itu][iplab]->GetYaxis()->SetNdivisions(505);

      hdummy_yield_cor[itu][iplab]->GetXaxis()->SetRangeUser(-1.,1.);
      //hdummy_yield_cor[itu][iplab]->SetTitle(Form("FIX THIS;%s[GeV^{2}];d#sigma_{J/#psi-#pi^{0}}/d%s[pb/GeV^{2}]",(itu==0?"t":"u"),(itu==0?"t":"u")));
      hdummy_yield_cor[itu][iplab]->SetTitle(Form(";cos(#theta);eff. corr. yield"));
      //if (iplab==1) {
      //	hdummy_yield_cor[itu][iplab]->GetXaxis()->SetNdivisions(505,false);
      //}
      //hdummy_yield_cor[itu][iplab]->SetMinimum(0);

      if (iplab==1)
	hdummy_yield_cor[itu][iplab]->SetMaximum(1990);

      tmg_yield_cor_pbp[itu][iplab]->SetMinimum(0.0);
      //TF1 *funkyfunk = get_func(iplab, tvalid[itu][iplab][0], tvalid[itu][iplab][nptok[itu][iplab]-1] );
      //funkyfunk->SetLineColor(col[iplab]);
      //funkyfunk->Draw("same");
      //legend2[itu][iplab]->AddEntry(funkyfunk,"TDA model","pl");

      h_yield_bg_cor[itu][iplab]->GetXaxis()->SetRangeUser(d_yield_bg_min[itu][iplab]-0.005,d_yield_bg_max[itu][iplab]+0.005);
      h_yield_bg_cor[itu][iplab]->Draw("same,hist");
      //h_yield_bg0_cor[itu][iplab]->GetXaxis()->SetRangeUser(d_yield_bg_min[itu][iplab]-0.005,d_yield_bg_max[itu][iplab]+0.005);
      //h_yield_bg0_cor[itu][iplab]->Draw("same,hist");

      legend2[itu][iplab]->Draw();
      fcosth_tu[itu][iplab]->Draw("same");

      tl[2][iplab]->DrawLatex(0.33,0.83,Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
      const char *txt = Form("(\"small %s\")",toru[itu/2]);
      //const char *txt = Form("(%.3f < %s [GeV^{2}] < %.2f)", tvalidmin[iplab], toru[itu/2], tvalidmax[iplab]);
      //if (iplab==1) txt = Form("(%.1f < %s [GeV^{2}] < %.2f)", tvalidmin[iplab], toru[itu/2], tvalidmax[iplab]);
      //if (iplab==2) txt = Form("(%.1f < %s [GeV^{2}] < %.1f)", tvalidmin[iplab], toru[itu/2], tvalidmax[iplab]);
      tl[5][iplab]->DrawLatex(0.45,0.76,txt);

      tl[4][iplab]->SetTextColor(col[iplab]);
      tl[4][iplab]->DrawLatexNDC(0.4,0.4,Form("A = %5.2f #pm %5.2f", fcosth_tu[itu][iplab]->GetParameter(1), fcosth_tu[itu][iplab]->GetParError(1)));

      cout << Form("itu= %d iplab= %d A = %5.2f #pm %5.2f", itu, iplab, fcosth_tu[itu][iplab]->GetParameter(1), fcosth_tu[itu][iplab]->GetParError(1)) << endl;

    }

    //tc_yield_cor_pbp[itu]->Print(Form("%s/figs/v2/%s/%s.pdf",bdir,(msv?"msv":"full"),tc_yield_cor_pbp[itu]->GetName()));
    tc_yield_cor_pbp[itu]->Print(Form("cth_%s_%s.pdf",tc_yield_cor_pbp[itu]->GetName(),(msv?"msv":"full")));

  }

}

void set_style_gen_dists(TH1* h, int col, int rebin=0, bool sumw2=false) {
  if (rebin>0)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetLabelSize(0.07);
  h->GetYaxis()->SetLabelOffset(0.02);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

double _Mpi = 0.135;
double _MN = 0.938;
double _Mj = 3.096;

double _s(double _plab) { return 2*_MN*_MN + 2*_MN*sqrt(_MN*_MN + _plab*_plab); }

double mirror(double _u, double _plab) {
  double msqTot = 2*_MN*_MN + _Mj*_Mj + _Mpi*_Mpi;
  return msqTot - _s(_plab) - _u;
}

double _pcm(double _plab) { return _plab * _MN / sqrt(_s(_plab)); }

double _q(double _plab) {
  double _num = sqrt( pow(_s(_plab),2) + pow(_Mj*_Mj - _Mpi*_Mpi, 2) - 2*_s(_plab)*(_Mj*_Mj + _Mpi*_Mpi) ) ;
  return _num / 2 / sqrt(_s(_plab));
}

double _tt(double _cth, double _plab) {
  double _det = sqrt( pow(_q(_plab),2) + _Mpi*_Mpi ) - sqrt( pow(_pcm(_plab),2) + _MN*_MN );
  return _det*_det - pow(_pcm(_plab),2) - pow(_q(_plab),2) + 2*_q(_plab)*_pcm(_plab)*_cth;
}

double _costhcm(double _t, double _plab) {
  double _det = sqrt( pow(_q(_plab),2) + _Mpi*_Mpi ) - sqrt( pow(_pcm(_plab),2) + _MN*_MN );
  double _num = _t + pow(_pcm(_plab),2) + pow(_q(_plab),2) - _det*_det;
  double _costhcm_ = _num / 2 / _q(_plab) / _pcm(_plab);
  if (_costhcm_>1) _costhcm_ = 1.0;
  if (_costhcm_<-1) _costhcm_ = -1.0;
  return _costhcm_;

}

double _thcm(double _t, double _plab) {
  double _costhcm_ = _costhcm(_t,_plab);
  return acos(_costhcm_);
}

double _beta_cm(double _plab) {
  double E_antip = TMath::Hypot(_MN, _plab);
  return _plab/(E_antip + _MN);
}

double _thlab(double _t, double _plab) {

  double _costhcm_ = _costhcm(_t, _plab);
  cout << "_costhcm_ = " << Form("%20.19f",_costhcm_) << endl;
  cout << "ACos(_costhcm_) = " << TMath::ACos(_costhcm_) << endl;
  double _thcm_ = _thcm(_t, _plab);
  cout << "_thcm_ = " << _thcm_ << endl;

  TVector3 boost_to_lab;
  boost_to_lab.SetZ(_beta_cm(_plab));
  //cout << "=============================" << endl;
  //cout << "t= " << _t << " plab= " << _plab << endl;

  TVector3 p3g;
  p3g.SetMagThetaPhi(1.0,_thcm_,0);
  TLorentzVector p4g(p3g,1.0);
  p4g.Boost(boost_to_lab);
  return p4g.Vect().Theta()*TMath::RadToDeg();

  //double _thelab_ = 0.0;
  //for (int ii=0; ii < 3; ++ii) {
  //  double _tmp_mom = 1.0*(ii+1);
  //
  //  TVector3 p3g;
  //  p3g.SetMagThetaPhi(_tmp_mom,_thcm_,0);
  //  TLorentzVector p4g(p3g,_tmp_mom);
  //  p4g.Boost(boost_to_lab);
  //  _thelab_ = p4g.Vect().Theta();
  //  cout << "Phot: mom= " << _tmp_mom << " philab= " << p4g.Vect().Phi() << " _costhlab= " << p4g.CosTheta() << " => the_lab(deg)= " << _thelab_*TMath::RadToDeg() << endl;
  //
  //}
  //
  //cout << "----------------------"<<endl;
  //for (int ii=0; ii < 5; ++ii) {
  //  double _tmp_mom = 1.0*(ii+1);
  //
  //  TVector3 p3pi;
  //  p3pi.SetMagThetaPhi(_tmp_mom,_thcm_,0.0);
  //  TLorentzVector p4pi;
  //  p4pi.SetPxPyPzE(p3pi.Px(), p3pi.Py(), p3pi.Pz(), sqrt(_Mpi*_Mpi + p3pi.Mag2() ) );
  //  p4pi.Boost(boost_to_lab);
  //  _thelab_ = p4pi.Vect().Theta();
  //  cout << "Pion: mom= " << _tmp_mom << " philab= " << p4pi.Vect().Phi() << " _costhlab= " << p4pi.CosTheta() << " => the_lab(deg)= " << _thelab_*TMath::RadToDeg() << endl;
  //
  //}
  //
  //cout << "----------------------"<<endl;
  //for (int ii=0; ii < 5; ++ii) {
  //  double _tmp_mom = 1.0*(ii+1);
  //
  //  TVector3 p3p;
  //  p3p.SetMagThetaPhi(_tmp_mom,_thcm_,0.0);
  //  TLorentzVector p4p(p3p,sqrt(_MN*_MN + p3p.Mag2()));
  //  //TLorentzVector p4p;
  //  //p4p.SetPxPyPzE(p3p.Px(), p3p.Py(), p3p.Pz(), sqrt(_MN*_MN + p3p.Mag2() ) );
  //  p4p.Boost(boost_to_lab);
  //  _thelab_ = p4p.Vect().Theta();
  //  cout << "Prot: mom= " << _tmp_mom << " philab= " << p4p.Vect().Phi() << " _costhlab= " << p4p.CosTheta() << " => the_lab(deg)= " << _thelab_*TMath::RadToDeg() << endl;
  //}
  //return _thelab_;
}

double costh_p0(double _t) {
  return _costhcm(_t, plab[0]);
}

double tt_p0(double _cth) {
  return _tt(_cth, plab[0]);
}

double _tmin_[3] = {-0.092, -1.0, -1.0};

void gen_dists() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  bool save = false;
  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";
  gROOT->LoadMacro(Form("%s/tda/tda.C",bdir));
  gROOT->LoadMacro(Form("%s/figs/pv0/ananote.C",bdir));

  //cout << "tmax= " << tmax[0] << " => costh_cm= " << costh_p0(tmax[0]) << endl;
  //cout << "tmin= " << mirror(tmax[0],plab[0]) << " => costh_cm= " << costh_p0(mirror(tmax[0],plab[0])) << endl;
  //TF1 *fth = new TF1("fth","costh_p0(x)",mirror(tmax[0],plab[0]),tmax[0]);
  //fth->Draw();
  //return;

  //TF1 *fth = new TF1("fth","tt_p0(x)",-1,1);
  //fth->Draw();
  //return;

  TFile *fsig[nplab], *fbg[nplab];
  TH1F* h_ttrue_sig[nplab], *h_ttrue_bg[nplab];
  TH1F* h_trec_sig[nplab], *h_trec_bg[nplab];

  TH1F* h_cthcm_true_sig[nplab], *h_cthcm_true_bg[nplab];
  TH1F* h_thlab_true_sig[nplab], *h_thlab_true_bg[nplab];

  TLatex *tl[6][nplab];
  for (int iplab=0; iplab<nplab; ++iplab) {
    fsig[iplab] = TFile::Open(Form("%s/hists/note.aug.2015/anav2_jpsi_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));
    fbg[iplab] = TFile::Open(Form("%s/hists/note.aug.2015/anav2_pip_pim_%s_plab%3.1f.root",bdir,(ibrem==0?"raw":"brem"), plab[iplab]));
    h_ttrue_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/httrumc")->Clone(Form("httruemc_sig_p%d",iplab));
    h_ttrue_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/httrumc")->Clone(Form("httruemc_bg_p%d",iplab));
    h_trec_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrecgg")->Clone(Form("htrecgg_sig_p%d",iplab));
    h_trec_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrecgg")->Clone(Form("htrecgg_bg_p%d",iplab));
    set_style(h_ttrue_sig[iplab], 1, 5);
    set_style(h_ttrue_bg[iplab], 1, 5);
    set_style(h_trec_sig[iplab], 2, 5);
    set_style(h_trec_bg[iplab], 2, 5);

    h_cthcm_true_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrupi0costhcm")->Clone(Form("h_cthcm_sig_p%d",iplab));
    h_thlab_true_sig[iplab] = (TH1F*) fsig[iplab]->Get("tu/htrupi0thlab")->Clone(Form("h_thlab_sig_p%d",iplab));
    h_cthcm_true_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrupi0costhcm_mc")->Clone(Form("h_cthcm_bg_p%d",iplab));
    h_thlab_true_bg[iplab] = (TH1F*) fbg[iplab]->Get("tu/htrupi0thlab_mc")->Clone(Form("h_thlab_bg_p%d",iplab));
    set_style(h_cthcm_true_sig[iplab], 1, 5);
    set_style(h_thlab_true_sig[iplab], 1, 5);
    set_style(h_cthcm_true_bg[iplab], 1, 5);
    set_style(h_thlab_true_bg[iplab], 1, 5);

    h_cthcm_true_sig[iplab]->SetTitle(";cos{#theta_{CM}}");
    h_thlab_true_sig[iplab]->SetTitle(";#theta_{LAB}");
    h_cthcm_true_bg[iplab]->SetTitle(";cos{#theta_{CM}}");
    h_thlab_true_bg[iplab]->SetTitle(";#theta_{LAB}");

    for (int i=0; i<6; ++i) {
      tl[i][iplab] =new TLatex();
      tl[i][iplab]->SetNDC(true);
      tl[i][iplab]->SetLineColor(i==0?2:1);
      tl[i][iplab]->SetTextSize(((i==2)?1.4:1.5)*tl[i][iplab]->GetTextSize());
    }
  }

  //TCanvas *tc_gen_t_sig = new TCanvas("gen_t_dist_sig","gen_t_dist_sig",1600,800);
  //tc_gen_t_sig->Divide(3);
  //TCanvas *tc_gen_t_bg = new TCanvas("gen_t_dist_bg","gen_t_dist_bg",1600,800);
  //tc_gen_t_bg->Divide(3);
  //for (int iplab = 0; iplab < nplab; ++iplab) {
  //  tc_gen_t_sig->cd(1+iplab);
  //  h_ttrue_sig[iplab]->SetTitle(";t[GeV^{2}]");
  //  h_ttrue_sig[iplab]->Draw();
  //  tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  //  tc_gen_t_bg->cd(1+iplab);
  //  h_ttrue_bg[iplab]->SetTitle(";t[GeV^{2}]");
  //  h_ttrue_bg[iplab]->Draw();
  //  tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  //}
  //
  //if (save) {
  //  tc_gen_t_sig->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_sig->GetName()));
  //  tc_gen_t_bg->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_bg->GetName()));
  //}
  //
  //TCanvas *tc_gen_cthcm_sig = new TCanvas("gen_cthcm_dist_sig","gen_cthcm_dist_sig",1600,800);
  //tc_gen_cthcm_sig->Divide(3);
  //TCanvas *tc_gen_cthcm_bg = new TCanvas("gen_cthcm_dist_bg","gen_cthcm_dist_bg",1600,800);
  //tc_gen_cthcm_bg->Divide(3);
  //for (int iplab = 0; iplab < nplab; ++iplab) {
  //  tc_gen_cthcm_sig->cd(1+iplab);
  //  h_cthcm_true_sig[iplab]->SetTitle(";cos(#theta_{CM})");
  //  h_cthcm_true_sig[iplab]->Draw();
  //  tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  //  tc_gen_cthcm_bg->cd(1+iplab);
  //  h_cthcm_true_bg[iplab]->SetTitle(";cos(#theta_{CM})");
  //  h_cthcm_true_bg[iplab]->Draw();
  //  tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  //}
  //
  //if (save) {
  //  tc_gen_cthcm_sig->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_cthcm_sig->GetName()));
  //  tc_gen_cthcm_bg->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_cthcm_bg->GetName()));
  //}

  TCanvas *tc_gen_thlab_sig = new TCanvas("gen_thlab_dist_sig","gen_thlab_dist_sig",1600,800);
  tc_gen_thlab_sig->Divide(3);

  //TCanvas *tc_gen_thlab_bg = new TCanvas("gen_thlab_dist_bg","gen_thlab_dist_bg",1600,800);
  //tc_gen_thlab_bg->Divide(3);

  TH1F *h_thlab_true_sig_deg[nplab];
  double thlab_min[nplab]={0};
  double thlab_max[nplab]={0};

  for (int iplab=0; iplab < nplab; ++iplab) {
    int _numbin = h_thlab_true_sig[iplab]->GetXaxis()->GetNbins();
    h_thlab_true_sig_deg[iplab] = new TH1F(Form("h_thlab_true_sig_deg_p%d",iplab),Form(";#theta_{LAB}[deg]; counts"),_numbin,0,180);
    set_style_gen_dists(h_thlab_true_sig_deg[iplab], 1, 0);
    for (int ibin=0; ibin <= _numbin; ++ibin) {
      h_thlab_true_sig_deg[iplab]->SetBinContent(ibin,h_thlab_true_sig[iplab]->GetBinContent(ibin));
    }
  }

  for (int iplab = 0; iplab < nplab; ++iplab) {
    tc_gen_thlab_sig->cd(1+iplab);

    double fwd_thmax = _thlab(_tmin_[iplab], plab[iplab]);
    //_thlab(tmax[iplab], plab[iplab]);
    double fwd_thmin = _thlab(_tt(1.0, plab[iplab]), plab[iplab]);

    double bwd_thmin = _thlab(mirror(_tmin_[iplab], plab[iplab]), plab[iplab]);
    //_thlab(mirror(tmax[iplab], plab[iplab]), plab[iplab]);
    double bwd_thmax = _thlab(_tt(-1.0, plab[iplab]), plab[iplab]);

    if (true){
      gPad->SetRightMargin(0.05);
      gPad->SetTopMargin(0.05);
      //TGaxis::SetMaxDigits(3);

      //h_thlab_true_sig_deg[iplab]->GetXaxis()->SetRange(0.0,179);
      h_thlab_true_sig_deg[iplab]->Draw();

      TBox *b_fwd = new TBox(fwd_thmin, 0, fwd_thmax, 1.05*h_thlab_true_sig_deg[iplab]->GetMaximum());
      b_fwd->SetFillColor(kRed-10);
      b_fwd->Draw();

      TBox *b_bwd = new TBox(bwd_thmin, 0, bwd_thmax, 1.05*h_thlab_true_sig_deg[iplab]->GetMaximum());
      b_bwd->SetFillColor(kCyan-10);
      b_bwd->Draw();

      h_thlab_true_sig_deg[iplab]->Draw("same");
      h_thlab_true_sig_deg[iplab]->GetXaxis()->SetNdivisions(603,false);
      h_thlab_true_sig_deg[iplab]->GetYaxis()->SetNdivisions(505);
      //tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
      tl[0][iplab]->DrawLatex(0.35,0.8,Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

      if (iplab==0) {
	TText *tt = new TText();
	tt->SetTextSize(0.08);
	//tt->DrawText(2.0, 0.15*h_thlab_true_sig_deg[iplab]->GetMaximum(),"Fwd Kin");
	//tt->DrawText(2.0, 0.07*h_thlab_true_sig_deg[iplab]->GetMaximum(),"Valid. range");
	tt->SetTextAngle(90);
	//tt->DrawText(iplab!=2?12.0:7.0, 0.02*h_thlab_true_sig_deg[iplab]->GetMaximum(),"Fwd Kin");
	tt->DrawText(18.0, 0.02*h_thlab_true_sig_deg[iplab]->GetMaximum(),"Fwd. Kin");

	tt->SetTextAngle(0);
	//TText *tt = new TText();
	tt->DrawText(90, 0.3*h_thlab_true_sig_deg[iplab]->GetMaximum(),"Bwd. Kin");
	//tt->DrawText(120, 0.12*h_thlab_true_sig_deg[iplab]->GetMaximum(),"Valid. range");
      }
      gPad->RedrawAxis();

    } else {
      h_thlab_true_sig[iplab]->SetTitle(";#theta_{LAB}[rad]");
      h_thlab_true_sig[iplab]->Draw();
      tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
    }

    //tc_gen_thlab_bg->cd(1+iplab);
    //h_thlab_true_bg[iplab]->SetTitle(";#theta_{LAB}[rad]");
    //h_thlab_true_bg[iplab]->Draw();
    //tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

  }

  if (save) {
    tc_gen_thlab_sig->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_thlab_sig->GetName()));
    //tc_gen_thlab_bg->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_thlab_bg->GetName()));
  }

  //TCanvas *tc_gen_t_sig2 = new TCanvas("gen_rec_t_dist_sig","gen_rec_t_dist_sig",1600,800);
  //tc_gen_t_sig2->Divide(3);
  //TCanvas *tc_gen_t_bg2 = new TCanvas("gen_rec_t_dist_bg2","gen_rec_t_dist_bg",1600,800);
  //tc_gen_t_bg2->Divide(3);
  //for (int iplab = 0; iplab < nplab; ++iplab) {
  //  tc_gen_t_sig2->cd(1+iplab);
  //  h_ttrue_sig[iplab]->Draw();
  //  h_trec_sig[iplab]->Draw("same");
  //  tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  //  tc_gen_t_bg2->cd(1+iplab);
  //  h_ttrue_bg[iplab]->Draw();
  //  h_trec_bg[iplab]->Draw("same");
  //  tl[0][iplab]->DrawLatex(0.25,0.8,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  //}
  //
  //if (save) {
  //  tc_gen_t_sig2->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_sig2->GetName()));
  //  tc_gen_t_bg2->Print(Form("%s/figs/2015.09.15/%s.pdf", bdir, tc_gen_t_bg2->GetName()));
  //}

  //TCanvas *tc_gen_rec_t_sig = new TCanvas("gen_rec_t_dist_sig","gen_rec_t_dist_sig");
  //tc_gen_rec_t_sig->Divide(3);
  //TCanvas *tc_gen_rec_t_bg = new TCanvas("gen_rec_t_dist_bg","gen_rec_t_dist_bg");
  //tc_gen_rec_t_sig->Divide(3);
  //for (int iplab = 0; iplab < nplab; ++iplab) {
  //  cout << iplab << endl;
  //  tc_gen_rec_t_sig->cd(1+iplab);
  //  h_ttrue_sig[iplab]->Draw();
  //  //h_trec_sig[iplab]->SetLineStyle(9);
  //  h_trec_sig[iplab]->Draw("same");
  //  tc_gen_rec_t_bg->cd(1+iplab);
  //  h_ttrue_bg[iplab]->Draw();
  //  //h_trec_bg[iplab]->SetLineStyle(9);
  //  h_trec_bg[iplab]->Draw("same");
  //}

}

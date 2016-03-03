void set_style(TH1* h, int col, int rebin=0, bool sumw2=false) {
  if (rebin>0)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void a(int iwt=0) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  //int ncth = 28;
  string d0,d1,d2,d3,d4;
  int i0,_itu,_iplab;
  double a, a_er;

  double plab[3] = {5.513,8.0,12.0};
  const char *toru[2] = {"t","u"};
  const char *_wt[3] = {"wt0","wt1","nowt"};

  double h_min[3] = {0,-0.6,-1};
  double h_max[3] = {2,1.4,1};
  double f_min[3] = {0.4,-0.2,-0.6};
  double f_max[3] = {1.4,0.8,0.6};
  double a_true[3] = {1.0, 0.4, 0.0};

  const int ntu = 2;
  const int nplab = 3;
  TH1F *h_a[ntu][nplab];
  for (int itu=0; itu < ntu; ++itu) {
    for (int iplab=0; iplab < nplab; ++iplab) {
      h_a[itu][iplab] = new TH1F(Form("h_a_itu%d_ip%d",itu,iplab),Form("%s, p_{#bar{p}} = %3.1f GeV/c;A",toru[itu],plab[iplab]),20,h_min[iwt],h_max[iwt]);
      h_a[itu][iplab]->GetXaxis()->SetNdivisions(505,false);
      set_style(h_a[itu][iplab], 4, 0, false);
    }
  }

  for (int icth=0; icth < 80; ++icth) {
    //cout << "icth= " << icth << endl;
    ifstream inf;
    inf.open(Form("full/%s/data_table_cth%d.dat", _wt[iwt], icth));
    while (true) {
      inf >> d0 >> i0 >> d1 >> _itu >> d2 >> _iplab >> d3 >> a >> d4 >> a_er;
      if (!inf.good()) break;
      //cout << " a= " << a << endl;
      int ncth = (_iplab==0)? 100:60;
      if (icth<ncth) h_a[_itu/2][_iplab]->Fill(a);
    }
    inf.close();
  }

  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(3,2);
  for (int itu=0; itu < ntu; ++itu) {
    for (int iplab=0; iplab < nplab; ++iplab) {
      tc->cd(1+3*itu+iplab);
      //h_a[itu][iplab]->Draw();
      const char *fname = Form("a_gaus_itu%d_ip%d",itu/2,iplab);
      TF1* f_a = new TF1(fname,"gaus",f_min[iwt],f_max[iwt]);

      h_a[itu][iplab]->Draw();

      TLatex *tlat = new TLatex();
      tlat->SetTextSize(0.06);
      //tlat->SetTextSize(0.07);
      tlat->DrawLatexNDC(0.55,0.8,Form("A_{true} = %3.1f", a_true[iwt]));
      tlat->DrawLatexNDC(0.53,0.7,Form("<M>=%4.2f", h_a[itu][iplab]->GetMean()));
      tlat->DrawLatexNDC(0.53,0.6,Form("  #sigma=%4.2f#pm%4.2f", h_a[itu][iplab]->GetRMS()));

      //h_a[itu][iplab]->Fit(fname,"RQO");
      //tlat->DrawLatexNDC(0.53,0.8,Form("<M>=%4.2f#pm%4.2f", f_a->GetParameter(1), f_a->GetParError(1) ));
      //tlat->DrawLatexNDC(0.53,0.7,Form("  #sigma=%4.2f#pm%4.2f", f_a->GetParameter(2), f_a->GetParError(2) ));
      //cout << plab[iplab] << " \& \$" << toru[itu] << "\$ \& "<<  a_true[iwt] << " \& "
      //	   << Form("%4.2f",f_a->GetParameter(1)) << " \& " << Form("%4.2f",f_a->GetParameter(2))
      //	   << " \& " << Form("%4.2f",f_a->GetParameter(2)/f_a->GetParameter(1)) <<  " \\\\ " << endl;

      cout << plab[iplab] << " \& \$" << toru[itu] << "\$ \& "<<  a_true[iwt] << " \& "
	   << Form("%4.2f",h_a[itu][iplab]->GetMean()) << " \& " << Form("%4.2f",h_a[itu][iplab]->GetRMS())
	   << " \& " << Form("%4.2f",h_a[itu][iplab]->GetRMS()/h_a[itu][iplab]->GetMean()) <<  " \\\\ " << endl;
    }
  }

  tc->Print(Form("a_fit_cumul_%s.pdf",_wt[iwt]));
}

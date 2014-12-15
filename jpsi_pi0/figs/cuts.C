void cuts(string fname) {

  double mmin = 1.5, mmax = 3.5;
  int immin = 0, immax = 0;
  double tot = 10000;
  TFile *f = new TFile(fname.c_str());
  const int nstep = 5;
  double x[nstep] = {1, 2 ,3, 4, 5};
  double eff[nstep]= {0.0};
  TH1F* hmep[nstep];
  TCanvas *tc = new TCanvas("tc","tc");
  TLegend *tl = new TLegend(0.15,0.5,0.5,0.7);
  for (int ii = 0; ii < nstep; ++ii) {
    hmep[ii] = (TH1F*) f->Get(Form("hmep_%d",ii));
    hmep[ii]->SetLineWidth(2);
    hmep[ii]->SetLineColor(1+ii);
    tl->AddEntry(hmep[ii],Form("step %d",ii), "l");
    //hmep[ii]->Draw((ii==0||ii==1)?"":"same");
    hmep[ii]->Draw((ii==0)?"":"same");
    if (immin==0||immax==0) {
      immin = hmep[ii]->GetXaxis()->FindBin(mmin);
      immax = hmep[ii]->GetXaxis()->FindBin(mmax);
    }
    double num = hmep[ii]->Integral(immin, immax);
    eff[ii] = num/tot;
    cout << "eff(step=" << ii << ")= " << eff[ii] << endl;
  }
  tl->Draw();
  TCanvas *tceff = new TCanvas("tceff","tceff");
  tceff->cd();
  TGraph *tg = new TGraph(nstep,x,eff);
  tg->SetMarkerSize(2);
  tg->SetMarkerStyle(20);
  tg->SetMarkerColor(2);
  tg->SetMinimum(0);
  tg->Draw();

  const int nbinth = 12;
  TH1F *hmep_pi0cost_cm[nbinth];
  double xcost[nbinth],xcost_e[nbinth];
  for (int ib = 0; ib < nbinth; ++ib) {
    double low = -1.0 + ( 2.0*ib/nbinth );
    double high = -1.0 + ( 2.0*(ib+1)/nbinth );
    xcost[ib] = (low+high)/2.0;
    xcost_e[ib] = 0.0;
  }
  double yield[nbinth], yield_e[nbinth];
  for (int ib = 0; ib < nbinth; ++ib) {
    hmep_pi0cost_cm[ib] = (TH1F*)f->Get(Form("pi0cost_cm_bins/hmep_pi0cost_cm_%d", ib));
    yield[ib] = hmep_pi0cost_cm[ib]->GetEntries();
    yield_e[ib] = sqrt(yield[ib]);
  }
  TCanvas *tcyield = new TCanvas("tcyield","tcyield");
  tcyield->cd();
  TGraphErrors *tg_yield = new TGraphErrors(nbinth,xcost,yield,xcost_e,yield_e);
  tg_yield->SetMarkerSize(2);
  tg_yield->SetMarkerStyle(20);
  tg_yield->SetMarkerColor(2);
  tg_yield->SetLineColor(2);
  tg_yield->SetLineWidth(2);
  tg_yield->SetMinimum(0);
  tg_yield->Draw();

}

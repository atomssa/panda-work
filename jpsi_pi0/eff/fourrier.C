double s(double x, double f, double i) {
  return TMath::Sin(f*x*i);
}

double c(double x, double f, double i) {
  return TMath::Cos(f*x*i);
}

double func(double *x, double *p) {
  double xx=x[0];
  double f1 = p[0]+  p[2]*s(xx,p[1],1)+   p[3]*s(xx,p[1],2)+  p[4]*c(xx,p[1],1) + p[5]*c(xx,p[1],2);
  double t1 = f1/(p[6]+(p[7]*TMath::Power(xx,p[8])));
  double p2 = p[9]+  p[10]*x[0]+  p[11]*x[0]*x[0];
  double t2 = p[12]*p2;
  return p[13]+t1+t2;
}

void set_style(TGraph* tmg) {
  TH1F* h = tmg->GetHistogram();
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  h->GetXaxis()->SetNdivisions(610);
}

double fourrier() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  double pars[] = { 6.35356e+00,
		    1.0,
		    4.13113e+00,
  		    -4.43669e+00,
		    0.1,
		    0.01,
		    0.1,
		    9.64513e+03,
		    1.22279e+00,
  		    4.66147e-04,
  		    2.96494e-05,
  		    -6.21090e-06,
		    1.0,
		    -3.23049e-06
  };

  TFile *f = TFile::Open("hadd_out/hadd.pi.root");
  TEfficiency* eff = (TEfficiency*) f->Get("prob_cut_9/eff1d_mom_e_id");

  TFile *f2 = TFile::Open("../new_test/eid90pct/anav2_pip_pim_brem_plab5.5.root");
  TEfficiency* eff2 = (TEfficiency*) f2->Get("pi_eff1d_total_smooth1d_clone");
  TList *list = (TList*) eff2->GetListOfFunctions();
  TF1* fsmth = (TF1*) list->First();
  TCanvas *tc = new TCanvas("tc","tc");
  tc->cd();
  eff->Draw();
  //
  tc->Update();
  //eff->GetPaintedGraph()->GetHistogram()->SetMinimum(0);
  //eff->GetPaintedGraph()->GetHistogram()->SetMaximum(0.006);
  eff->GetPaintedGraph()->SetMinimum(5e-6);
  eff->GetPaintedGraph()->SetMaximum(0.007);
  eff->Draw();
  eff->SetTitle(Form("%s;p_{MC}[GeV/c];#varepsilon(#pi^{#pm})^{EID}",eff->GetTitle()));
  set_style(eff->GetPaintedGraph());
  fsmth->Draw("same");
  gPad->SetLogy();
  tc->Update();
  TLegend *tl = new TLegend(0.25,0.2,0.65,0.4);
  tl->AddEntry(eff,"#pi^{#pm} mis-id eff");
  tl->AddEntry(fsmth,"parametrization");
  tl->SetBorderSize(0);
  tl->Draw();

  //TF1* f1 = new TF1("f1",func,0.0001,12,14);
  //for (int ii=0; ii < 14; ++ii) {
  //  if (ii==9||ii==10||ii==11)
  //    f1->FixParameter(ii,pars[ii]);
  //  else
  //    f1->SetParameter(ii,pars[ii]);
  //}
  //eff->Fit(f1,"+RME");
  //eff->Draw();
  //cout << "Eff (0.01)= " << f1->Eval(0.01) << endl;
  //cout << "Eff (0.1)= " << f1->Eval(0.1) << endl;
  //cout << "Eff (0.2)= " << f1->Eval(0.2) << endl;
  //cout << "Eff (0.5)= " << f1->Eval(0.5) << endl;
  //cout << "Eff (1.0)= " << f1->Eval(1.0) << endl;
  //TF1* f2 = new TF1("f1",func,0.0001,12,13);
  //for (int ii=0; ii < 13; ++ii) {
  //  if (ii!=12)
  //    f2->SetParameter(ii,f1->GetParameter(ii));
  //  else
  //    f2->SetParameter(ii,100*f1->GetParameter(ii));
  //}
  //f2->SetLineColor(4);
  //f2->Draw("same");

}

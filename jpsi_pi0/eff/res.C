res() {

  //TFile *f = TFile::Open("pim_res_effhists.root");
  //TFile *f = TFile::Open("hadd_out/hadd.pi.root");
  TFile *f = TFile::Open("hadd_out/hadd.e.root");

  TCanvas *tc[3];
  for (int d=0; d<3; ++d) {
    tc[d] = new TCanvas(Form("tc%d",d),Form("tc%d",d));
    tc[d]->Divide(8,5);
  }

  const int npbin = 40;

  TF1* ff[3][npbin];
  TH1F* res[3][npbin];
  TH2F* h_p[3];
  h_p[0] = (TH2F*) f->Get("h_dpx");
  h_p[1] = (TH2F*) f->Get("h_dpy");
  h_p[2] = (TH2F*) f->Get("h_dpz");

  int col[3] = {1,2,4};
  double x[npbin] = {0.}, xe[npbin] = {0.};
  double sig[3][npbin] = {{0.}};
  double sige[3][npbin] = {{0.}};
  TGraphErrors *tge[3];
  TMultiGraph *tmg = new TMultiGraph("tmg","tmg");
  TLine *tl1 = new TLine();
  TLine *tl2 = new TLine();
  for (int d=0; d<3; ++d) {
    for (int p=0; p<npbin; ++p) {
      res[d][p] = (TH1F*) h_p[d]->ProjectionX(Form("h_d%d_p%d",d,p), p*5, (p+1)*5);
      double range = p<21||d<2?0.2:0.8;
      res[d][p]->GetXaxis()->SetRangeUser(-range,range);
      double fitmax = p<21||d<2?0.05:0.15;
      ff[d][p] = new TF1(Form("f_d%d_p%d",d,p),"gaus",-fitmax,fitmax);
      tc[d]->cd(p+1);
      res[d][p]->Fit(ff[d][p], "+RIME");
      sig[d][p] = ff[d][p]->GetParameter(2)*5;
      sige[d][p] = ff[d][p]->GetParError(2)*5;
      x[p] = 0.25+(0.5*p);
      xe[p] = 0;
      double min= ff[d][p]->GetParameter(1) - sig[d][p];
      double max= ff[d][p]->GetParameter(1) + sig[d][p];
      tl1->DrawLine(min,0,min,res[d][p]->GetMaximum());
      tl2->DrawLine(max,0,max,res[d][p]->GetMaximum());
    }
    tge[d] = new TGraphErrors(npbin,x,sig[d],xe,sige[d]);
    tge[d]->SetMarkerStyle(20);
    tge[d]->SetMarkerColor(col[d]);
    tge[d]->SetLineColor(col[d]);
    tmg->Add(tge[d],"pl");
  }
  TCanvas *tcg = new TCanvas("tcg","tcg");
  tmg->Draw("a");

  cout << "static const double dpx_min["<< npbin << "] = { ";
  for (int p=0; p<npbin; ++p) {
    cout << ff[0][p]->GetParameter(1) - sig[0][p];
    if (p!=npbin-1) cout << ", ";
  }
  cout << "};" << endl;
  cout << "static const double dpx_max["<< npbin << "] = { ";
  for (int p=0; p<npbin; ++p) {
    cout << ff[0][p]->GetParameter(1) + sig[0][p];
    if (p!=npbin-1) cout << ", ";
  }
  cout << "};" << endl;
  cout << "static const double dpy_min["<< npbin << "] = { ";
  for (int p=0; p<npbin; ++p) {
    cout << ff[1][p]->GetParameter(1) - sig[1][p];
    if (p!=npbin-1) cout << ", ";
  }
  cout << "};" << endl;
  cout << "static const double dpy_max["<< npbin << "] = { ";
  for (int p=0; p<npbin; ++p) {
    cout << ff[1][p]->GetParameter(1) + sig[1][p];
    if (p!=npbin-1) cout << ", ";
  }
  cout << "};" << endl;
  cout << "static const double dpz_min["<< npbin << "] = { ";
  for (int p=0; p<npbin; ++p) {
    cout << ff[2][p]->GetParameter(1) - sig[2][p];
    if (p!=npbin-1) cout << ", ";
  }
  cout << "};" << endl;
  cout << "static const double dpz_max["<< npbin << "] = { ";
  for (int p=0; p<npbin; ++p) {
    cout << ff[2][p]->GetParameter(1) + sig[2][p];
    if (p!=npbin-1) cout << ", ";
  }
  cout << "};" << endl;


}

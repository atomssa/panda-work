void draw_pad4(TCanvas *tc, int ic, const char* fname_mu, const char* fname_e, double pt, double thmin, double thmax) {

  int ity = 1;
  const char* histn;
  const char* histt;
  if (ity==0) {
    //rec_mu->SetTitle(Form("p_{T}= %4.2f GeV/c, %3.0f < #theta < %3.0f; (p_{MC} - p_{KF})/p_{MC};counts",pt,thmin,thmax));
    histn = "h_rec";
    histt = Form("; (p_{MC} - p_{KF})/p_{MC};counts",pt,thmin,thmax);
  } else if (ity==1) {
    histn = "h_res_phi";
    histt = Form("; #phi_{MC} - #phi_{KF}(deg);counts",pt,thmin,thmax);
  } else {
    histn = "h_res_the";
    histt = Form("; (p_{MC} - p_{KF})/p_{MC};counts",pt,thmin,thmax);
  }

  ////const char *histn = "h_rec";
  //const char *histn = "h_res_phi";
  ////const char *histn = "h_res_the";

  int ipad  = ic+1;
  TFile *fmu = TFile::Open(fname_mu);
  TH1F* rec_mu = ((TH1F*) fmu->Get(histn))->Clone(Form("%s_mu",histn));
  rec_mu->Sumw2();
  rec_mu->SetLineWidth(2);
  rec_mu->SetLineColor(1);
  //rec_mu->SetLineStyle(6);

  cout << "Opening file " << fname_e << endl;
  TFile *f = TFile::Open(fname_e);

  TH1F* rec_e = ((TH1F*) f->Get(histn))->Clone(Form("%s_e",histn));
  rec_e->SetLineWidth(2);
  rec_e->SetLineColor(2);
  //rec_e->SetLineStyle(6);
  rec_e->Sumw2();

  tc->cd(1+ic);
  //rec_mu->Scale(double(rec_e->Integral())/double(rec_mu->Integral()));
  rec_mu->Scale(double(rec_e->GetEntries())/double(rec_mu->GetEntries()));
  rec_mu->SetTitle("");
  rec_mu->GetXaxis()->SetRangeUser(-0.5,1.0);

  rec_mu->SetTitle(histt);
  //rec_mu->SetTitle(Form("; #theta_{MC} - #theta_{KF}(deg);counts",pt,thmin,thmax));
  //rec_mu->GetYaxis()->SetTitleOffset(1);
  TGaxis::SetMaxDigits(3);
  rec_mu->Draw("hist");
  rec_e->Draw("hist,same");
  gPad->SetGridx();

  TLatex *tt = new TLatex();
  tt->SetNDC(kTRUE);
  tt->SetTextSize(0.08);

  tt->SetTextColor(2);
  tt->DrawLatex(0.5,0.75,Form("Electron"));

  tt->SetTextColor(1);
  tt->DrawLatex(0.55,0.65,Form("Muon"));

  tt->SetTextColor(1);
  tt->SetTextSize(0.05);
  tt->DrawLatex(0.5,0.56,Form("p_{T} = %3.1f GeV/c",pt));
  tt->DrawLatex(0.5,0.50,Form(" %3.0f#circ < #theta < %3.0f#circ",thmin,thmax));

}

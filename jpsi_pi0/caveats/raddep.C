void raddep(){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  //int iconfig = 8; double max = 0.2;
  int iconfig = 2; double max = 0.2;

  Double_t config[22][3] = {
    {0.15, 30.0, 45.0}, {0.25, 30.0, 45.0}, {0.5, 30.0, 45.0}, {0.5, 45.0, 60.0},
    {0.5, 60.0, 75.0}, {0.5, 75.0, 90.0}, {0.5, 90.0, 105.0}, {0.5, 105.0, 120.0},
    {1.0, 30.0, 45.0}, {1.0, 45.0, 60.0}, {1.0, 60.0, 75.0}, {1.0, 75.0, 90.0},
    {1.5, 30.0, 45.0}, {1.5, 45.0, 60.0}, {2.0, 30.0, 45.0}, {0.15, 10.0, 20.0},
    {0.25, 10.0, 20.0}, {0.5, 10.0, 20.0}, {1.0, 10.0, 20.0}, {1.5, 10.0, 20.0},
    {2.0, 10.0, 20.0}, {2.5, 10.0, 20.0}};

  TFile *f = new TFile(Form("../grid.out/psim_oct14_binsong_config%d/bremcorr_merged_old_hists.root",iconfig));
  TH1F* h_rec_1brem[6];
  TCanvas *tc[6];
  TLatex *tt = new TLatex();
  tt->SetTextSize(0.06);
  tt->SetNDC(kTRUE);

  for (int i=5; i>=0; --i) {
    tc[i] = new TCanvas(Form("tc_rad_dep_%d",i),Form("tc_rad_dep_%d",i),1000,1000);
    h_rec_1brem[i] = (TH1F*) f->Get(Form("h_rec_1brem_%d",i));
    h_rec_1brem[i]->Scale(1./h_rec_1brem[i]->GetEntries());

    h_rec_1brem[i]->GetXaxis()->SetTitle(h_rec_1brem[i]->GetXaxis()->GetTitle());
    h_rec_1brem[i]->SetMaximum(max);
    //h_rec_1brem[i]->Draw(i==5?"":"same");
    h_rec_1brem[i]->SetLineColor(1);
    h_rec_1brem[i]->Draw();

    tt->DrawLatex(0.53,0.75,Form("%3.0f < #theta < %3.0f",config[iconfig][1],config[iconfig][2]));
    tt->DrawLatex(0.63,0.65,Form("p_{T}= %4.1f",config[iconfig][0]));
    tt->DrawLatex(0.36,0.85,Form("%d < R_{True}(cm) < %d", i*7, (i+1)*7));
    gPad->SetLogy();
    tc[i]->Print(Form("rad_dep_config%d_r%d.pdf",iconfig,i));
  }

}

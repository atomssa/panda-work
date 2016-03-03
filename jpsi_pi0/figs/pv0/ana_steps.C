void set_style(TH1F* h) {
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

void ana_steps(int ip=0 ){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  TGaxis::SetMaxDigits(3);

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";
  //TFile *fsim,*fbg;
  //if (ip==0){
  //  fsim = TFile::Open(Form("%s/hists/note.aug.2015//anav2_jpsi_brem_plab5.5.root",bdir));
  //  fbg = TFile::Open(Form("%s/hists/note.aug.2015//anav2_pip_pim_brem_plab5.5.root",bdir));
  //} else if (ip==1) {
  //  fsim = TFile::Open(Form("%s/hists/note.aug.2015//anav2_jpsi_brem_plab8.0.root",bdir));
  //  fbg = TFile::Open(Form("%s/hists/note.aug.2015//anav2_pip_pim_brem_plab8.0.root",bdir));
  //} else {
  //  fsim = TFile::Open(Form("%s/hists/note.aug.2015//anav2_jpsi_brem_plab12.0.root",bdir));
  //  fbg = TFile::Open(Form("%s/hists/note.aug.2015//anav2_pip_pim_brem_plab12.0.root",bdir));
  //}
  fsim = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0jpsi_brem4eff_p%d_pass18.root",bdir,ip));
  fbg = TFile::Open(Form("%s/hists/paper.v0.feb.2016/anav2_pi0pipm_brem_p%d_pass18.root",bdir));

  double nevt[3] = {992780, 990142, 991860};
  //double nevt[3] = {32780.0,50142.0,51860.0};
  double f_valid[3] = {24561.88/32780.0,33327.48/50142.0, 23742.89/51860.0};
  for (int ii=0; ii < 3; ++ii) nevt[ii]*=f_valid[ii];

  TLegend *tl = new TLegend(0.35,0.5,0.85,0.9);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  tl->SetTextSize(0.05);

  TH1F* hmep[2][9];
  for (int i=0; i<9; ++i) {
    hmep[0][i] = (TH1F*)fsim->Get(Form("hmep_valid_%d",i));
    set_style(hmep[0][i]);
    hmep[1][i] = (TH1F*)fbg->Get(Form("hmep_%d",i));
    set_style(hmep[1][i]);
  }

  TCanvas *tc_vert = new TCanvas("tc_vert","tc_vert",800,1600);
  tc_vert->Divide(1,2);
  TCanvas *tc_vert = new TCanvas("tc_vert","tc_vert",1600,700);
  tc_vert->Divide(2,1);
  tc_vert->GetPad(1)->SetRightMargin(0.05);
  tc_vert->GetPad(2)->SetRightMargin(0.05);

  int immin = hmep[0][0]->GetXaxis()->FindBin(2.8);
  int immax = hmep[0][0]->GetXaxis()->FindBin(3.3);

  for (int ii=0; ii < 2; ++ii) {

    tc_vert->cd(1+ii);

    hmep[ii][0]->Draw();
    hmep[ii][0]->SetTitle("");
    if (ii==0) hmep[ii][0]->SetTitle("Signal;M_{e^{+}e^{-}}[GeV/c];counts");
    if (ii==1) hmep[ii][0]->SetTitle("Background;M_{#pi^{+}#pi^{-}}[GeV/c];counts");
    hmep[ii][0]->SetLineColor(1);
    if (ii==1) { hmep[ii][0]->Scale(2e-7); }
    if (ii==0) tl->AddEntry(hmep[ii][0],"All pairs (x 2x10^{-7} for bg)","l");
    if (ii==0) cout << "Nnorm = " << hmep[0][0]->Integral(immin,immax)<< endl;

    if (ii==0) cout << "Eff0. Acc = " << 100*hmep[0][0]->Integral(immin,immax)/nevt[ip] << endl;
    hmep[ii][1]->SetLineColor(2);
    hmep[ii][1]->Draw("same");
    //if (ii==0) hmep[ii][1]->Scale(0.9);
    if (ii==0) tl->AddEntry(hmep[ii][1],"Pairs after eid (wt for bg)","l");
    //if (ii==0) cout << "Eff1. Eid = " << 100*hmep[0][1]->Integral(immin,immax)/hmep[0][0]->Integral(immin,immax) << endl;
    //if (ii==0) cout << "Eff2. (All) = " << 100*hmep[0][2]->Integral(immin,immax)/hmep[0][0]->Integral(immin,immax) <<  " should be equal to Eff1 " << endl;
    if (ii==0) cout << "Eff1. Eid = " << 100*hmep[0][1]->Integral(immin,immax)/nevt[ip] << endl;
    if (ii==0) cout << "Eff2. (All) = " << 100*hmep[0][2]->Integral(immin,immax)/nevt[ip] <<  " should be equal to Eff1 " << endl;


    hmep[ii][3]->SetLineColor(4);
    hmep[ii][3]->SetLineWidth(3);
    hmep[ii][3]->Draw("same");
    if (ii==0) tl->AddEntry(hmep[ii][3],"N_{#pi^{0}}>0 after #pi^{0} cuts","l");
    //if (ii==0) cout << "Eff3. (Npi0>0) = " << 100*hmep[0][3]->Integral(immin,immax)/hmep[0][0]->Integral(immin,immax)<< endl;
    if (ii==0) cout << "Eff3. (Npi0>0) = " << 100*hmep[0][3]->Integral(immin,immax)/nevt[ip]<< endl;

    hmep[ii][4]->SetLineColor(3);
    hmep[ii][4]->Draw("same");
    if (ii==0) tl->AddEntry(hmep[ii][4],"Most b-to-b #pi^{0}-(e^{+}e^{-}) pair","l");
    if (ii==1) tl->Draw();

    //if (ii==0) cout << "Eff4. (MostBtoB) = "     << 100*hmep[0][4]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << " should be equal to Eff3" << endl;
    //if (ii==0) cout << "Eff5. (MostBtoB.bis) = " << 100*hmep[0][5]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << " should be equal to Eff4" << endl;
    //if (ii==0) cout << "Eff6. (chi2sig) = "      << 100*hmep[0][6]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << endl;
    //if (ii==0) cout << "Eff7. (chi2bg) = "       << 100*hmep[0][7]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << endl;
    //if (ii==0) cout << "Eff8. (ngam>0) = "       << 100*hmep[0][8]->Integral(immin,immax)/hmep[ii][0]->Integral(immin,immax) << endl;
    if (ii==0) cout << "Eff4. (MostBtoB) = "     << 100*hmep[0][4]->Integral(immin,immax)/nevt[ip] << " should be equal to Eff3" << endl;
    if (ii==0) cout << "Eff5. (MostBtoB.bis) = " << 100*hmep[0][5]->Integral(immin,immax)/nevt[ip] << " should be equal to Eff4" << endl;
    if (ii==0) cout << "Eff6. (chi2sig) = "      << 100*hmep[0][6]->Integral(immin,immax)/nevt[ip] << endl;
    if (ii==0) cout << "Eff7. (chi2bg) = "       << 100*hmep[0][7]->Integral(immin,immax)/nevt[ip] << endl;
    if (ii==0) cout << "Eff8. (ngam>0) = "       << 100*hmep[0][8]->Integral(immin,immax)/nevt[ip] << endl;

  }

  //tc_vert->Print(Form("%s/figs/2015.06.08/ana_steps.pdf",bdir));

}

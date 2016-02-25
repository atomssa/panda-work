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

void pull(){

  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.13);
  //gStyle->SetTitleOffset(0.0,"X");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleAlign(33);
  //TGaxis::SetMaxDigits(3);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111111);

  //TFile *_file0 = TFile::Open("anav2_pi0jpsi_brem4eff_p0_pass17.root");
  TFile *_file0 = TFile::Open("anav2_pi0jpsi_raw_p1_pass17.root");
  //TFile *_file1 = TFile::Open("anav2_pi0pi0jpsi_brem_p0_pass17.root");
  //TFile *_file0 = TFile::Open("anav2_pi0pi0jpsi_brem_p0_pass17.root");

  TF1 *fgaus = new TF1("fgaus","gaus",-2.3,2.3);
  TCanvas *tc1 = new TCanvas("tc1","tc1",1400,1400);
  tc1->Divide(2,2);

  tc1->cd(1);
  TH1F* hpx_pull_ep_r = (TH1F*) _file0->Get("pull/hpx_pull_ep_r");
  hpx_pull_ep_r->SetTitle("#Delta_{x} = (p_{x,REC}-p_{x,MC})/#sigma_{xx}; #Delta_{x}");
  //hpx_pull_ep_r->Sumw2();
  set_style(hpx_pull_ep_r);
  hpx_pull_ep_r->Fit("fgaus","r");

  tc1->cd(2);
  TH1F* hpy_pull_ep_r = (TH1F*) _file0->Get("pull/hpy_pull_ep_r");
  hpy_pull_ep_r->SetTitle("#Delta_{y} = (p_{y,REC}-p_{y,MC})/#sigma_{yy}; #Delta_{y}");
  //hpy_pull_ep_r->Draw();
  set_style(hpy_pull_ep_r);
  hpy_pull_ep_r->Fit("fgaus","r");

  tc1->cd(3);
  TH1F* hpz_pull_ep_r = (TH1F*) _file0->Get("pull/hpz_pull_ep_r");
  hpz_pull_ep_r->SetTitle("#Delta_{z} = (p_{z,REC}-p_{y,MC})/#sigma_{zz}; #Delta_{z}");
  //hpz_pull_ep_r->Draw();
  set_style(hpz_pull_ep_r);
  hpz_pull_ep_r->Fit("fgaus","r");

  TCanvas *tc2 = new TCanvas("tc2","tc2",1400,1400);
  tc2->cd();

  TH1F* pullpx = (TH1F*) _file0->Get("hpi0jpsi_pull4c");
  pullpx->SetTitle("p_{x} pull of kin. fit; (p_{x,NEW}-p_{x,OLD})/#sigma_{xx}");
  pullpx->Rebin(2);
  //pullpx->Draw();
  pullpx->Fit("fgaus","r");
}

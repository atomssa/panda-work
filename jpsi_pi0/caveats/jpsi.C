void set_style(TH1* h, int col) {
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelSize(0.05);

  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelSize(0.05);

  h->SetTitleFont(22,"t");
  h->SetTitleSize(0.08,"t");
  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void jpsi(int itype,  int id = 1) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  const char *tt[2] = {"pi0jpsi_oct14", "pi0pippim_oct14"};
  TCanvas *tc = new TCanvas("tc_jpsi","tc_jpsi",1600,900);
  tc->Divide(2,1);
  for (int itype = 0; itype < 2; ++itype) {
    TFile *f_raw = new TFile(Form("../grid.out/%s/anav2_plab5.5_raw_merged.root",tt[itype]));
    TFile *f_brem = new TFile(Form("../grid.out/%s/anav2_plab5.5_brem_merged.root",tt[itype]));

    TH1F *m_epm_raw_all, *m_epm_brem_all;
    TH1F *m_epm_raw_eid, *m_epm_brem_eid;

    TString title = itype==0? "J/#psi inv. mass;M_{e^{+}e^{-}}":"#pi^{+}#pi^{-} inv. mass;M_{#pi^{+}#pi^{-}}";

    //m_epm_raw_all = (TH1F*)f_raw->Get("hmep_0");
    //m_epm_raw_all->SetTitle(title);
    //m_epm_raw_all->SetLineWidth(2);
    //m_epm_raw_all->SetLineColor(2);
    m_epm_raw_eid = (TH1F*)f_raw->Get(Form("hmep_%d",id));
    m_epm_raw_eid->SetTitle(title);
    set_style(m_epm_raw_eid,2);
    //m_epm_brem_all = (TH1F*)f_brem->Get("hmep_0");
    //m_epm_brem_all->SetTitle(title);
    //m_epm_brem_all->SetLineWidth(2);
    //m_epm_brem_all->SetLineColor(4);
    m_epm_brem_eid = (TH1F*)f_brem->Get(Form("hmep_%d",id));
    m_epm_brem_eid->SetTitle(title);
    m_epm_brem_eid->SetLineWidth(2);
    m_epm_brem_eid->SetLineColor(4);
    set_style(m_epm_brem_eid,4);

    tc->cd(1+itype);
    m_epm_brem_eid->Draw();
    m_epm_raw_eid->Draw("same");

  }

  TLegend *tl =new TLegend(0.45,0.6,0.75,0.8);
  tl->AddEntry(m_epm_raw_eid,"Raw (eid)","pl");
  tl->AddEntry(m_epm_brem_eid,"Corrected (eid)","pl");
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);
  tl->SetTextSize(0.06);
  tl->Draw();
  tc->Print("~/grive/physics/PANDA/collab_meeting/2015.03.17/figs/full_event_sim.pdf");
}

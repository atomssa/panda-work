{

  gStyle->SetOptStat(0);

  TFile *_file0 = TFile::Open("anav2_pi0jpsi_brem4eff_p0_pass17.root");
  //TFile *_file0 = TFile::Open("anav2_pi0pi0jpsi_brem_p0_pass17.root");

  bool dth = true;
  if (dth) {
    hpi0jpsi_chi24c_vs_cm_dth_r->ProjectionY("pass",1,4);
    hpi0jpsi_chi24c_vs_cm_dth_r->ProjectionY("ratio",1,4);
    hpi0jpsi_chi24c_vs_cm_dth_r->ProjectionY("all",1,2000);
  } else {
    hpi0jpsi_chi24c_vs_cm_dph_r->ProjectionY("pass",1,4);
    hpi0jpsi_chi24c_vs_cm_dph_r->ProjectionY("ratio",1,4);
    hpi0jpsi_chi24c_vs_cm_dph_r->ProjectionY("all",1,2000);
  }

  TH1F* all = (TH1F*) gDirectory->Get("all");
  all->SetTitle(Form("blue=all, red=pass #chi^{2} cut;#Delta#%s_{cm}", (dth?"theta":"phi")));
  all->SetLineWidth(2);

  TH1F* pass = (TH1F*) gDirectory->Get("pass");
  pass->SetLineColor(2);
  pass->SetLineWidth(2);

  TH1F* ratio = (TH1F*) gDirectory->Get("ratio");
  ratio->Divide(all);
  ratio->SetTitle(Form("pass/all;#Delta#theta_{cm}",(dth?"theta":"phi")));
  ratio->SetLineColor(1);
  ratio->SetLineWidth(2);

  TCanvas *tc0 = new TCanvas("tc0","tc0");
  tc0->cd();
  ratio->Draw();

  TCanvas *tc1 = new TCanvas("tc1","tc1");
  all->Draw();
  pass->Draw("same");

  TLatex *tt = new TLatex();
  tt->SetNDC(true);
  tt->DrawLatex(0.15,0.6,Form("pass/all = %4.1f%%",100.0*double(pass->GetEntries())/double(all->GetEntries())));

}

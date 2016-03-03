{

  //anav2_pi0jpsi_brem4eff_p0_pass18.root
  gStyle->SetOptStat(0);
  epeff->cd();
  TCanvas *tc =new TCanvas("tc","tc");
  tc->cd();
  hepcosth_jpsi_mc_all->Draw();
  hepcosth_jpsi_mc_all->SetTitle(";cos(#theta_{e^{+}});counts");

  TCanvas *tc2 =new TCanvas("tc2","tc2");
  tc2->cd();

  hepcosth_jpsi_mc_all_wt0->SetTitle(";cos(#theta_{e^{+}});counts");
  hepcosth_jpsi_mc_all_wt0->Draw();
  hepcosth_jpsi_mc_all_wt1->Draw("same");
  hepcosth_jpsi_mc_all_wt1->SetLineColor(2);
  hepcosth_jpsi_mc_all_wt0->SetLineColor(1);
  TLegend *tl = new TLegend(0.3,0.7,0.7,0.88);
  //TLegend *tl = new TLegend(0.3,0.7,0.7,0.9);
  tl->AddEntry(hepcosth_jpsi_mc_all_wt0,"1+cos^{2}(#theta_{e^{+}})","l");
  tl->AddEntry(hepcosth_jpsi_mc_all_wt1,"1+0.4cos^{2}(#theta_{e^{+}})","l");
  tl->Draw();

}

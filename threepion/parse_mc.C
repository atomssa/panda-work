#include "Riostream.h"
#include <string>
#include <iostream>

void parse_mc(int model = 0 /*0=>Dpm, 1=>FTF*/ ) {

  gStyle->SetOptStat(0);

  string fname;
  string tag;
  if (model == 0) {
    fname="DPM_PipPimPi0.dat";
    tag="DPM";
  } else {
    fname="FTF_PipPimPi0.dat";
    tag="FTF";
  }

  ifstream inf;
  inf.open(fname.c_str());

  int evtnum, npart;
  string dummy_s;
  int partnum, pdgid, dummy_i;
  double px,py,pz,E;
  double t,x,y,z;

  double m_prot= 0.938;
  double p_antip = 5.513;
  double E_antip = TMath::Hypot(m_prot, p_antip);
  double beta_cm = p_antip/(E_antip + m_prot);
  cout << "betac_cm = " << beta_cm << endl;

  TVector3 boost_vector(0,0,-beta_cm);
  boost_vector.Print();

  //TH1F* dthe_cm = new TH1F("dthe_cm","dthe_cm",200,179.99,180.01);
  TH1F* dthe_cm = new TH1F(Form("dthe_cm_%s",tag.c_str()),Form("dthe_cm_%s",tag.c_str()),200,0,200);
  TH1F* dthe_lab = new TH1F(Form("dthe_lab_%s",tag.c_str()),Form("dthe_lab_%s",tag.c_str()),200,0,185);
  TH1F* dphi_cm = new TH1F(Form("dphi_cm_%s",tag.c_str()),Form("dphi_cm_%s",tag.c_str()),200,179.99,180.01);

  TH1F* the_pi0_cm = new TH1F(Form("the_pi0_cm_%s",tag.c_str()),Form("%s #theta_{#pi^{0}}^{cm};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pi0_cm->SetLineColor(2);
  the_pi0_cm->SetLineWidth(2);
  TH1F* the_pippim_cm = new TH1F(Form("the_pippim_cm_%s",tag.c_str()),Form("%s #theta_{#pi^{+}#pi^{-}}^{cm};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pippim_cm->SetLineColor(4);
  the_pippim_cm->SetLineWidth(2);
  TH1F* the_pi0_lab = new TH1F(Form("the_pi0_lab_%s",tag.c_str()),Form("%s #theta_{#pi^{0}}^{lab};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pi0_lab->SetLineColor(2);
  the_pi0_lab->SetLineWidth(2);
  TH1F* the_pippim_lab = new TH1F(Form("the_pippim_lab_%s",tag.c_str()),Form("%s #theta_{#pi^{+}#pi^{-}}^{lab};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pippim_lab->SetLineColor(4);
  the_pippim_lab->SetLineWidth(2);
  TH1F* phi_pi0_cm = new TH1F(Form("phi_pi0_cm_%s",tag.c_str()),Form("%s #phi_{#pi^{0}}^{cm};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pi0_cm->SetLineColor(2);
  phi_pi0_cm->SetLineWidth(2);
  TH1F* phi_pippim_cm = new TH1F(Form("phi_pippim_cm_%s",tag.c_str()),Form("%s #phi_{#pi^{+}#pi^{-}}^{cm};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pippim_cm->SetLineColor(4);
  phi_pippim_cm->SetLineWidth(2);
  TH1F* phi_pi0_lab = new TH1F(Form("phi_pi0_lab_%s",tag.c_str()),Form("%s #phi_{#pi^{0}}^{lab};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pi0_lab->SetLineColor(2);
  phi_pi0_lab->SetLineWidth(2);
  TH1F* phi_pippim_lab = new TH1F(Form("phi_pippim_lab_%s",tag.c_str()),Form("%s #phi_{#pi^{+}#pi^{-}}^{lab};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pippim_lab->SetLineColor(4);
  phi_pippim_lab->SetLineWidth(2);

  TH1F* minv_pippim_lab = new TH1F(Form("minv_pippim_lab_%s",tag.c_str()),Form("%s minv_pippim_lab; M_{inv} [GeV/c^2]; dN/dM_{inv}",tag.c_str()),200,0,4);
  minv_pippim_lab->SetLineColor(2);
  minv_pippim_lab->SetLineWidth(2);
  TH1F* minv_pippim_cm = new TH1F(Form("minv_pippim_cm_%s",tag.c_str()),Form("%s minv_pippim_cm; M_{inv} [GeV/c^2]; dN/dM_{inv}",tag.c_str()),200,0,4);
  minv_pippim_cm->SetLineColor(4);
  minv_pippim_cm->SetLineWidth(2);

  int in2Sig=0, tot=0;

  while(1) {
    inf >> evtnum >> npart;
    //if (evtnum>=10) break;
    if (npart!=3) {
      cout << "Something wrong, npart= " << npart << " != 3" << endl;
      return;
    }
    if (!inf.good()) break;
    if (evtnum%1000==0) cout << " ================  EventNo= " << evtnum << " ================= " << endl;

    for (int i=0; i<15; ++i) inf>>dummy_s;
    TLorentzVector mom_pi0,mom_cm_pi0;
    TLorentzVector mom_pip,mom_cm_pip;
    TLorentzVector mom_pim,mom_cm_pim;
    TLorentzVector mom_pich,mom_cm_pich;

    for (int ipart=0; ipart<3; ++ipart) {
      inf >> partnum >> pdgid;
      for (int j=0; j<5; ++j) inf>>dummy_i;
      inf>>px>>py>>pz>>E;
      inf>>t>>x>>y;
      if (model==1) inf >>z; // (for some reason z is not tabulate...)
      if (pdgid==111) { mom_pi0.SetPxPyPzE(px,py,pz,E); mom_cm_pi0.SetPxPyPzE(px,py,pz,E); mom_cm_pi0.Boost(boost_vector); }
      if (pdgid==211) { mom_pip.SetPxPyPzE(px,py,pz,E); mom_cm_pip.SetPxPyPzE(px,py,pz,E); mom_cm_pip.Boost(boost_vector); }
      if (pdgid==-211) { mom_pim.SetPxPyPzE(px,py,pz,E); mom_cm_pim.SetPxPyPzE(px,py,pz,E); mom_cm_pim.Boost(boost_vector); }
    }
    mom_pich = mom_pip + mom_pim;
    mom_cm_pich = mom_cm_pip + mom_cm_pim;

    TLorentzVector mom_tot_cm = mom_cm_pi0 + mom_cm_pip + mom_cm_pim;
    TLorentzVector mom_tot = mom_pi0 + mom_pip + mom_pim;

    double rtd= TMath::RadToDeg();

    //mom_cm_pip.Print();

    //cout << "Mpi0 = " << mom_pi0.M() << " Mpip= " << mom_pip.M() << " Mpim= " << mom_pim.M() << endl;
    //cout << "Minv = " << mom_tot.M() << "  Minv(cm) = " << mom_tot_cm.M() << endl;
    //cout << "Minv(pi+pi-)= " << mom_pich.M() << " Minv(cm,pi+pi-)= " << mom_pich.M() << endl;
    //cout << "CM Total Mom = "; mom_tot_cm.Print();
    //cout << "LAB Tot Mom= "; mom_tot.Print();
    //
    //cout << "THETA LAB: pi0= " << rtd*mom_pi0.Theta() << " pip= " << rtd*mom_pip.Theta() << " pim= " << rtd*mom_pim.Theta()
    //	 << " pich= " << rtd*mom_pich.Theta() << " D0-ch= " << rtd*(mom_pi0.Theta()+mom_pich.Theta()) <<endl;
    //
    //cout << "THETA  CM: pi0= " << rtd*mom_cm_pi0.Theta() << " pip= " << rtd*mom_cm_pip.Theta() << " pim= " << rtd*mom_cm_pim.Theta()
    //	 << " pich= " << rtd*mom_cm_pich.Theta() << " D0-ch= " << rtd*(mom_cm_pi0.Theta()+mom_cm_pich.Theta()) << endl;
    //
    //cout << "PHI   LAB: pi0= " << rtd*mom_pi0.Phi() << " pip= " << rtd*mom_pip.Phi() << " pim= " << rtd*mom_pim.Phi()
    //	 << " pich= " << rtd*mom_pich.Phi() << " D0-ch= " << rtd*(mom_pi0.Phi()-mom_pich.Phi()) << endl;
    //
    //cout << "PHI    CM: pi0= " << rtd*mom_cm_pi0.Phi() << " pip= " << rtd*mom_cm_pip.Phi() << " pim= " << rtd*mom_cm_pim.Phi()
    //	 << " pich= " << rtd*mom_cm_pich.Phi() << " D0-ch= " << rtd*(mom_cm_pi0.Phi()-mom_cm_pich.Phi()) << endl;

    dthe_cm->Fill( rtd*(mom_cm_pi0.Theta()+mom_cm_pich.Theta()) );
    dthe_lab->Fill( rtd*(mom_pi0.Theta()+mom_pich.Theta()) );
    dphi_cm->Fill( rtd*(mom_cm_pi0.Phi()-mom_cm_pich.Phi()) );

    the_pi0_cm->Fill( rtd*mom_cm_pi0.Theta() );
    the_pippim_cm->Fill( rtd*mom_cm_pich.Theta() );
    the_pi0_lab->Fill( rtd*mom_pi0.Theta() );
    the_pippim_lab->Fill( rtd*mom_pich.Theta() );
    phi_pi0_cm->Fill( rtd*mom_cm_pi0.Phi() );
    phi_pippim_cm->Fill( rtd*mom_cm_pich.Phi() );
    phi_pi0_lab->Fill( rtd*mom_pi0.Phi() );
    phi_pippim_lab->Fill( rtd*mom_pich.Phi() );

    minv_pippim_cm->Fill( mom_cm_pich.M() );
    minv_pippim_lab->Fill( mom_pich.M() );

    tot++;
    if ( 2.96 < mom_pich.M() && mom_pich.M() < 3.22 ) in2Sig++;

  }

  double frac = (double)in2Sig/(double)tot;
  cout << "Frac ( %): " << frac*100 << endl;

  return;

  TCanvas *c0 = new TCanvas("c0","c0");
  dphi_cm->Draw();
  TCanvas *c1 = new TCanvas("c1","c1");
  dthe_cm->Draw();
  //TCanvas *c2 = new TCanvas("c2","c2");
  //dthe_lab->Draw();
  TCanvas *c_minv = new TCanvas("cminv","cminv");
  c_minv->cd();
  minv_pippim_cm->Draw();
  minv_pippim_lab->Draw("same");

  TLegend *tl = new TLegend(0.35,0.65,0.65,0.89);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  TCanvas *c3 = new TCanvas("c3","c3",1400,1000);
  c3->Divide(2,2);

  c3->cd(1);
  the_pi0_cm->SetTitle(Form("%s, %s", the_pi0_cm->GetTitle(), the_pippim_cm->GetTitle() ));
  the_pi0_cm->Draw();
  tl->AddEntry(the_pi0_cm,Form("#pi^{0}, %s",tag.c_str()),"l");
  the_pippim_cm->Draw("same");
  tl->AddEntry(the_pippim_cm,Form("#pi^{+}#pi^{-}, %s",tag.c_str()),"l");
  tl->Draw();

  c3->cd(2);
  gPad->SetLogy();
  the_pippim_lab->SetTitle(Form("%s, %s", the_pi0_lab->GetTitle(), the_pippim_lab->GetTitle() ));
  the_pippim_lab->Draw();
  the_pi0_lab->Draw("same");
  tl->Draw();

  c3->cd(3);
  phi_pi0_cm->SetTitle(Form("%s, %s", phi_pi0_cm->GetTitle(), phi_pippim_cm->GetTitle() ));
  phi_pi0_cm->SetMinimum(0.0);
  phi_pi0_cm->Draw();
  phi_pippim_cm->Draw("same");
  tl->Draw();

  c3->cd(4);
  phi_pi0_lab->SetTitle(Form("%s, %s", phi_pi0_lab->GetTitle(), phi_pippim_lab->GetTitle() ));
  phi_pi0_lab->SetMinimum(0.0);
  phi_pi0_lab->Draw();
  phi_pippim_lab->Draw("same");
  tl->Draw();

  c3->Print(Form("ang_dist_pbarp_pi0pippim_%s.png",tag.c_str()));
  //if (model==0) {
  //  c3->Print("ang_dist_pbarp_pi0pippim_dpm.png");
  //} else {
  //  c3->Print("ang_dist_pbarp_pi0pippim_ftf.png");
  //}


  //TFile *fout = TFile::Open(Form("output_%s.root",fname.c_str()),"RECREATE");
  //fout->cd();
  //dthe_cm->Write();
  //dthe_lab->Write();
  //dphi_cm->Write();
  //the_pi0_cm->Write();
  //the_pippim_cm->Write();
  //the_pi0_lab->Write();
  //the_pippim_lab->Write();
  //phi_pi0_cm->Write();
  //phi_pippim_cm->Write();
  //phi_pi0_lab->Write();
  //phi_pippim_lab->Write();
  //minv_pippim_lab->Write();
  //minv_pippim_cm->Write();





}

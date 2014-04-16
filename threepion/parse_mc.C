#include "Riostream.h"
#include <string>
#include <iostream>

void parse_mc(int model = 0 /*0=>Dpm, 1=>FTF*/ ) {

  gStyle->SetOptStat(0);
  
  string fname;
  if (model == 0) {
    fname="Dpm_PipPimPi0.dat";
  } else {
    fname="FTF_PipPimPi0.dat";
  }
  
  ifstream inf;
  inf.open(fname.c_str());

  int evtnum, npart;
  string dummy_s;
  int partnum, pdgid, dummy_i;
  double px,py,pz,E;
  double t,x,y,z;
  TLorentzVector mom_pbar;

  //mom_pbar.SetXYZM(0.0,0.0,3.5/2.0,0.938);
  if (model==0) {
    mom_pbar.SetXYZM(0.0,0.0,1.477425,0.938);   // Seems to work for Dpm file (trial and error)
  } else {
    mom_pbar.SetXYZM(0.0,0.0,1.4771898,0.938); // Seems to work for FTF file (trial and error)
  }

  TVector3 boost_vector = -((mom_pbar).BoostVector());

  TH1F* dthe_cm = new TH1F("dthe_cm","dthe_cm",200,179.99,180.01);
  //TH1F* dthe_cm = new TH1F("dthe_cm","dthe_cm",200,0,200);
  TH1F* dthe_lab = new TH1F("dthe_lab","dthe_lab",200,0,185);
  TH1F* dphi_cm = new TH1F("dphi_cm","dphi_cm",200,179.99,180.01);  

  TH1F* the_pi0_cm = new TH1F("the_pi0_cm","#theta_{#pi^{0}}^{cm};#theta[deg];dN/d#theta",200,0,180);
  the_pi0_cm->SetLineColor(2);
  the_pi0_cm->SetLineWidth(2);
  TH1F* the_pippim_cm = new TH1F("the_pippim_cm","#theta_{#pi^{+}#pi^{-}}^{cm};#theta[deg];dN/d#theta",200,0,180);
  the_pippim_cm->SetLineColor(4);
  the_pippim_cm->SetLineWidth(2);
  TH1F* the_pi0_lab = new TH1F("the_pi0_lab","#theta_{#pi^{0}}^{lab};#theta[deg];dN/d#theta",200,0,180);
  the_pi0_lab->SetLineColor(2);
  the_pi0_lab->SetLineWidth(2);
  TH1F* the_pippim_lab = new TH1F("the_pippim_lab","#theta_{#pi^{+}#pi^{-}}^{lab};#theta[deg];dN/d#theta",200,0,180);
  the_pippim_lab->SetLineColor(4);
  the_pippim_lab->SetLineWidth(2);
  TH1F* phi_pi0_cm = new TH1F("phi_pi0_cm","#phi_{#pi^{0}}^{cm};#phi[deg];dN/d#phi",200,0,180);
  phi_pi0_cm->SetLineColor(2);
  phi_pi0_cm->SetLineWidth(2);
  TH1F* phi_pippim_cm = new TH1F("phi_pippim_cm","#phi_{#pi^{+}#pi^{-}}^{cm};#phi[deg];dN/d#phi",200,0,180);
  phi_pippim_cm->SetLineColor(4);
  phi_pippim_cm->SetLineWidth(2);
  TH1F* phi_pi0_lab = new TH1F("phi_pi0_lab","#phi_{#pi^{0}}^{lab};#phi[deg];dN/d#phi",200,0,180);
  phi_pi0_lab->SetLineColor(2);
  phi_pi0_lab->SetLineWidth(2);
  TH1F* phi_pippim_lab = new TH1F("phi_pippim_lab","#phi_{#pi^{+}#pi^{-}}^{lab};#phi[deg];dN/d#phi",200,0,180);
  phi_pippim_lab->SetLineColor(4);
  phi_pippim_lab->SetLineWidth(2);
  
  while(1) {
    inf >> evtnum >> npart;
    //if (evtnum>=10) break;
    if (npart!=3) {
      cout << "Something wrong, npart= " << npart << " != 3" << endl;
      return;
    }
    if (!inf.good()) break;
    cout << " ================  EventNo= " << evtnum << " ================= " << endl;

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
    
    mom_cm_pip.Print();

    cout << "Mpi0 = " << mom_pi0.M() << " Mpip= " << mom_pip.M() << " Mpim= " << mom_pim.M() << endl;
    cout << "Minv = " << mom_tot.M() << "  Minv(cm) = " << mom_tot_cm.M() << endl;
    cout << "Minv(pi+pi-)= " << mom_pich.M() << " Minv(cm,pi+pi-)= " << mom_pich.M() << endl;
    cout << "CM Total Mom = "; mom_tot_cm.Print();
    cout << "LAB Tot Mom= "; mom_tot.Print();

    cout << "THETA LAB: pi0= " << rtd*mom_pi0.Theta() << " pip= " << rtd*mom_pip.Theta() << " pim= " << rtd*mom_pim.Theta()
	 << " pich= " << rtd*mom_pich.Theta() << " D0-ch= " << rtd*(mom_pi0.Theta()+mom_pich.Theta()) <<endl;

    cout << "THETA  CM: pi0= " << rtd*mom_cm_pi0.Theta() << " pip= " << rtd*mom_cm_pip.Theta() << " pim= " << rtd*mom_cm_pim.Theta()
	 << " pich= " << rtd*mom_cm_pich.Theta() << " D0-ch= " << rtd*(mom_cm_pi0.Theta()+mom_cm_pich.Theta()) << endl;

    cout << "PHI   LAB: pi0= " << rtd*mom_pi0.Phi() << " pip= " << rtd*mom_pip.Phi() << " pim= " << rtd*mom_pim.Phi()
	 << " pich= " << rtd*mom_pich.Phi() << " D0-ch= " << rtd*(mom_pi0.Phi()-mom_pich.Phi()) << endl;

    cout << "PHI    CM: pi0= " << rtd*mom_cm_pi0.Phi() << " pip= " << rtd*mom_cm_pip.Phi() << " pim= " << rtd*mom_cm_pim.Phi()
	 << " pich= " << rtd*mom_cm_pich.Phi() << " D0-ch= " << rtd*(mom_cm_pi0.Phi()-mom_cm_pich.Phi()) << endl;


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
    
  }

  TCanvas *c0 = new TCanvas("c0","c0");
  dphi_cm->Draw();  
  TCanvas *c1 = new TCanvas("c1","c1");
  dthe_cm->Draw();
  //TCanvas *c2 = new TCanvas("c2","c2");  
  //dthe_lab->Draw();

  TLegend *tl = new TLegend(0.35,0.65,0.65,0.89);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  TCanvas *c3 = new TCanvas("c3","c3",1400,1000);
  c3->Divide(2,2);
  
  c3->cd(1);
  the_pi0_cm->SetTitle(Form("%s, %s", the_pi0_cm->GetTitle(), the_pippim_cm->GetTitle() ));
  the_pi0_cm->Draw();
  tl->AddEntry(the_pi0_cm,Form("#pi^{0}, %s",(model==0?"DPM":"FTF")),"l");
  the_pippim_cm->Draw("same");
  tl->AddEntry(the_pippim_cm,Form("#pi^{+}#pi^{-}, %s",(model==0?"DPM":"FTF")),"l");  
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

  if (model==0) {
    c3->Print("ang_dist_pbarp_pi0pippim_dpm.png");
  } else {
    c3->Print("ang_dist_pbarp_pi0pippim_ftf.png");    
  }
  
  //  ================  EventNo= 9999 ================= (boost = pbarmom/2 )
  //(x,y,z,t)=(0.217625,0.090804,12.220965,12.224041) (P,eta,phi,E)=(12.223240,4.641125,0.395290,12.224041)
  //Mpi0 = 0.134003 Mpip= 0.139953 Mpim= 0.139905
  //Minv = 3.49916  Minv(cm) = 3.49916
  //Minv(pi+pi-)= 1.31048 Minv(cm,pi+pi-)= 1.31048
  //THETA LAB: pi0= 90.8613 pip= 4.40236 pim= 21.9122 pich= 8.10117
  //THETA  CM: pi0= 28.2517 pip= 1.10541 pim= 5.56212 pich= 2.01164
  //PHI   LAB: pi0= 57.1132 pip= 22.6484 pim= -130.636 pich= -122.887
  //PHI    CM: pi0= 57.1132 pip= 22.6484 pim= -130.636 pich= -122.887

  // ================  EventNo= 9999 ================= (boost = pbarmom)
  //(x,y,z,t)=(0.217625,0.090804,23.307016,23.308629) (P,eta,phi,E)=(23.308208,5.286659,0.395290,23.308629)
  //Mpi0 = 0.134003 Mpip= 0.139953 Mpim= 0.139905
  //Minv = 3.49916  Minv(cm) = 3.49916
  //Minv(pi+pi-)= 1.31048 Minv(cm,pi+pi-)= 1.31048
  //THETA LAB: pi0= 90.8613 pip= 4.40236 pim= 21.9122 pich= 8.10117
  //THETA  CM: pi0= 15.0175 pip= 0.579672 pim= 2.91831 pich= 1.05432
  //PHI   LAB: pi0= 57.1132 pip= 22.6484 pim= -130.636 pich= -122.887
  //PHI    CM: pi0= 57.1132 pip= 22.6484 pim= -130.636 pich= -122.887

		       
}

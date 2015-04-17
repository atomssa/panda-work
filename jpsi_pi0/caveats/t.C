#define t_cxx
#include "t.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
//#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>

#include <iostream>

using namespace std;
void t::Loop()
{

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  int nbin = 200;
  double xmin = -0.2, xmax = 0.2;
  TH2F* h_rad_calc_vs_true = new TH2F("h_rad_calc_vs_true","Calculated radius vs. true radius; R_{True}; R_{Calc}", 200, 0, 40, 200, 0, 40);
  set_style(h_rad_calc_vs_true);
  TH2F* h_zed_calc_vs_true = new TH2F("h_zed_calc_vs_true","Calculated Z vs. true Z; Z_{True}; Z_{Calc}", 500, 0, 500, 500, 0, 500);
  set_style(h_zed_calc_vs_true);
  double zdet[7] = {0,6,8.5,12,17,25,500};
  TH1F* h_zed_calc[6];
  for (int ii = 0; ii < 6; ++ii) {
    h_zed_calc[ii] = new TH1F(Form("h_zed_calc_%d",ii),Form("h_zed_calc_%d",ii),200,0,60);
  }
  TH2F* h_cor_vs_rec = new TH2F("h_cor_vs_rec","Reconstructed vs. corrected momentum",nbin,xmin,xmax,nbin,xmin,xmax);
  TH2F* prec_vs_rec = new TH2F("prec_vs_rec","prec_vs_rec",nbin,-0.3,1.0,nbin,0,10.0);
  TH1F* h_rec = th1f("h_rec","Resolution of reconstructed momentum;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,1);
  TH1F* h_wtd = th1f("h_wtd","Resolution of fully corrected (with weight) momentum;(p_{MC} - p_{WTD})/p_{MC}",nbin,xmin,xmax,9);
  TH1F* h_cor = th1f("h_cor","Resolution of fully corrected momentum;(p_{MC} - p_{COR})/p_{MC}",nbin,xmin,xmax,2);
  TH1F* h_sep = th1f("h_sep","Resolution of sep. clust corrected momentum;(p_{MC} - p_{SEP})/p_{MC}",nbin,xmin,xmax,3);
  TH1F* h_mrg = th1f("h_mrg","Resolution of mrg. clust momentum;(p_{MC} - p_{mrg})/p_{MC}",nbin,xmin,xmax,4);
  TH1F* h_sto = th1f("h_sto","Resolution of stored momentum;(p_{MC} - p_{STO})/p_{MC}",nbin,xmin,xmax,7);

  TH1F* h_out = th1f("h_out","Resolution W.R.T p_{MC}-E_{#gamma} #approx p_{OUT} ... (Single Brem < 42cm);((p_{MC}-E_{#gamma})-p_{KF})/(p_{MC}-E_{#gamma})",nbin,xmin,xmax,4);

  TH1F* h_nmcb_gt1mev = th1f("h_nmcb_gt1mev","Percentage of tracks with N Brem. #gamma (E>1MeV) R<42cm;N;%-age",10,0,10,1);
  TH1F* h_dphi_all = th1f("h_dphi_all","#Delta#phi (all bumps)",nbin,-60,40,4);
  TH1F* h_dthe_all = th1f("h_dthe_all","#Delta#theta (all bumps)",nbin,-20,20,4);
  TH2F* h_dphi_dthe_all = new TH2F("h_dphi_dthe_all","#Delta#phi vs #Delta#theta (all bumps)",nbin,-20,20,nbin,-60,40);
  TH1F* h_dphi_brem = th1f("h_dphi_brem","#Delta#phi (brem bumps)",nbin,-60,40,2);
  TH1F* h_dthe_brem = th1f("h_dthe_brem","#Delta#theta (brem bumps)",nbin,-20,20,2);
  TH2F* h_dphi_dthe_brem = new TH2F("h_dphi_dthe_brem","#Delta#phi vs #Delta#theta (brem bumps)",nbin,-20,20,nbin,-60,40);
  TH1F* h_dphi_brem_cut = th1f("h_dphi_brem_cut","#Delta#phi (brem bumps)",nbin,-60,40,2);
  TH1F* h_dthe_brem_cut = th1f("h_dthe_brem_cut","#Delta#theta (brem bumps)",nbin,-20,20,2);
  TH2F* h_dphi_dthe_brem_cut = new TH2F("h_dphi_dthe_brem_cut","#Delta#phi vs #Delta#theta (brem bumps)",nbin,-20,20,nbin,-60,40);
  TH1F* h_rec_1brem[6];
  TH1F* h_rec_1brem_found[6];
  int col_rec_1brem[6] = {1,2,4,6,7,9};
  for (int ii = 0; ii < 6; ++ii) {
    h_rec_1brem_found[ii] = th1f(Form("h_rec_1brem_found_%d",ii),"Resolution of reconstructed momentum;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
    h_rec_1brem[ii] = th1f(Form("h_rec_1brem_%d",ii),"Resolution of reconstructed momentum;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
  }

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (ientry%10000==0) cout << "Loop i: " << ientry << "/" << nentries << endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //cout << "----------- New Event " << jentry << " nch= " << nch << " nsb= " << nsb << " nmcb= " << nmcb << " --------------- " << endl;
    if (nch>1) continue;
    for (int ich = 0; ich < nch; ++ich) {
      //cout << ">>> track " << ich << " <<<< " << endl;
      //cout << " q= " << charge[ich] << " pmc= " << mom_mc[ich] << " pr= " << mom_rec[ich] <<  " pc= " << mom_cor[ich] << endl;
      //cout << " _nsb= " << _nsb[ich] << " _nmcb= " << _nmcb[ich] << endl;

      if ( !is_prim[ich] ) {
      	cout << "Non Primary track  in event " << jentry << endl;
      	continue;
      }

      //if ( the_mc[ich] > 15 ) continue;

      float res_cor = (mom_mc[ich]-mom_cor[ich])/mom_mc[ich];
      float res_wtd = (mom_mc[ich]-mom_wcor[ich])/mom_mc[ich];
      float res_rec = (mom_mc[ich]-mom_rec[ich])/mom_mc[ich];
      float res_sep = (mom_mc[ich]-mom_sep[ich])/mom_mc[ich];
      float res_mrg = (mom_mc[ich]-mom_mrg[ich])/mom_mc[ich];
      float res_sto = (mom_mc[ich]-mom_stored[ich])/mom_mc[ich];

      h_cor->Fill(res_cor);
      h_wtd->Fill(res_wtd);
      h_sep->Fill(res_sep);
      h_mrg->Fill(res_mrg);
      h_rec->Fill(res_rec);
      h_sto->Fill(res_sto);
      h_cor_vs_rec->Fill(res_rec,res_cor);
      prec_vs_rec->Fill(res_cor,mom_rec[ich]);

      //if ( phi[ich] <23 ) continue;
      //if (nsb==0) continue;
      //if (charge[ich] < 0 ) continue;
      //if (nch != 2) continue;

      // Tracks that emitted only one bremstrahlung photon
      //  and have only one separated photon bump found
      //  and the two have a good match
      if ( _nmcb[ich]==1 && _nsb[ich]==1) {
	if (sb_match[isb_s[ich]] != mcb_match[imcb_s[ich]] ) {
	  cout << "Single phton events non matching mcb and sb. This shoudlnt happen" << endl;
	}
	double rad_mc = mcb_rad[imcb_s[ich]];
	if (rad_mc<42) {
	  int irad_mc = int(rad_mc/7);
	  h_rec_1brem_found[irad_mc]->Fill(res_rec);
	}
      }

      // Tracks that emitted only one bremstrahlung photon
      if ( _nmcb[ich]==1 ) {
	double rad_mc = mcb_rad[imcb_s[ich]];
	if (rad_mc<42) {
	  int irad_mc = int(rad_mc/7);
	  h_rec_1brem[irad_mc]->Fill(res_rec);
	}
	double mom_out = (mom_mc[ich]-mcb_ene[imcb_s[ich]]);
	h_out->Fill((mom_out-mom_rec[ich])/mom_out);
      }

      for (int isb=isb_s[ich]; isb<isb_e[ich]; ++isb) {
      	if ( sb_match[isb] >= 0) {
      	  if (sb_score[isb]< 1) continue;
      	  if (sb_match[isb] >= nmcb ) {
      	    cout << "ient= " << jentry << " ich= " << ich << " sb_match= " << sb_match[isb] << " nmcb= " << nmcb;
      	    cout << "!! sb_match>nmcb. Shouldnt happen."<< endl;
	    continue;
      	  }
      	  //if (nmcb > 4) continue;
      	  h_rad_calc_vs_true->Fill(mcb_rad[sb_match[isb]], sb_rcalc[isb]*42./30.);
      	  //h_rad_calc_vs_true->Fill(mcb_rad[sb_match[isb]], sb_rcalc[isb]);

      	  h_zed_calc_vs_true->Fill(mcb_zed[sb_match[isb]], sb_rcalc[isb]/TMath::Tan(TMath::DegToRad()*sb_the[isb]));
	  int ized = 0;
	  for (int ii=0; ii<6;++ii) {
	    if (mcb_zed[sb_match[isb]]>zdet[ii] && mcb_zed[sb_match[isb]]<zdet[ii+1]) {
	      ized = ii;
	      break;
	    }
	  }
      	  h_zed_calc[ized]->Fill(sb_rcalc[isb]/TMath::Tan(TMath::DegToRad()*sb_the[isb]));
	}
      }

      for (int iab = 0; iab < nab; ++iab) {
	h_dphi_all->Fill(ab_phi[iab]-phi[ich]);
	h_dthe_all->Fill(ab_the[iab]-the[ich]);
	h_dphi_dthe_all->Fill(ab_the[iab]-the[ich],ab_phi[iab]-phi[ich]);
	if ( ab_match[iab] > -1 ) {
	  h_dphi_brem->Fill(ab_phi[iab]-phi[ich]);
	  h_dthe_brem->Fill(ab_the[iab]-the[ich]);
	  h_dphi_dthe_brem->Fill(ab_the[iab]-the[ich],ab_phi[iab]-phi[ich]);
	}
	if ( ab_isb[iab] == ich ) {
	  h_dphi_brem_cut->Fill(ab_phi[iab]-phi[ich]);
	  h_dthe_brem_cut->Fill(ab_the[iab]-the[ich]);
	  h_dphi_dthe_brem_cut->Fill(ab_the[iab]-the[ich],ab_phi[iab]-phi[ich]);
	}
      }

      // Count number of MC brem photons (R<42cm) emitted from this track
      int nmcb_gt1mev = 0;
      for (int imcb=imcb_s[ich]; imcb<imcb_e[ich]; ++imcb) {
	if ( mcb_ene[ich] > 0.001 ) {
	  nmcb_gt1mev++;
	}
      }
      h_nmcb_gt1mev->Fill(nmcb_gt1mev);
    }
  }
  h_nmcb_gt1mev->Scale(100./h_nmcb_gt1mev->GetEntries());

  TFile *fout = TFile::Open(fout_name,"RECREATE");
  fout->cd();
  h_rad_calc_vs_true->Write();
  h_zed_calc_vs_true->Write();
  for (int ii = 0; ii < 6; ++ii) {
    h_zed_calc[ii]->Write();
  }
  h_cor_vs_rec->Write();
  prec_vs_rec->Write();
  h_rec->Write();
  h_wtd->Write();
  h_cor->Write();
  h_sep->Write();
  h_mrg->Write();
  h_sto->Write();
  h_out->Write();
  h_nmcb_gt1mev->Write();
  h_dphi_all->Write();
  h_dthe_all->Write();
  h_dphi_dthe_all->Write();
  h_dphi_brem->Write();
  h_dthe_brem->Write();
  h_dphi_dthe_brem->Write();
  h_dphi_brem_cut->Write();
  h_dthe_brem_cut->Write();
  h_dphi_dthe_brem_cut->Write();
  for (int ii = 0; ii < 6; ++ii) {
    h_rec_1brem[ii]->Write();
    h_rec_1brem_found[ii]->Write();
  }
  fout->Close();

}




//   In a ROOT session, you can do:
//      Root > .L t.C
//      Root > t t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

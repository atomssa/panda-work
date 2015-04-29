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

  //TH1F* h_rec = th1f("h_rec","Resolution of reconstructed momentum;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,1);
  //TH1F* h_wtd = th1f("h_wtd","Resolution of fully corrected (with weight) momentum;(p_{MC} - p_{WTD})/p_{MC}",nbin,xmin,xmax,9);
  //TH1F* h_cor = th1f("h_cor","Resolution of fully corrected momentum;(p_{MC} - p_{COR})/p_{MC}",nbin,xmin,xmax,2);
  //TH1F* h_sep = th1f("h_sep","Resolution of sep. clust corrected momentum;(p_{MC} - p_{SEP})/p_{MC}",nbin,xmin,xmax,3);
  //TH1F* h_mrg = th1f("h_mrg","Resolution of mrg. clust momentum;(p_{MC} - p_{mrg})/p_{MC}",nbin,xmin,xmax,4);
  //TH1F* h_sto = th1f("h_sto","Resolution of stored momentum;(p_{MC} - p_{STO})/p_{MC}",nbin,xmin,xmax,7);

  TH1F* h_rec = th1f("h_rec","(p_{MC} - p_{KF})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_cor = th1f("h_cor","(p_{MC} - p_{COR})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_wcor = th1f("h_wcor","(p_{MC} - p_{WTD})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_mrg = th1f("h_mrg","(p_{MC} - p_{MRG})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_sep = th1f("h_sep","(p_{MC} - p_{SEP})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_stored = th1f("h_stored","(p_{MC} - p_{STO})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_sep_w = th1f("h_sep_w","(p_{MC} - p_{SEP_W})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_sep_w_bf = th1f("h_sep_w_bf","(p_{MC} - p_{SEP_WBF})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_mrg_w = th1f("h_mrg_w","(p_{MC} - p_{MRG_W})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_mrg_w_bf = th1f("h_mrg_w_bf","(p_{MC} - p_{MRG_WBF})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_mrg_pc = th1f("h_mrg_pc","(p_{MC} - p_{MRG_PC})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_mrg_w_pc = th1f("h_mrg_w_pc","(p_{MC} - p_{MRG_WPC})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_mrg_w_bf_pc = th1f("h_mrg_w_bf_pc","(p_{MC} - p_{MRG_WBFPC})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_wbfcor = th1f("h_wbfcor","(p_{MC} - p_{SEP_WBF_MRG})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_wbfcor_mw_bf = th1f("h_wbfcor_mw_bf","(p_{MC} - p_{SEP_WBF_MRG_WBF})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_wbfcor_mw_bf_pc = th1f("h_wbfcor_mw_bf_pc","(p_{MC} - p_{SEP_WBF_MRG_WBFPC})/p_{MC}", nbin,xmin,xmax,1);

  TH1F* h_out = th1f("h_out","Resolution W.R.T p_{MC}-E_{#gamma} #approx p_{OUT} ... (Single Brem < 42cm);((p_{MC}-E_{#gamma})-p_{KF})/(p_{MC}-E_{#gamma})",nbin,xmin,xmax,4);

  TH1F* h_nmcb_gt1mev = th1f("h_nmcb_gt1mev","Percentage of tracks with N Brem. #gamma (E>1MeV) R<42cm;N;%-age",10,0,10,1);

  const double dth_min = -20;
  const double dth_max = 20;
  const double dphi_min = -100;
  const double dphi_max = 20;

  // dphi_dtheta for all bumps
  TH2F* h_dphi_dthe_all = new TH2F("h_dphi_dthe_all","#Delta#phi vs #Delta#theta (all bumps)",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);
  // dphi vs. dtheta for MC brem photons
  TH2F* h_dphi_dthe_brem = new TH2F("h_dphi_dthe_brem","#Delta#phi vs #Delta#theta (brem bumps)",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);
  // dphi vs. dtheta for bumps that pass acceptance (filled from ab lists)
  TH2F* h_dphi_dthe_brem_cut = new TH2F("h_dphi_dthe_brem_cut","#Delta#phi vs #Delta#theta (brem cut (ab))",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);
  // dphi vs. dtheta for bumps that pass acceptance (filled from sb lists)
  TH2F* h_dphi_dthe_brem_cut2 = new TH2F("h_dphi_dthe_brem_cut2","#Delta#phi vs #Delta#theta (brem cut (sb))",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);

  TH1F* h_mcb_match_rad_of_ab = new TH1F("h_mcb_match_rad_of_ab", "R_{True} of MC brem photons that match to a bump a priori", nbin, 0, 100);
  TH1F* h_dphi_max_rec = new TH1F("h_dphi_max_rec", "#Delta#phi maximum cut from reco pT", nbin, -20, 200);
  TH1F* h_dphi_max_mc = new TH1F("h_dphi_max_mc", "#Delta#phi maximum cut from mc pT", nbin, -20, 200);

  // dphi_dtheta for all merged phi-bumps
  TH2F* h_dphi_dthe_mrg = new TH2F("h_dphi_dthe_mrg","#Delta#phi vs #Delta#theta (all phi bumps)",nbin,dth_min,dth_max,nbin,dphi_min*8,dphi_max*12);
  // dphi_dtheta for accepted merged phi-bumps
  TH2F* h_dphi_dthe_mrg_acc = new TH2F("h_dphi_dthe_mrg_acc","#Delta#phi vs #Delta#theta (accepted phi bumps)",nbin,dth_min,dth_max,nbin,dphi_min*8,dphi_max*12);

  TH1F* h_rec_1brem[6];
  TH1F* h_rec_1brem_found[6];
  int col_rec_1brem[6] = {1,2,4,6,7,9};
  for (int ii = 0; ii < 6; ++ii) {
    h_rec_1brem_found[ii] = th1f(Form("h_rec_1brem_found_%d",ii),"Resolution of reconstructed momentum;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
    h_rec_1brem[ii] = th1f(Form("h_rec_1brem_%d",ii),"Resolution of reconstructed momentum;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
  }

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  int nskip_nch = 0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;
    if (ientry%10000==0) cout << "Loop i: " << ientry << "/" << nentries << endl;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if (!check_invariants(ientry, true)) {
      cout << "Error on input event number " << ientry << endl;
      continue;
    }

    if (nch>1) { nskip_nch++; continue; }

    for (int ich = 0; ich < nch; ++ich) {

      if ( !is_prim[ich] ) {
      	cout << "Non Primary track  in event " << jentry << endl;
      	continue;
      }

      if ( charge[ich] <0 ) {
	//cout << "Not considering tracks with ch < 0. for positron MC" << endl;
	continue;
      }

      h_rec->Fill(res_func(mom_rec[ich],mom_mc[ich]));
      h_cor->Fill(res_func(mom_cor[ich],mom_mc[ich]));
      h_wcor->Fill(res_func(mom_wcor[ich],mom_mc[ich]));
      h_mrg->Fill(res_func(mom_mrg[ich],mom_mc[ich]));
      h_sep->Fill(res_func(mom_sep[ich],mom_mc[ich]));
      h_stored->Fill(res_func(mom_stored[ich],mom_mc[ich]));
      h_sep_w->Fill(res_func(mom_sep_w[ich],mom_mc[ich]));
      h_sep_w_bf->Fill(res_func(mom_sep_w_bf[ich],mom_mc[ich]));
      h_mrg_w->Fill(res_func(mom_mrg_w[ich],mom_mc[ich]));
      h_mrg_w_bf->Fill(res_func(mom_mrg_w_bf[ich],mom_mc[ich]));
      h_mrg_pc->Fill(res_func(mom_mrg_pc[ich],mom_mc[ich]));
      h_mrg_w_pc->Fill(res_func(mom_mrg_w_pc[ich],mom_mc[ich]));
      h_mrg_w_bf_pc->Fill(res_func(mom_mrg_w_bf_pc[ich],mom_mc[ich]));
      h_wbfcor->Fill(res_func(mom_wbfcor[ich],mom_mc[ich]));
      h_wbfcor_mw_bf->Fill(res_func(mom_wbfcor_mw_bf[ich],mom_mc[ich]));
      h_wbfcor_mw_bf_pc->Fill(res_func(mom_wbfcor_mw_bf_pc[ich],mom_mc[ich]));

      //if ( the_mc[ich] > 15 ) continue;

      float res_cor = (mom_mc[ich]-mom_cor[ich])/mom_mc[ich];
      //float res_wtd = (mom_mc[ich]-mom_wcor[ich])/mom_mc[ich];
      float res_rec = (mom_mc[ich]-mom_rec[ich])/mom_mc[ich];
      //float res_sep = (mom_mc[ich]-mom_sep[ich])/mom_mc[ich];
      //float res_mrg = (mom_mc[ich]-mom_mrg[ich])/mom_mc[ich];
      //float res_sto = (mom_mc[ich]-mom_stored[ich])/mom_mc[ich];
      //h_cor->Fill(res_cor);
      //h_wtd->Fill(res_wtd);
      //h_sep->Fill(res_sep);
      //h_mrg->Fill(res_mrg);
      //h_rec->Fill(res_rec);
      //h_sto->Fill(res_sto);

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
	if (sb_score[isb_s[ich]] != mcb_score[imcb_s[ich]] ) {
	  cout << "Single phton events non matching mcb and sb scores. This shoudlnt happen" << endl;
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

      double ptrec = mom_rec[ich]*TMath::Sin(the[ich]*TMath::DegToRad());
      double ptmc = mom_mc[ich]*TMath::Sin(the_mc[ich]*TMath::DegToRad());
      double _dph_max_rec = 2*TMath::ASin(0.12/ptrec)*TMath::RadToDeg();
      double _dph_max_mc = 2*TMath::ASin(0.12/ptmc)*TMath::RadToDeg();
      if (_dph_max_rec>175) {
	cout << "dph_max_rec " << _dph_max_rec  << " > 150: dph_max_mc= " << _dph_max_mc << " ptmc= " << ptmc << " ptrec= " << ptrec << endl;
      }
      h_dphi_max_rec->Fill(_dph_max_rec);
      h_dphi_max_mc->Fill(_dph_max_mc);
      if (_nmcb[ich]==1) { // only one brem photon
	for (int iab = 0; iab < nab; ++iab) {
	  if (ab_ich[iab]>=0) continue; // skip bumps associated with track
	  if (ab_ene[iab]<0.01*mom_rec[ich]) continue;
	  h_dphi_dthe_all->Fill(ab_the[iab]-the[ich],ab_phi[iab]-phi[ich]);
	  if ( ab_score[iab] > 0 ) {
	    h_mcb_match_rad_of_ab->Fill(mcb_rad[ab_match[iab]]);
	    if (mcb_rad[ab_match[iab]]<13)
	      h_dphi_dthe_brem->Fill(ab_the[iab]-the[ich],ab_phi[iab]-phi[ich]);
	  }
	  if ( ab_isb[iab] == ich ) {
	    h_dphi_dthe_brem_cut->Fill(ab_the[iab]-the[ich],ab_phi[iab]-phi[ich]);
	  }
	}
      }

      for (int isb = isb_s[ich]; isb < isb_e[ich]; ++isb) {
	h_dphi_dthe_brem_cut2->Fill(sb_the[isb]-the[ich],sb_phi[isb]-phi[ich]);
      }

      for (int ipb = ipb_s[ich]; ipb < ipb_e[ich]; ++ipb) {
	h_dphi_dthe_mrg->Fill(pb_the[ipb]-the[ich], pb_phi[ipb]-phi[ich]);
	if (pb_acc[ipb])
	  h_dphi_dthe_mrg_acc->Fill(pb_the[ipb]-the[ich], pb_phi[ipb]-phi[ich]);
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

  cout << "Skipped " << nskip_nch << "/" << nentries << endl;

  h_nmcb_gt1mev->Scale(100./h_nmcb_gt1mev->GetEntries());

  cout << "Writing histograms to output file: " << fout_name << endl;
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
  h_cor->Write();
  h_wcor->Write();
  h_mrg->Write();
  h_sep->Write();
  h_stored->Write();
  h_sep_w->Write();
  h_sep_w_bf->Write();
  h_mrg_w->Write();
  h_mrg_w_bf->Write();
  h_mrg_pc->Write();
  h_mrg_w_pc->Write();
  h_mrg_w_bf_pc->Write();
  h_wbfcor->Write();
  h_wbfcor_mw_bf->Write();
  h_wbfcor_mw_bf_pc->Write();

  //h_rec->Write();
  //h_wtd->Write();
  //h_cor->Write();
  //h_sep->Write();
  //h_mrg->Write();
  //h_sto->Write();

  h_out->Write();
  h_nmcb_gt1mev->Write();
  h_dphi_dthe_all->Write();
  h_dphi_dthe_brem->Write();
  h_dphi_dthe_brem_cut->Write();
  h_dphi_dthe_brem_cut2->Write();

  h_dphi_max_rec->Write();
  h_dphi_max_mc->Write();
  h_mcb_match_rad_of_ab->Write();

  h_dphi_dthe_mrg->Write();
  h_dphi_dthe_mrg_acc->Write();

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

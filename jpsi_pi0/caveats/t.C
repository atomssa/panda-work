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

static const double epbs_binSize[8]={2.25,1,1.5,1.8,2.0,2.5,3.0,4.0};

using namespace std;

void t::Loop()
{

  //gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(62);

  int nbin = 2000;
  //double xmin = -0.1, xmax = 0.3;
  double xmin = -0.5, xmax = 1.0;
  TH2F* h_rad_calc_vs_true = new TH2F("h_rad_calc_vs_true","Calculated radius vs. true radius; R_{True}; R_{Calc}", 200, 0, 40, 200, 0, 40);
  set_style(h_rad_calc_vs_true);
  TH2F* h_zed_calc_vs_true = new TH2F("h_zed_calc_vs_true","Calculated Z vs. true Z; Z_{True}; Z_{Calc}", 500, 0, 180, 500, 0, 200);
  set_style(h_zed_calc_vs_true);

  TH1F* h_rad_true = new TH1F("h_rad_true","Brem #gamma emission point R (Barrel); R_{True}[cm]", 500, 0, 45);
  TH1F* h_zed_true = new TH1F("h_zed_true","Brem #gamma emission point Z (FWD endcap); Z_{True}[cm]", 500, 0, 180);

  TH1F* h_eloss_all = new TH1F("h_eloss_all","h_eloss_all",200, 0, 1.0);
  TH2F* h_eloss_vs_the_all = new TH2F("h_eloss_vs_the_all","h_eloss_vs_the_all",200, 0, 1.0, 200, 30, 45);
  TH2F* h_eloss_vs_phi_all = new TH2F("h_eloss_vs_phi_all","h_eloss_vs_phi_all",200, 0, 1.0, 200, -180, 180);

  TH1F* h_eloss_1brem = new TH1F("h_eloss_1brem","h_eloss_1brem",200, 0, 1.0);

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

  TH1F* h_res_phi = th1f("h_res_phi",";#phi_{MC} - #phi_{KF} (deg)", nbin,-0.5,1.5,1);
  TH1F* h_res_the = th1f("h_res_the",";#theta_{MC} - #theta_{KF} (deg)", nbin,-0.5,1.5,1);
  TH1F* h_res_no_brem = th1f("h_res_no_brem","No #gamma_{Brem}>1~MeV;(p_{MC} - p_{KF})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_res_gt1_brem = th1f("h_res_gt1_brem","1+ #gamma_{Brem}>1~MeV;(p_{MC} - p_{KF})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_rec = th1f("h_rec","(p_{MC} - p_{KF})/p_{MC}", nbin,xmin,xmax,1);

  TH1F* h_sep = th1f("h_sep","(p_{MC} - p_{SEP})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_stored = th1f("h_stored","(p_{MC} - p_{STO})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_sep_w = th1f("h_sep_w","(p_{MC} - p_{SEP_W})/p_{MC}", nbin,xmin,xmax,1);
  TH1F* h_sep_w_bf = th1f("h_sep_w_bf","(p_{MC} - p_{SEP_WBF})/p_{MC}", nbin,xmin,xmax,1);

  TH1F* h_cor[nbs];
  TH1F* h_wcor[nbs];
  TH1F* h_mrg[nbs];
  TH1F* h_mrg_w[nbs];
  TH1F* h_mrg_w_bf[nbs];
  TH1F* h_mrg_pc[nbs];
  TH1F* h_mrg_w_pc[nbs];
  TH1F* h_mrg_w_bf_pc[nbs];
  TH1F* h_wbfcor[nbs];
  TH1F* h_wbfcor_mw_bf[nbs];
  TH1F* h_wbfcor_mw_bf_pc[nbs];

  for (int ibs=0; ibs < nbs; ++ibs) {
    h_cor[ibs] = th1f(Form("b%d_h_cor",ibs),"(p_{MC} - p_{COR})/p_{MC}", nbin,xmin,xmax,1);
    h_wcor[ibs] = th1f(Form("b%d_h_wcor",ibs),"(p_{MC} - p_{WTD})/p_{MC}", nbin,xmin,xmax,1);
    h_mrg[ibs] = th1f(Form("b%d_h_mrg",ibs),"(p_{MC} - p_{MRG})/p_{MC}", nbin,xmin,xmax,1);
    h_mrg_w[ibs] = th1f(Form("b%d_h_mrg_w",ibs),"(p_{MC} - p_{MRG_W})/p_{MC}", nbin,xmin,xmax,1);
    h_mrg_w_bf[ibs] = th1f(Form("b%d_h_mrg_w_bf",ibs),"(p_{MC} - p_{MRG_WBF})/p_{MC}", nbin,xmin,xmax,1);
    h_mrg_pc[ibs] = th1f(Form("b%d_h_mrg_pc",ibs),"(p_{MC} - p_{MRG_PC})/p_{MC}", nbin,xmin,xmax,1);
    h_mrg_w_pc[ibs] = th1f(Form("b%d_h_mrg_w_pc",ibs),"(p_{MC} - p_{MRG_WPC})/p_{MC}", nbin,xmin,xmax,1);
    h_mrg_w_bf_pc[ibs] = th1f(Form("b%d_h_mrg_w_bf_pc",ibs),"(p_{MC} - p_{MRG_WBFPC})/p_{MC}", nbin,xmin,xmax,1);
    h_wbfcor[ibs] = th1f(Form("b%d_h_wbfcor",ibs),"(p_{MC} - p_{SEP_WBF_MRG})/p_{MC}", nbin,xmin,xmax,1);
    h_wbfcor_mw_bf[ibs] = th1f(Form("b%d_h_wbfcor_mw_bf",ibs),"(p_{MC} - p_{SEP_WBF_MRG_WBF})/p_{MC}", nbin,xmin,xmax,1);
    h_wbfcor_mw_bf_pc[ibs] = th1f(Form("b%d_h_wbfcor_mw_bf_pc",ibs),"(p_{MC} - p_{SEP_WBF_MRG_WBFPC})/p_{MC}", nbin,xmin,xmax,1);
  }

  TH1F* h_out = th1f("h_out","Resolution W.R.T p_{MC}-E_{#gamma} #approx p_{OUT} ... (Single Brem < 42cm);((p_{MC}-E_{#gamma})-p_{KF})/(p_{MC}-E_{#gamma})",nbin,xmin,xmax,4);
  TH1F* h_nmcb_gt1mev = th1f("h_nmcb_gt1mev","Percentage of tracks with N Brem. #gamma (E>1MeV) R<42cm;N;%-age",10,0,10,1);

  const double dth_min = -20;
  const double dth_max = 20;
  const double dphi_min = -20;
  const double dphi_max = 100;

  // dphi_dtheta for all bumps
  TH2F* h_dphi_dthe_all = new TH2F("h_dphi_dthe_all","#Delta#phi vs #Delta#theta (all bumps)",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);
  // dphi vs. dtheta for MC brem photons
  TH2F* h_dphi_dthe_brem = new TH2F("h_dphi_dthe_brem","#Delta#phi vs #Delta#theta (brem bumps)",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);
  // dphi vs. dtheta for bumps that pass acceptance (filled from ab lists)
  TH2F* h_dphi_dthe_brem_cut = new TH2F("h_dphi_dthe_brem_cut","#Delta#phi vs #Delta#theta (brem cut (ab))",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);
  // dphi vs. dtheta for bumps that pass acceptance (filled from sb lists)
  TH2F* h_dphi_dthe_brem_cut2 = new TH2F("h_dphi_dthe_brem_cut2","#Delta#phi vs #Delta#theta (brem cut (sb))",nbin,dth_min,dth_max,nbin,dphi_min,dphi_max);

  TH1F* h_mcb_match_rad_of_ab = new TH1F("h_mcb_match_rad_of_ab", "R_{True} of MC brem photons that match to a bump a priori", nbin, 0, 100);
  TH1F* h_dphi_max_rec = new TH1F("h_dphi_max_rec", "#Delta#phi maximum cut from reco pT", nbin, -10, 60);
  TH1F* h_dphi_max_mc = new TH1F("h_dphi_max_mc", "#Delta#phi maximum cut from mc pT", nbin, -10, 60);
  TH1F* h_dphi_all_sep = new TH1F("h_dphi_all_sep", "#Delta#phi of all non electron bumps in event", nbin, -10, 60);
  TH1F* h_dphi_brem_sep = new TH1F("h_dphi_brem_sep", "#Delta#phi of all non electron bumps in event", nbin, -10, 60);
  TH1F* h_dphi_cut_sep = new TH1F("h_dphi_cut_sep", "#Delta#phi of all non electron bumps in event", nbin, -10, 60);

  TH1F* h_dphi_all_mrg[nbs];
  TH1F* h_dphi_cut_mrg[nbs];
  TH2F* h_dphi_dthe_mrg[nbs];
  TH2F* h_dphi_dthe_mrg_acc[nbs];
  for (int ibs=0; ibs < nbs; ++ibs) {
    h_dphi_all_mrg[ibs]  = new TH1F(Form("b%d_h_dphi_all_mrg",ibs), Form("#Delta#phi of all phi-bumps (#Delta#phi binsize= %4.2f)",epbs_binSize[ibs]), nbin, -10, 60);
    h_dphi_cut_mrg[ibs]  = new TH1F(Form("b%d_h_dphi_cut_mrg",ibs), Form("#Delta#phi of accepted phi-bumps (#Delta#phi binsize= %4.2f)",epbs_binSize[ibs]), nbin, -10, 60);
    h_dphi_dthe_mrg[ibs]  = new TH2F(Form("b%d_h_dphi_dthe_mrg",ibs),
				     Form("#Delta#phi vs #Delta#theta all phi-bumps (#Delta#phi binsize= %4.2f)",epbs_binSize[ibs]),
				     nbin,dth_min,dth_max,nbin,dphi_min*8,dphi_max*12);
    h_dphi_dthe_mrg_acc[ibs]  = new TH2F(Form("b%d_h_dphi_dthe_mrg_acc",ibs),
					 Form("#Delta#phi vs #Delta#theta accepted phi-bumps (#Delta#phi binsize= %4.2f)",epbs_binSize[ibs]),
					 nbin,dth_min,dth_max,nbin,dphi_min*8,dphi_max*12);
  }

  TH1F* h_rec_1brem[6];
  TH1F* h_out_1brem[6];
  TH1F* h_rec_1brem_found[6];
  int col_rec_1brem[6] = {1,2,4,6,7,9};
  for (int ii = 0; ii < 6; ++ii) {
    h_rec_1brem[ii] = th1f(Form("h_rec_1brem_%d",ii),"1 #gamma_{Brem};(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
    h_out_1brem[ii] = th1f(Form("h_out_1brem_%d",ii),"1 #gamma_{Brem};(p_{OUT} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
    h_rec_1brem_found[ii] = th1f(Form("h_rec_1brem_found_%d",ii),"1 #gamma_{Brem} + found;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
  }

  TH1F* h_rec_1brem_zed[6];
  TH1F* h_out_1brem_zed[6];
  TH1F* h_rec_1brem_found_zed[6];
  //int col_rec_1brem[6] = {1,2,4,6,7,9};
  for (int ii = 0; ii < 6; ++ii) {
    h_rec_1brem_zed[ii] = th1f(Form("h_rec_1brem_zed_%d",ii),"1 #gamma_{Brem};(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
    h_out_1brem_zed[ii] = th1f(Form("h_out_1brem_zed_%d",ii),"1 #gamma_{Brem};(p_{OUT} - p_{KF})/p_{OUT}",nbin,xmin,xmax,col_rec_1brem[ii]);
    h_rec_1brem_found_zed[ii] = th1f(Form("h_rec_1brem_found_zed_%d",ii),"1 #gamma_{Brem} + found;(p_{MC} - p_{KF})/p_{MC}",nbin,xmin,xmax,col_rec_1brem[ii]);
  }

  if (fChain == 0) {
    cout << "fChain is null" << endl;
    return;
  }
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

      //if (the_mc[ich] > 15) continue;

      if ( ((sim_type==t::esim)&&charge[ich]>0) ||
	   ((sim_type==t::mum)&&charge[ich]<0) ||
	   ((sim_type==t::ntype)&&charge[ich]>0) ) continue;

      h_res_phi->Fill(phi_mc[ich]-phi[ich]);
      h_res_the->Fill(the_mc[ich]-the[ich]);
      h_rec->Fill(res_func(mom_rec[ich],mom_mc[ich]));
      h_sep->Fill(res_func(mom_sep[ich],mom_mc[ich]));
      h_stored->Fill(res_func(mom_stored[ich],mom_mc[ich]));
      h_sep_w->Fill(res_func(mom_sep_w[ich],mom_mc[ich]));
      h_sep_w_bf->Fill(res_func(mom_sep_w_bf[ich],mom_mc[ich]));

      // ones that depend on binning of phi bump analysis
      for (int ibs=0; ibs < nbs; ++ibs) {
	h_cor[ibs]->Fill(res_func(mom_cor[ibs][ich],mom_mc[ich]));
	h_wcor[ibs]->Fill(res_func(mom_wcor[ibs][ich],mom_mc[ich]));
	h_mrg[ibs]->Fill(res_func(mom_mrg[ibs][ich],mom_mc[ich]));
	h_mrg_w[ibs]->Fill(res_func(mom_mrg_w[ibs][ich],mom_mc[ich]));
	h_mrg_w_bf[ibs]->Fill(res_func(mom_mrg_w_bf[ibs][ich],mom_mc[ich]));
	h_mrg_pc[ibs]->Fill(res_func(mom_mrg_pc[ibs][ich],mom_mc[ich]));
	h_mrg_w_pc[ibs]->Fill(res_func(mom_mrg_w_pc[ibs][ich],mom_mc[ich]));
	h_mrg_w_bf_pc[ibs]->Fill(res_func(mom_mrg_w_bf_pc[ibs][ich],mom_mc[ich]));
	h_wbfcor[ibs]->Fill(res_func(mom_wbfcor[ibs][ich],mom_mc[ich]));
	h_wbfcor_mw_bf[ibs]->Fill(res_func(mom_wbfcor_mw_bf[ibs][ich],mom_mc[ich]));
	h_wbfcor_mw_bf_pc[ibs]->Fill(res_func(mom_wbfcor_mw_bf_pc[ibs][ich],mom_mc[ich]));
      }

      //h_cor_vs_rec->Fill(res_recres_func(mom_rec[ich],mom_mc[ich]),res_func(mom_cor[ibs][ich],mom_mc[ich]));
      prec_vs_rec->Fill(res_func(mom_rec[ich],mom_mc[ich]),mom_rec[ich]);

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
	  h_rec_1brem_found[irad_mc]->Fill(res_func(mom_rec[ich],mom_mc[ich]));
	}
	int ized_mc = calc_ized(mcb_zed[imcb_s[ich]]);
	if (ized_mc>=0&&ized_mc<6) {
	  h_rec_1brem_found_zed[ized_mc]->Fill(res_func(mom_rec[ich],mom_mc[ich]));
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

      if (nmcb_gt1mev==0) {
	h_res_no_brem->Fill(res_func(mom_rec[ich],mom_mc[ich]));
      } else {
	h_res_gt1_brem->Fill(res_func(mom_rec[ich],mom_mc[ich]));
      }

      // Tracks that emitted only one bremstrahlung photon
      if ( _nmcb[ich]==1 ) {
	double rad_mc = mcb_rad[imcb_s[ich]];
	double zed_mc = mcb_zed[imcb_s[ich]];
	double res_rec = res_func(mom_rec[ich],mom_mc[ich]);
	double mom_out = mom_mc[ich]-mcb_ene[imcb_s[ich]];
	double res_out = res_func(mom_rec[ich],mom_out);
	// this what was meant?
	//double mom_out = mom_rec[ich]+mcb_ene[imcb_s[ich]];
	//double res_out = res_func(mom_rec[ich],mom_mc[ich]);
	if (rad_mc<42) {
	  int irad_mc = int(rad_mc/7);
	  h_rec_1brem[irad_mc]->Fill(res_rec);
	  h_out_1brem[irad_mc]->Fill(res_out);
	}
	if (zed_mc<190) {
	  int ized_mc = calc_ized(mcb_zed[imcb_s[ich]]);
	  if (ized_mc>=0&&ized_mc<6)
	    h_rec_1brem_zed[ized_mc]->Fill(res_rec);
	    h_out_1brem_zed[ized_mc]->Fill(res_out);
	}
	h_out->Fill(res_out);
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
      	  h_zed_calc_vs_true->Fill(mcb_zed[sb_match[isb]], sb_rcalc[isb]*42./TMath::Tan(TMath::DegToRad()*sb_the[isb])/30.);
	  int ized = calc_ized(mcb_zed[sb_match[isb]]);
	  if (ized>=0&&ized<6)
	    h_zed_calc[ized]->Fill(sb_rcalc[isb]/TMath::Tan(TMath::DegToRad()*sb_the[isb]));
	}
      }

      double eloss = 0;
      for (int imcb=imcb_s[ich]; imcb<imcb_e[ich]; ++imcb) {
	eloss += mcb_ene[imcb];
	if (mcb_ene[imcb]>0.001) {
	  h_rad_true->Fill(mcb_rad[imcb]);
	  h_zed_true->Fill(mcb_zed[imcb]);
	}
      }
      h_eloss_all->Fill(eloss);
      h_eloss_vs_the_all->Fill(eloss,the_mc[ich]);
      h_eloss_vs_phi_all->Fill(eloss,phi_mc[ich]);
      if (nmcb_gt1mev==1) h_eloss_1brem->Fill(eloss);

      double ptrec = mom_rec[ich]*TMath::Sin(the[ich]*TMath::DegToRad());
      double ptmc = mom_mc[ich]*TMath::Sin(the_mc[ich]*TMath::DegToRad());
      double _dph_max_rec = 2*TMath::ASin(0.12/ptrec)*TMath::RadToDeg();
      double _dph_max_mc = 2*TMath::ASin(0.12/ptmc)*TMath::RadToDeg();
      //if (_dph_max_rec>175) {
      //	cout << "dph_max_rec " << _dph_max_rec  << " > 150: dph_max_mc= " << _dph_max_mc << " ptmc= " << ptmc << " ptrec= " << ptrec << endl;
      //}
      h_dphi_max_rec->Fill(_dph_max_rec);
      h_dphi_max_mc->Fill(_dph_max_mc);
      if (_nmcb[ich]==1) { // only one brem photon
	for (int iab = 0; iab < nab; ++iab) {
	  if (ab_ich[iab]>=0) continue; // skip bumps associated with track
	  if (ab_ene[iab]<0.05*mom_rec[ich]) continue; // skip photons with energy < 1% or the electron

	  h_dphi_dthe_all->Fill(ab_the[iab]-the[ich],dphi(ab_phi[iab],ich));
	  h_dphi_all_sep->Fill(dphi(ab_phi[iab],ich));
	  if ( ab_score[iab] > 10 ) {
	    h_mcb_match_rad_of_ab->Fill(mcb_rad[ab_match[iab]]);
	    h_dphi_dthe_brem->Fill(ab_the[iab]-the[ich],dphi(ab_phi[iab],ich));
	    h_dphi_brem_sep->Fill(dphi(ab_phi[iab],ich));
	  }
	  if ( ab_isb[iab] == ich ) {
	    h_dphi_dthe_brem_cut->Fill(ab_the[iab]-the[ich],dphi(ab_phi[iab],ich));
	    h_dphi_cut_sep->Fill(dphi(ab_phi[iab],ich));
	  }
	}
      }

      for (int isb = isb_s[ich]; isb < isb_e[ich]; ++isb) {
	h_dphi_dthe_brem_cut2->Fill(sb_the[isb]-the[ich],sb_phi[isb]-phi[ich]);
      }

      for (int ibs=0; ibs < nbs; ++ibs) {
	for (int ipb = ipb_s[ibs][ich]; ipb < ipb_e[ibs][ich]; ++ipb) {
	  h_dphi_dthe_mrg[ibs]->Fill(pb_the[ibs][ipb]-the[ich], dphi(pb_phi[ibs][ipb],ich));
	  h_dphi_all_mrg[ibs]->Fill(dphi(pb_phi[ibs][ipb],ich));
	  if (pb_acc[ibs][ipb]) {
	    h_dphi_dthe_mrg_acc[ibs]->Fill(pb_the[ibs][ipb]-the[ich], dphi(pb_phi[ibs][ipb],ich));
	    h_dphi_cut_mrg[ibs]->Fill(dphi(pb_phi[ibs][ipb],ich));
	  }
	}
      }

    }
  }

  cout << "Skipped " << nskip_nch << "/" << nentries << endl;

  h_nmcb_gt1mev->Scale(100./h_nmcb_gt1mev->GetEntries());

  cout << "Writing histograms to output file: " << fout_name << endl;
  TFile *fout = TFile::Open(fout_name,"RECREATE");
  fout->cd();
  h_rad_calc_vs_true->Write();
  h_zed_calc_vs_true->Write();
  h_rad_true->Write();
  h_zed_true->Write();
  h_eloss_1brem->Write();
  h_eloss_all->Write();
  h_eloss_vs_the_all->Write();
  h_eloss_vs_phi_all->Write();
  for (int ii = 0; ii < 6; ++ii) {
    h_zed_calc[ii]->Write();
  }
  h_cor_vs_rec->Write();
  prec_vs_rec->Write();

  h_res_no_brem->Write();
  h_res_gt1_brem->Write();
  h_res_phi->Write();
  h_res_the->Write();
  h_rec->Write();
  h_sep->Write();
  h_stored->Write();
  h_sep_w->Write();
  h_sep_w_bf->Write();

  for (int ibs=0; ibs < nbs; ++ibs) {
    h_cor[ibs]->Write();
    h_wcor[ibs]->Write();
    h_mrg[ibs]->Write();
    h_mrg_w[ibs]->Write();
    h_mrg_w_bf[ibs]->Write();
    h_mrg_pc[ibs]->Write();
    h_mrg_w_pc[ibs]->Write();
    h_mrg_w_bf_pc[ibs]->Write();
    h_wbfcor[ibs]->Write();
    h_wbfcor_mw_bf[ibs]->Write();
    h_wbfcor_mw_bf_pc[ibs]->Write();
  }

  h_out->Write();
  h_nmcb_gt1mev->Write();
  h_dphi_dthe_all->Write();
  h_dphi_dthe_brem->Write();
  h_dphi_dthe_brem_cut->Write();
  h_dphi_dthe_brem_cut2->Write();

  h_mcb_match_rad_of_ab->Write();
  h_dphi_max_rec->Write();
  h_dphi_max_mc->Write();
  h_dphi_all_sep->Write();
  h_dphi_brem_sep->Write();
  h_dphi_cut_sep->Write();

  for (int ibs=0; ibs < nbs; ++ibs) {
    h_dphi_all_mrg[ibs]->Write();
    h_dphi_cut_mrg[ibs]->Write();
    h_dphi_dthe_mrg[ibs]->Write();
    h_dphi_dthe_mrg_acc[ibs]->Write();
  }

  for (int ii = 0; ii < 6; ++ii) {
    h_rec_1brem[ii]->Write();
    h_out_1brem[ii]->Write();
    h_rec_1brem_found[ii]->Write();

    h_rec_1brem_zed[ii]->Write();
    h_out_1brem_zed[ii]->Write();
    h_rec_1brem_found_zed[ii]->Write();
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

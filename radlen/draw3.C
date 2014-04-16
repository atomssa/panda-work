#include <vector>

draw3() {
  
  gStyle->SetOptStat(0);
  //gStyle->SetTitleOffset("x",0.01);
  //gStyle->SetTitleOffset("y",0.01);
  gStyle->SetLabelSize(0.2,"xyz");  
  //gStyle->SetLabelSize(0.08,"y");

  int ndet = 20;
  int start_det = 0;
  int stop_det = 6;
  bool logy = false;
  
  std::vector<TH1D*> vh_cumul_upto_emc_end_with_mvd;
  std::vector<TProfile*> vect_tprof_upto_emc_end_with_mvd;
  TFile* fin_upto_emc_end_with_mvd = TFile::Open("radlen_20k_upto_emc_end_with_mvd.root");
  for (int idet=0; idet<20; ++idet) {
    TProfile *tmp = (TProfile*)fin_upto_emc_end_with_mvd->Get(Form("RadLenProf_Det%d",idet));
    cout << "tmp= " << tmp << " tmp.entries()== " << tmp->GetEntries() << endl;
    if (tmp->GetEntries()>0) {
      vect_tprof_upto_emc_end_with_mvd.push_back( (TProfile*) fin_upto_emc_end_with_mvd->Get(Form("RadLenProf_Det%d",idet)) );
    }
  }
  std::vector<TH1D*> vh_cumul_upto_emc_end_with_mvd;
  //  for (int idet=start_det; idet<vect_tprof.size(); ++idet) {
  for (int idet=start_det; idet<stop_det; ++idet) {    
    cout << "idet = " << idet;
    vh_cumul_upto_emc_end_with_mvd.push_back(vect_tprof_upto_emc_end_with_mvd[idet]->ProjectionX(Form("RadLenProf_CumulUpToDet%d",idet)));
    cout << " vh_cumul[" << idet-start_det-1 << "]= " << vh_cumul_upto_emc_end_with_mvd[idet-start_det-1];
    cout << " vh_cumul[" << idet-start_det << "]= " << vh_cumul_upto_emc_end_with_mvd[idet-start_det] << endl;    
    if (idet>start_det)
      vh_cumul_upto_emc_end_with_mvd[idet-start_det]->Add(vh_cumul_upto_emc_end_with_mvd[idet-start_det-1]);
  }

  TPaveText *pt1= new TPaveText(0.7,0.6,0.89,0.97,"NDC");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(1);
  pt1->AddText("Registered Modules:" );
  //pt1->AddText("Cave, Magnet, Dipole, Pipe, STT, MVD, GEM");
  pt1->AddText("Cave");
  pt1->AddText("Magnet");
  pt1->AddText("Dipole");
  pt1->AddText("Pipe");
  pt1->AddText("STT");
  pt1->AddText("GEM");
  pt1->AddText("#color[3]{EMC}");
  pt1->AddText("#color[42]{MVD}");

  

  std::vector<TH1D*> vh_cumul_upto_emc_no_mvd;
  std::vector<TProfile*> vect_tprof_upto_emc_no_mvd;
  TFile* fin_upto_emc_no_mvd = TFile::Open("radlen_20k_upto_emc_no_mvd.root");
  for (int idet=0; idet<20; ++idet) {
    TProfile *tmp = (TProfile*)fin_upto_emc_no_mvd->Get(Form("RadLenProf_Det%d",idet));
    cout << "tmp= " << tmp << " tmp.entries()== " << tmp->GetEntries() << endl;
    if (tmp->GetEntries()>0) {
      vect_tprof_upto_emc_no_mvd.push_back( (TProfile*) fin_upto_emc_no_mvd->Get(Form("RadLenProf_Det%d",idet)) );
    }
  }
  std::vector<TH1D*> vh_cumul_upto_emc_no_mvd;
  //  for (int idet=start_det; idet<vect_tprof.size(); ++idet) {
  for (int idet=start_det; idet<stop_det; ++idet) {    
    cout << "idet = " << idet;
    vh_cumul_upto_emc_no_mvd.push_back(vect_tprof_upto_emc_no_mvd[idet]->ProjectionX(Form("RadLenProf_CumulUpToDet%d",idet)));
    cout << " vh_cumul[" << idet-start_det-1 << "]= " << vh_cumul_upto_emc_no_mvd[idet-start_det-1];
    cout << " vh_cumul[" << idet-start_det << "]= " << vh_cumul_upto_emc_no_mvd[idet-start_det] << endl;    
    if (idet>start_det)
      vh_cumul_upto_emc_no_mvd[idet-start_det]->Add(vh_cumul_upto_emc_no_mvd[idet-start_det-1]);
  }

  TPaveText *pt2= new TPaveText(0.7,0.55,0.89,0.97,"NDC");
  pt2->SetBorderSize(0);
  pt2->SetFillStyle(1);
  pt2->AddText("Registered Modules:" );
  pt2->AddText("Cave");
  pt2->AddText("Magnet");
  pt2->AddText("Dipole");
  pt2->AddText("Pipe");
  pt2->AddText("STT");
  pt2->AddText("GEM");  
  pt2->AddText("#color[3]{EMC}");
  

  

  
  TCanvas *tc = new TCanvas("tc","tc",1400,1000);
  tc->Divide(2,2);

  TH1D *h_gem_upto_emc_end_with_mvd = vect_tprof_upto_emc_end_with_mvd[5]->ProjectionX("RadLenProf_GEM_upto_emc_end_with_mvd");
  h_gem_upto_emc_end_with_mvd->SetLineColor(4);
  h_gem_upto_emc_end_with_mvd->SetLineWidth(3);  
  h_gem_upto_emc_end_with_mvd->SetTitle("Rad. Len. Profile of GEM;#theta[deg];X/X_{0}");    
  tc->cd(1);
  //if (logy)
  gPad->SetLogy();
  h_gem_upto_emc_end_with_mvd->SetMaximum(16);
  h_gem_upto_emc_end_with_mvd->DrawCopy("hist");
  pt1->Draw();
  TH1D *h_emc_upto_emc_end_with_mvd = vect_tprof_upto_emc_end_with_mvd[6]->ProjectionX("RadLenProf_EMC_upto_emc_end_with_mvd");
  h_emc_upto_emc_end_with_mvd->SetLineColor(2);
  h_emc_upto_emc_end_with_mvd->SetLineWidth(3);  
  h_emc_upto_emc_end_with_mvd->SetTitle("Rad. Len. Profile of EMC;#theta[deg];X/X_{0}");    
  tc->cd(2);
  //if (logy)
  gPad->SetLogy();
  //h_emc_upto_emc_end_with_mvd->SetMaximum(16);
  h_emc_upto_emc_end_with_mvd->DrawCopy("hist");
  pt1->Draw();



  TH1D *h_gem_upto_emc_no_mvd = vect_tprof_upto_emc_no_mvd[5]->ProjectionX("RadLenProf_GEM_upto_emc_no_mvd");
  h_gem_upto_emc_no_mvd->SetLineColor(4);
  h_gem_upto_emc_no_mvd->SetLineWidth(3);  
  h_gem_upto_emc_no_mvd->SetTitle("Rad. Len. Profile of GEM;#theta[deg];X/X_{0}");    
  tc->cd(3);
  //if (logy)
  gPad->SetLogy();
  h_gem_upto_emc_no_mvd->SetMaximum(16);
  h_gem_upto_emc_no_mvd->DrawCopy("hist");
  pt2->Draw();

  TH1D *h_emc_upto_emc_no_mvd = vect_tprof_upto_emc_no_mvd[6]->ProjectionX("RadLenProf_EMC_upto_emc_no_mvd");
  h_emc_upto_emc_no_mvd->SetLineColor(2);
  h_emc_upto_emc_no_mvd->SetLineWidth(3);  
  h_emc_upto_emc_no_mvd->SetTitle("Rad. Len. Profile of EMC;#theta[deg];X/X_{0}");    
  tc->cd(4);
  if (logy) gPad->SetLogy();
  h_emc_upto_emc_no_mvd->DrawCopy("hist");
  pt2->Draw();

  
//TH1D *h_gem_emc_upto_emc_no_mvd = h_gem_upto_emc_no_mvd->Clone("RadLenProf_GEM_EMC_upto_emc_no_mvd");
//h_gem_emc_upto_emc_no_mvd->Add(h_emc_upto_emc_no_mvd);
//tc->cd(4);
//if (logy) gPad->SetLogy();
//h_gem_emc_upto_emc_no_mvd->SetLineColor(1);
//h_gem_emc_upto_emc_no_mvd->SetLineWidth(6);    
//h_gem_emc_upto_emc_no_mvd->SetTitle("Rad. Len. Profile of EMC and GEM;#theta[deg];X/X_{0}");
//TLegend *tl = new TLegend(0.35,0.67,0.65,0.87);
//
//h_gem_emc_upto_emc_no_mvd->DrawCopy("hist");
//
//h_emc_upto_emc_no_mvd->DrawCopy("hist same");
//tl->AddEntry(h_emc_upto_emc_no_mvd,"EMC","pl");
//h_gem_upto_emc_no_mvd->DrawCopy("hist same");
//tl->AddEntry(h_gem_upto_emc_no_mvd,"GEM","pl");  
////h_gem_emc_upto_emc_no_mvd->DrawCopy("hist same");
//tl->AddEntry(h_gem_emc_upto_emc_no_mvd,"EMC+GEM","pl");
//tl->Draw();
//pt2->Draw();
//
  
}

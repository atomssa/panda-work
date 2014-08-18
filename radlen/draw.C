#include <vector>

draw(string fname="radlen_20k_upto_emc.root") {

  gStyle->SetOptStat(0);
  //gStyle->SetTitleOffset("x",0.01);
  //gStyle->SetTitleOffset("y",0.01);
  gStyle->SetLabelSize(0.2,"xyz");  
  //gStyle->SetLabelSize(0.08,"y");
  
  int ndet = 20;
  TFile* fin = TFile::Open(fname.c_str());
  std::vector<TH1D*> vh_cumul;
  std::vector<TProfile*> vect_tprof;

  for (int idet=0; idet<20; ++idet) {
    TProfile *tmp = (TProfile*)fin->Get(Form("RadLenProf_Det%d",idet));
    cout << "tmp= " << tmp << " tmp.entries()==" << tmp->GetEntries() << endl;
    if (tmp->GetEntries()>0) {
      vect_tprof.push_back( (TProfile*) fin->Get(Form("RadLenProf_Det%d",idet)) );
    }
  }

  std::vector<TH1D*> vh_cumul;
  int start_det = 0;
  int stop_det = 6;
  //  for (int idet=start_det; idet<vect_tprof.size(); ++idet) {
  for (int idet=start_det; idet<stop_det; ++idet) {    
    cout << "idet = " << idet;
    vh_cumul.push_back(vect_tprof[idet]->ProjectionX(Form("RadLenProf_CumulUpToDet%d",idet)));
    cout << " vh_cumul[" << idet-start_det-1 << "]= " << vh_cumul[idet-start_det-1];
    cout << " vh_cumul[" << idet-start_det << "]= " << vh_cumul[idet-start_det] << endl;    
    if (idet>start_det)
      vh_cumul[idet-start_det]->Add(vh_cumul[idet-start_det-1]);
  }
  
  //for (int idet=19; idet>=0; --idet) {
  //  TProfile *tmp = (TProfile*)fin->Get(Form("RadLenProf_Det%d",idet));
  //  cout << "tmp= " << tmp << " tmp.entries()==" << tmp->GetEntries() << endl;
  //  if (tmp->GetEntries()>0) {
  //    vect_h_cumul.push_back( (TH1D*) fin->Get(Form("RadLenProf_CumulUpToDet%d",idet)) );
  //  }
  //}

  cout <<  "Drawing vh_cumul" << endl;
  int col[] = {40,41,42,43,44,45,46,47,48,49,1,2,3,4,5,6,7,8,9,11};
  TCanvas *tc = new TCanvas("tc","tc");
  tc->cd();
  gPad->SetLogy();
  bool first = true;
  for (int i=vh_cumul.size()-1; i>=0; --i) {
    cout << i << endl;
    vh_cumul[i]->SetLineColor(col[i]);
    vh_cumul[i]->SetFillColor(col[i]);
    vh_cumul[i]->SetFillStyle(1001);
    vh_cumul[i]->Draw((first?"hist":"hist same"));
    bool first = false;
  }

  TCanvas *tc_pipe = new TCanvas("tc_pipe","tc_pipe",1400,1000);
  tc_pipe->cd(1);
  TH1D *h_pipe = vect_tprof[3]->ProjectionX("RadLenProf_PIPE");
  h_pipe->SetLineColor(48);
  h_pipe->SetLineWidth(3);  
  h_pipe->SetTitle("Radiation Length Profile of PIPE;#theta[deg];X/X_{0}");    
  gPad->SetLogy();
  h_pipe->DrawCopy("hist");
  tc_pipe->Print("pipe.eps");  
  
  TCanvas *tc_stt = new TCanvas("tc_stt","tc_stt",1400,1000);
  tc_stt->cd(1);
  TH1D *h_stt = vect_tprof[4]->ProjectionX("RadLenProf_STT");
  h_stt->SetLineColor(2);
  h_stt->SetLineWidth(3);  
  h_stt->SetTitle("Radiation Length Profile of STT;#theta[deg];X/X_{0}");    
  gPad->SetLogy();
  h_stt->DrawCopy("hist");
  tc_stt->Print("stt.eps");  

  TCanvas *tc_mvd = new TCanvas("tc_mvd","tc_mvd",1400,1000);
  tc_mvd->cd(2);
  TH1D *h_mvd = vect_tprof[7]->ProjectionX("RadLenProf_MVD");  
  h_mvd->SetLineColor(4);
  h_mvd->SetLineWidth(3);
  h_mvd->SetTitle("Radiation Length Profile of MVD;#theta[deg];X/X_{0}");
  gPad->SetLogy();
  h_mvd->DrawCopy("hist");
  tc_mvd->Print("mvd.eps");    

  TCanvas *tc_gem = new TCanvas("tc_gem","tc_gem",1400,1000);
  tc_gem->cd(1);
  TH1D *h_gem = vect_tprof[5]->ProjectionX("RadLenProf_GEM");
  h_gem->SetLineColor(42);
  h_gem->SetLineWidth(3);  
  h_gem->SetTitle("Radiation Length Profile of GEM;#theta[deg];X/X_{0}");    
  gPad->SetLogy();
  h_gem->DrawCopy("hist");
  tc_gem->Print("gem.eps");  

  TCanvas *tc_emc = new TCanvas("tc_emc","tc_emc",1400,1000);
  tc_emc->cd(1);
  TH1D *h_emc = vect_tprof[6]->ProjectionX("RadLenProf_EMC");
  h_emc->SetLineColor(45);
  h_emc->SetLineWidth(3);  
  h_emc->SetTitle("Radiation Length Profile of EMC;#theta[deg];X/X_{0}");    
  gPad->SetLogy();
  h_emc->DrawCopy("hist");
  tc_emc->Print("emc.eps");  


  TCanvas *tc_gem_emc = new TCanvas("tc_gem_emc","tc_gem_emc",1400,1000);
  TH1D *h_gem_emc = h_gem->Clone("RadLenProf_GEM_EMC");
  h_gem_emc->Add(h_emc);
  tc_gem_emc->cd(3);
  gPad->SetLogy();
  h_gem_emc->SetLineColor(1);
  h_gem_emc->SetLineWidth(3);    
  h_gem_emc->SetTitle("Radiation Length Profile of EMC and GEM;#theta[deg];X/X_{0}");
  TLegend *tl = new TLegend(0.35,0.67,0.65,0.87);
  h_gem_emc->DrawCopy("hist");
  h_gem->DrawCopy("hist same");
  tl->AddEntry(h_emc,"EMC","pl");
  h_gem->DrawCopy("hist same");
  tl->AddEntry(h_gem,"GEM","pl");  
  h_gem_emc->DrawCopy("hist same");
  tl->AddEntry(h_gem_emc,"EMC+GEM","pl");
  tl->Draw();
  tc_gem_emc->Print("gem_emc.eps");

  TCanvas *tc_stt_mvd = new TCanvas("tc_stt_mvd","tc_stt_mvd",1400,1000);
  TH1D *h_stt_mvd = h_stt->Clone("RadLenProf_STT_MVD");
  h_stt_mvd->Add(h_mvd);
  tc_stt_mvd->cd(3);
  gPad->SetLogy();
  h_stt_mvd->SetLineColor(1);
  h_stt_mvd->SetLineWidth(3);    
  h_stt_mvd->SetTitle("Radiation Length Profile of MVD and STT;#theta[deg];X/X_{0}");
  TLegend *tl = new TLegend(0.35,0.67,0.65,0.87);
  h_stt_mvd->DrawCopy("hist");
  h_mvd->DrawCopy("hist same");
  tl->AddEntry(h_mvd,"MVD","pl");
  h_stt->DrawCopy("hist same");
  tl->AddEntry(h_stt,"STT","pl");  
  h_stt_mvd->DrawCopy("hist same");
  tl->AddEntry(h_stt_mvd,"MVD+STT","pl");
  tl->Draw();
  tc_stt_mvd->Print("stt_mvd.eps");


  TCanvas *tc_stt_mvd_pipe = new TCanvas("tc_stt_mvd_pipe","tc_stt_mvd_pipe",1400,1000);
  TH1D *h_stt_mvd_pipe = h_stt->Clone("RadLenProf_STT_MVD_PIPE");
  h_stt_mvd_pipe->Add(h_mvd);
  h_stt_mvd_pipe->Add(h_pipe);
  tc_stt_mvd_pipe->cd(3);
  gPad->SetLogy();
  h_stt_mvd_pipe->SetMinimum(5e-3);
  h_stt_mvd_pipe->SetLineColor(1);
  h_stt_mvd_pipe->SetLineWidth(3);    
  h_stt_mvd_pipe->SetTitle("Radiation Length Profile of MVD, STT and PIPE;#theta[deg];X/X_{0}");
  TLegend *tl = new TLegend(0.35,0.67,0.65,0.87);
  h_stt_mvd_pipe->DrawCopy("hist");
  h_mvd->DrawCopy("hist same");
  tl->AddEntry(h_mvd,"MVD","pl");
  h_stt->DrawCopy("hist same");
  tl->AddEntry(h_stt,"STT","pl");
  h_pipe->DrawCopy("hist same");
  tl->AddEntry(h_pipe,"PIPE","pl");
  h_stt_mvd_pipe->DrawCopy("hist same");
  tl->AddEntry(h_stt_mvd_pipe,"MVD+STT+PIPE","pl");
  tl->Draw();
  tc_stt_mvd_pipe->Print("stt_mvd_pipe.eps");

  
  TCanvas *tc_stt_mvd_gem_pipe = new TCanvas("tc_stt_mvd_gem_pipe","tc_stt_mvd_gem_pipe",1400,1000);
  TH1D *h_stt_mvd_gem_pipe = h_stt->Clone("RadLenProf_STT_MVD_GEM_PIPE");
  h_stt_mvd_gem_pipe->Add(h_mvd);
  h_stt_mvd_gem_pipe->Add(h_gem);
  h_stt_mvd_gem_pipe->Add(h_pipe);
  tc_stt_mvd_gem_pipe->cd(3);
  gPad->SetLogy();
  h_stt_mvd_gem_pipe->SetMinimum(5e-3);
  h_stt_mvd_gem_pipe->SetLineColor(1);
  h_stt_mvd_gem_pipe->SetLineWidth(3);    
  h_stt_mvd_gem_pipe->SetTitle("Radiation Length Profile of MVD, STT, GEM and PIPE;#theta[deg];X/X_{0}");
  TLegend *tl = new TLegend(0.55,0.7,0.85,0.9);
  h_stt_mvd_gem_pipe->DrawCopy("hist");
  h_mvd->DrawCopy("hist same");
  tl->AddEntry(h_mvd,"MVD","pl");
  h_stt->DrawCopy("hist same");
  tl->AddEntry(h_stt,"STT","pl");
  h_gem->DrawCopy("hist same");
  tl->AddEntry(h_gem,"GEM","pl");
  h_pipe->DrawCopy("hist same");
  tl->AddEntry(h_pipe,"PIPE","pl");
  h_stt_mvd_gem_pipe->DrawCopy("hist same");
  tl->AddEntry(h_stt_mvd_gem_pipe,"MVD+STT+GEM+PIPE","pl");
  tl->Draw();
  tc_stt_mvd_gem_pipe->Print("stt_mvd_gem_pipe.eps");


  
}

void plot_pvsm(int itag) {

  gStyle->SetOptStat(0);
  
  vector<string> tags;
  tags.push_back("G3_HADR3_NUCRIN");
  tags.push_back("QGSP_BERT_EMV");
  tags.push_back("QGSP_BIC_EMV");
  tags.push_back("QGSP_BERT_EMV_OPTICAL");
  tags.push_back("QGSP_BIC_EMV_OPTICAL");
  tags.push_back("QGSC_BERT_EMV");
  tags.push_back("QGSP_BERT_HP");

  cout << "Running for tag " << tags[itag] << endl;
  
  static const int nmom = 6;
  double mom_array[nmom] = {0.5, 0.8, 1.0, 1.5, 3.0, 5.0};


  TH1F* h_pip_fp;
  TH1F* h_pip[nmom];
  TH1F* h_pim_fp;
  TH1F* h_pim[nmom];
  TFile* f_pip_fp;
  TFile* f_pip[nmom];
  TFile* f_pim_fp;
  TFile* f_pim[nmom];
  TLegend* tl_pm[nmom];
  TCanvas *tc_pm = new TCanvas("tc_pm","tc_pm");
  tc_pm->Divide(3,2);
  
  for (int imom=0; imom<nmom; ++imom) {

    //tl_pm[imom] = new TLegend(0.45,0.7,0.89,0.89);
    if (imom == 2 && (itag==0||itag==1||itag==2)) {
      tl_pm[imom] = new TLegend(0.12,0.12,0.75,0.42);
    } else {
      tl_pm[imom] = new TLegend(0.12,0.12,0.75,0.32);
    }
    tl_pm[imom]->SetFillStyle(0);
    tl_pm[imom]->SetFillColor(0);
    //tl_pm[imom]->SetBorderSize(0);
    tl_pm[imom]->SetHeader("      #color[1]{E_{reco}/E_{true}}");
    
    string fname_pip = Form("output/hist_%s_MOM_%3.1f_pip.root",tags[itag].c_str(),mom_array[imom]);
    if (itag==3||itag==4||itag==6) {
      fname_pip = Form("output/_hist_%s_MOM_%3.1f_pip.root",tags[itag].c_str(),mom_array[imom]);
    }
    cout << "Opening hists file " << fname_pip << endl;
    f_pip[imom] = TFile::Open(fname_pip.c_str());
    h_pip[imom] = (TH1F*) f_pip[imom]->Get("h_energy");
    h_pip[imom]->SetTitle(Form("%s; E_{reco}/E_{true}; Yield",tags[itag].c_str()));
    h_pip[imom]->SetLineWidth(2);
    h_pip[imom]->SetLineColor(2);
    h_pip[imom]->GetXaxis()->SetRangeUser(0.0,1.05);
    tl_pm[imom]->AddEntry(h_pip[imom],Form("#pi^{+}, p = %3.1f GeV/c",mom_array[imom]),"l");
      
    string fname_pim = Form("output/hist_%s_MOM_%3.1f_pim.root",tags[itag].c_str(),mom_array[imom]);
    cout << "Opening hists file " << fname_pim << endl;
    f_pim[imom] = TFile::Open(fname_pim.c_str());
    h_pim[imom] = (TH1F*) f_pim[imom]->Get("h_energy");
    h_pim[imom]->SetTitle(Form("%s; E_{reco}/E_{true}; Yield",tags[itag].c_str()));
    h_pim[imom]->SetLineWidth(2);
    h_pim[imom]->SetLineColor(4);
    h_pim[imom]->GetXaxis()->SetRangeUser(0.0,1.05);
    tl_pm[imom]->AddEntry(h_pim[imom],Form("#pi^{-}, p = %3.1f GeV/c",mom_array[imom]),"l");    


    if (imom == 2) {
      if (itag==0||itag==1||itag==2) {
	string fname_pip_fp = Form("output/hist_FULLPANDA_%s_MOM_1.0_pip.root",tags[itag].c_str());
	cout << "Opening hists file " << fname_pip_fp << endl;
	f_pip_fp = TFile::Open(fname_pip_fp.c_str());
	h_pip_fp = (TH1F*) f_pip_fp->Get("h_energy");
	h_pip_fp->SetTitle(Form("%s Full PANDA; E_{reco}/E_{true}; Yield",tags[itag].c_str()));
	h_pip_fp->SetLineWidth(2);
	h_pip_fp->SetLineColor(3);
	h_pip_fp->GetXaxis()->SetRangeUser(0.0,1.05);
	tl_pm[imom]->AddEntry(h_pip_fp,Form("#pi^{+}, p = 1.0 GeV/c (FP)"),"l");
      
	string fname_pim_fp = Form("output/hist_FULLPANDA_%s_MOM_1.0_pim.root",tags[itag].c_str());
	cout << "Opening hists file " << fname_pim_fp << endl;
	f_pim_fp = TFile::Open(fname_pim_fp.c_str());
	h_pim_fp = (TH1F*) f_pim_fp->Get("h_energy");
	h_pim_fp->SetTitle(Form("%s Full PANDA; E_{reco}/E_{true}; Yield",tags[itag].c_str()));
	h_pim_fp->SetLineWidth(2);
	h_pim_fp->SetLineColor(5);
	h_pim_fp->GetXaxis()->SetRangeUser(0.0,1.05);
	tl_pm[imom]->AddEntry(h_pim_fp,Form("#pi^{-}, p = 1.0 GeV/c (FP)"),"l");    
      }
    }
    
    
    tc_pm->cd(1+imom);
    gPad->SetLogy();
    h_pim[imom]->Draw();
    h_pip[imom]->Draw("same");
    if (imom==2) {
      if (itag==0||itag==1||itag==2) {
	h_pim_fp->Draw("same");
	h_pip_fp->Draw("same");
      }
    }
    tl_pm[imom]->Draw();

  }


  tc_pm->Print(Form("pvsm_%s.eps",tags[itag].c_str()));

  
}

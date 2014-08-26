#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "Riostream.h"
#include <iostream>

void plot_pippim() {

  gStyle->SetLabelSize(0.06,"X");
  gStyle->SetLabelSize(0.05,"Y");
  //gStyle->SetLabelOffset(0.01);

  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");  
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1.3,"Y");  

  gStyle->SetPadLeftMargin(0.125);
  
  vector<string> pathes;
  pathes.push_back("fevrier2014/data_pbarp_v_pippim/pbarp_pippim_079_123.dat");
  pathes.push_back("fevrier2014/data_pbarp_v_pippim/pbarp_pippim_130_186.dat");
  pathes.push_back("fevrier2014/data_pbarp_v_pippim/pbarp_pippim_191_223.dat");
  pathes.push_back("fevrier2014/data_pbarp_v_pippim/pbarp_pippim_233_243.dat"); 
  vector<int> nene;
  nene.push_back(6);
  nene.push_back(8);
  nene.push_back(4);
  nene.push_back(2);

  double ene[20] = {0.79,0.86,0.99,1.09,1.14,1.23,1.3,1.36,1.43,1.5,1.6,1.71,1.81,1.86,1.91,2.01,2.12,2.23,2.33,2.43};
  double sigma[20][48];
  double sigma_e[20][48];
  double costh[48];

  ifstream inf;
  int offset = 0;
  for (int ipath=0; ipath<4; ++ipath) {
    inf.open(pathes[ipath].c_str());
    //cout << "Opening file " << pathes[ipath] << endl;
    if (ipath>0) offset += nene[ipath-1];
    for (int idat=0; idat<48; ++idat) {
      inf >> costh[idat];
      if (!inf.good()) cout << "TROUBLE costh reading ipath= " << ipath << " idat= " << idat << endl;
      int iene = 0;
      for (int icol=0; icol<nene[ipath]; ++icol) {
	inf >> sigma[offset+iene][idat] >> sigma_e[offset+iene][idat];
	if (!inf.good()) cout << "sigma reading ipath " << ipath << " icol= " << icol << " sigma= "
			      << sigma[offset+iene][idat] << " sigma_e= " << sigma_e[offset+iene][idat] << endl;	
	++iene;
      }
    }
    inf.close();
  }

  //for (int i=0; i<48; ++i) {
  //  for (int j=0; j<20; ++j) {
  //    cout << Form("%7.2f",sigma[j][i]) << " ";
  //  }
  //  cout << endl;
  //}

  const int ncanv = 2;
  TCanvas *tc[ncanv];
  for (int ic=0; ic<ncanv; ++ic) {
    //tc[ic] = new TCanvas(Form("tc%d",ic),Form("tc%d",ic));    
    tc[ic] = new TCanvas(Form("tc%d",ic),Form("tc%d",ic),1400,1000);
    tc[ic]->Divide(4,3);
  }
  
  TGraphErrors *tge[20];

  for (int iene=0; iene<20; ++iene) {
    tge[iene] = new TGraphErrors(48,costh,sigma[iene],0,sigma_e[iene]);
    tge[iene]->SetTitle(Form("#bar{p}p#rightarrow#pi^{+}#pi{-} at E = %4.2f GeV;cos(#theta^{*});d#sigma/d(cos(#theta^{*})) [#mu b/sr]",ene[iene]));
    tge[iene]->SetMarkerStyle(20);
    tge[iene]->SetMarkerColor(2);
    tge[iene]->SetLineWidth(2);
    tge[iene]->SetLineColor(2);
    cout << "canv= " << iene/12 << " pad= " << 1+iene%12 << endl;
    tc[iene/12]->cd(1+iene%12);
    tge[iene]->Draw("ap");
  }  

  for (int ic=0; ic<ncanv; ++ic) {
    tc[ic]->Print( Form("pbarp_pippim_%d.gif",ic ));
  }
  
}


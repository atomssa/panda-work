#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "Riostream.h"
#include <iostream>

void plot(const char* config) {

  gStyle->SetLabelSize(0.06,"X");
  gStyle->SetLabelSize(0.05,"Y");
  //gStyle->SetLabelOffset(0.01);

  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");  
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1.3,"Y");  

  gStyle->SetPadLeftMargin(0.125);

  //int ene[nene] = {2911,2950,2975,2979,2981,2985,2990,2994,3005,3097,3526,3592,3617,4274};
  //int asym[nene] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1};
  
  ifstream inf;
  inf.open(config);
  int nfile = 0;
  vector<string> tag;
  vector<string> title;
  vector<string> path;
  vector<int> asym;
  vector<double> ene;
  while (1) {
    string stmp1,stmp2,stmp3;
    int itmp1, itmp2;
    inf >> stmp1 >> stmp2 >> stmp3 >> itmp1 >> itmp2;
    if (!inf.good()) break;
    tag.push_back(stmp1);
    title.push_back(stmp2);
    path.push_back(stmp3);
    asym.push_back(itmp1-3);
    ene.push_back(itmp2/1000.);
    //cout << "test: " << tag[tag.size()-1] << "  " << path[path.size()-1] << "  " << asym[asym.size()-1] << " " << ene[ene.size()-1] << endl;
  }
  inf.close();

  TCanvas *tc0 = new TCanvas("tc0","tc0",1400,1000);
  tc0->Divide(4,3);
  TCanvas *tc1 = new TCanvas("tc1","tc1",1400,1000);  
  tc1->Divide(4,3);
  

  TGraphErrors *tge;
  //cout << "nene = " << tag.size() << endl;
  for (int iene=0; iene<tag.size(); ++iene) {

    double costh[50];
    double sigma[50];
    double sigma_e[50];
    double sigma_eh[50];  
    double sigma_el[50];

    int idat=0;
    inf.open(path[iene].c_str());
    while(1) {
      if (!inf.good()) break;
      if (asym[iene]) {
	inf >> costh[idat] >> sigma[idat] >> sigma_el[idat] >> sigma_eh[idat];;
	sigma_el[idat] *= -1.0;
      } else {
	inf >> costh[idat] >> sigma[idat] >> sigma_e[idat];
      }
      ++idat;
    }
    inf.close();
    //cout << "Plotting E= " << ene[iene] << " from file " << path[iene] << " idat= " << idat << " asym= " << asym[iene] <<endl;

    if (asym[iene]) {
      tge = new TGraphAsymmErrors(idat-1,costh,sigma,0,0,sigma_el,sigma_eh);
    } else {
      tge = new TGraphErrors(idat-1,costh,sigma,0,sigma_e);
    }
    tge->SetTitle(Form("%s at E = %4.2f GeV;cos(#theta);d#sigma/d(cos(#theta)) [nb]",title[iene].c_str(),ene[iene]));
    tge->SetMarkerStyle(20);
    tge->SetMarkerColor(2);
    tge->SetLineWidth(2);
    tge->SetLineColor(2);
    if (iene<12) {
      tc0->cd(iene+1);
    } else {
      tc1->cd(iene-11);
    }
    tge->Draw("ap");
    
  }  

  tc0->Print( Form("%s_1.gif",tag[0].c_str()) );
  tc1->Print( Form("%s_2.gif",tag[0].c_str()) );

}

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"

void extractHisto(TString& histoname, TFile& storagefile, TCanvas* mycanvas, TString& myprefix)
{
	if (storagefile.Get(histoname.Data()) != 0) {
		TH1D* myhisto = (TH1D*)storagefile.Get(histoname.Data());
		TString outfilename = myprefix + histoname + ".pdf";
		mycanvas->cd();
		myhisto->Draw();
		mycanvas->SaveAs(outfilename.Data());
	}
}

int extract_histos(TString infilename, TString prefix="psi3770plot_")
{

//TStopwatch timer;
//timer.Start();

//gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");rootlogon();

TFile inputfile(infilename, "READ");

TCanvas* outputcanvas = new TCanvas();

extractHisto("hniceevents", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexxfit", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexyfit", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexzfit", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexxfit", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexyfit", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexzfit", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexxreco", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexyreco", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexzreco", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexxreco", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexyreco", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexzreco", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexxmc", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexymc", inputfile, outputcanvas, prefix);
extractHisto("hdpvertexzmc", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexxmc", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexymc", inputfile, outputcanvas, prefix);
extractHisto("hdmvertexzmc", inputfile, outputcanvas, prefix);
extractHisto("hfinalpsimass", inputfile, outputcanvas, prefix);

//extractHisto("", inputfile, outputcanvas, prefix);

return 0;

//timer.Stop();
//Double_t rtime = timer.RealTime();
//Double_t ctime = timer.CpuTime();
    
//printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
    
} // end macro

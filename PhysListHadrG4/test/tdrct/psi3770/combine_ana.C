#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"

void SaveAndUpdateHisto(TH1D* currenthisto, TFile& storagefile)
{
	if (storagefile.Get(currenthisto->GetName()) != 0) {
		currenthisto->Add((TH1D*)storagefile.Get(currenthisto->GetName()));
	}
	//cout << currenthisto->GetName() << ": " << currenthisto->GetEntries() << endl;
	currenthisto->Write();
}

int combine_ana(TString infilename, TString outfilename="combined_psi3770.root")
{

//TStopwatch timer;
//timer.Start();

//gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");rootlogon();

TFile inputfile(infilename, "READ");
TList* mylist = inputfile.GetListOfKeys();
TListIter myiter(mylist);
TFile outputstorage(outfilename, "UPDATE");
TKey* mykey;

//cout << mylist->GetEntries() << endl;

int nHistos = 149;
TString histoClassName = "TH1D";

if (mylist->GetEntries() != nHistos) {
	cout << "Fatal Entry Mismatch: " << mylist->GetEntries() << endl;
	return 1;
}

outputstorage->cd();

while (mykey = (TKey*)myiter()) {
	if (histoClassName == TString(mykey->GetClassName())) {
		TH1D* myhisto = (TH1D*)mykey->ReadObj();
		SaveAndUpdateHisto(myhisto, outputstorage);
	} else {
		cout << "Ignored, not a histogram: " << mykey->GetName() << ", " << mykey->GetClassName() << endl;
	}
}

return 0;

//timer.Stop();
//Double_t rtime = timer.RealTime();
//Double_t ctime = timer.CpuTime();
    
//printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
    
} // end macro

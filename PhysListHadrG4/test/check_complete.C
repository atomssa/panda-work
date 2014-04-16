#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TKey.h"
#include <iostream>
#include "TH1F.h"
#include "TRandom.h"

using std::cout;
using std::endl;

bool check_complete(TString fn="output_ana.root", TString fn2="ana_target.root", double minP = 0.03, int minev = 3, int maxfail=3)
{
	bool fTest=kFALSE;
        TString templateFile = gSystem->Getenv("VMCWORKDIR");
        templateFile += "/macro/run/";
        templateFile += fn2;
	
	TFile *f=new TFile(fn,"READ");
	if (!f->IsZombie())
	{
		TFile *f2=new TFile(templateFile,"READ");
			
		TKey *key;
		TIter next(f->GetListOfKeys());
		
		int failcount = 0;
		
		while ( (key = (TKey*)next()) )
		{
			TObject *obj = key->ReadObj();
			
			// only check TH1Fs
			if (!obj->InheritsFrom("TH1F")) continue;
			
			TString name = obj->GetName();
			TH1F* h  = (TH1F*) obj; 
			TH1F* h2 = (TH1F*) f2->Get(name);
			
			if ( h->GetEntries()<minev ) 
			{
				cout << "Histogram (almost) empty : " << name << " \"" << h2->GetTitle() << "\":  N = " <<  h->GetEntries() << endl;
				failcount++;
			}
			else  
			{
				double P = h2->KolmogorovTest(h);
				if ( P<minP )
				{
					cout << "Incompatible distribution: " << name << " \"" << h2->GetTitle() << "\":  P = " << P << endl;
					failcount++;
				}
			}
		
		}
		
		if (failcount<maxfail) fTest = kTRUE;
	}
	
	if (fTest){
		cout << " Test passed" << endl;
		cout << " All ok " << endl;  
	}else{
		cout << " Test Failed" << endl;
		cout << " Not Ok " << endl;         
	}
    
    exit(fTest);
}

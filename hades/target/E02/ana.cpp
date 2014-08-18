#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TSystem.h"

//#include "Ex02MCParticle.h"
#include "include/Ex02MCStack.h"

#include <iostream>

using namespace std;

int main() {

	//	gSystem->Load("lib/tgt_linuxx8664gcc/libexample02");

	//TFile *f = TFile::Open("example02.root");

	TFile *f = new TFile("example02.root","READ");

	TTree *t = (TTree*) f->Get("example02");

	//TClonesArray *part_array = new TClonesArray("Ex02MCStack");
	Ex02MCStack *part_array;
	t->SetBranchAddress("stack",&part_array);

	//TClonesArray *_array = new TClonesArray("Ex02MCParticle");
	//t->SetBranchAddress("??",&part_array);

	const int nent = t->GetEntriesFast();

	for (int ient = 0; ient < nent; ++ient) {
		
		t->GetEntry(ient);

		const int npart = part_array->GetNtrack();

		cout << "npart = " << npart << endl;

	}
	
	return 0;
}
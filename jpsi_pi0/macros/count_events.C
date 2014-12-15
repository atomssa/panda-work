#include "TFile.h"
#include "TTree.h"
#include "Riostream.h"

void count_events(){
  const int nplab = 3;
  const int ntype = 2;
  const char* tt[ntype]= {"pip_pim", "jpsi"};
  const double plab[nplab] = {5.513, 8., 12.};
  TFile *ff[ntype][nplab];
  TTree *tree[ntype][nplab];
  int nevt[ntype][nplab][20];
  cout << "[" << endl;
  for (int itype = 0; itype < ntype; ++itype) {
    cout << "[" << endl;
    for (int iplab = 0; iplab < nplab; ++iplab) {
      const char* listfname = Form("lists/pid_%s_%3.1f.list",tt[itype],plab[iplab]);
      ifstream listf;
      listf.open(listfname);
      string fname;
      int fid;
      int ifile = 0;
      cout << "dict( [";
      while(true) {
	listf >> fid >>fname;
	if (!listf.good()) break;
	ff[itype][iplab] = TFile::Open(fname.c_str());
	tree[itype][iplab] = (TTree*) ff[itype][iplab]->Get("cbmsim");
	nevt[itype][iplab][ifile] = tree[itype][iplab]->GetEntries();
	cout << " (" << fid << ", "<< nevt[itype][iplab][ifile] << "), ";
      }
      cout << " ]),";
      listf.close();
    }
    cout << "], " << endl;
  }
  cout << "]" << endl;


}

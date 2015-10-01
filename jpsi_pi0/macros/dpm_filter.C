#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include <iostream>
#include <string>
#include <algorithm>

enum {
  phot=22,
  ep=-11,
  em=11,
  pi0=111,
  pip=211,
  pim=-211,
  k0s=310,
  kp=321,
  km=-321,
  lambda=3122,
  alambda=-3122,
  phi=333,
  jpsi=443,
  prot=2212,
  aprot=-2212,
  neut=2112,
  aneut=-2112,
  wtf=130
};

struct evt_ref {
  evt_ref(string n, string t, int i1, int i2):count(0), name(n), title(t) { pdg.push_back(i1); pdg.push_back(i2); sort_evt(); }
  evt_ref(string n, string t, int i1, int i2, int i3): evt_ref(n,t,i1,i2) { pdg.push_back(i3); sort_evt(); }
  evt_ref(string n, string t, int i1, int i2, int i3, int i4): evt_ref(n,t,i1,i2,i3) { pdg.push_back(i4); sort_evt(); }
  evt_ref(string n, string t, int i1, int i2, int i3, int i4, int i5): evt_ref(n,t,i1,i2,i3,i4) { pdg.push_back(i5); sort_evt(); }
  void sort_evt() {std::sort(pdg.begin(),pdg.end());}
  std::vector<int> pdg;
  int count;
  string name;
  string title;
};

void print_event(const vector<int> &evt) {
    cout << "---> Event pdg codes : (";
    for (unsigned int ipdg=0; ipdg<evt.size(); ++ipdg) {
      cout << evt[ipdg];
      if (ipdg!=(evt.size()-1)) cout << ",";
    }
    cout << ")" << endl;
}

void print_log(int ievt, const vector<int> &pdg_codes){
  std::cout << "Event(Type3Pi) " << ievt << " Tca.entries= " << pdg_codes.size() << std::endl;
  print_event(pdg_codes);
  cout << endl;
}

bool compare_unsorted(const vector<int> &evt, const vector<int> &_ref) {

  if (evt.size()!=_ref.size()) return false;
  vector<int> visited;
  for (unsigned int ievt=0; ievt<evt.size(); ++ievt) {

    // search in _reference, and mark index visited if foun
    bool found = false;
    for (unsigned int i_ref=0; i_ref<_ref.size(); ++i_ref) {
      if ( evt[ievt] == _ref[i_ref]) {
	found = true;

	// check if this index of _ref was already
	bool was_visited = false;
	for (unsigned int ivis=0; ivis<visited.size(); ++ivis) {
	  if (int(i_ref) == visited[ivis]) {
	    was_visited = true;
	    break;
	  }
	}
	if (was_visited) continue;

	visited.push_back(i_ref);
	break;
      }
    }

    if (!found) return false;

  }
  return true;
}

bool compare_sorted(const vector<int> &evt, const vector<int> &_ref) {
  if (evt.size()!=_ref.size()) return false;
  for (unsigned int ievt=0; ievt<evt.size(); ++ievt) {
    if (evt[ievt]!=_ref[ievt]) return false;
  }
  return true;
}

bool compare_sorted(const vector<int> &evt, evt_ref &_ref) {
  if (compare_sorted(evt, _ref.pdg)) {
    _ref.count++;
    return true;
  } else {
    return false;
  }
}

bool interesting_event(const vector<int> &evt) {
  for (unsigned int ipdg=0; ipdg<evt.size(); ++ipdg) {
    if (evt[ipdg]!=phot&&
	evt[ipdg]!=pi0&&
	evt[ipdg]!=pip&&
	evt[ipdg]!=pim&&
	evt[ipdg]!=prot&&
	evt[ipdg]!=aprot&&
	evt[ipdg]!=neut&&
	evt[ipdg]!=aneut&&
	evt[ipdg]!=k0s&&
	evt[ipdg]!=kp&&
	evt[ipdg]!=km&&
	evt[ipdg]!=lambda&&
	evt[ipdg]!=alambda&&
	evt[ipdg]!=phi&&
	evt[ipdg]!=wtf
	) return true;
  }
  return false;
}

void print_interesting_events(int ievt, const vector<int> &evt) {
  if (interesting_event(evt)) print_log(ievt, evt);
}

void dpm_filter(string rxn="pi0pippim") {

  bool verb = false;

  std::vector<evt_ref> refs;
  refs.push_back(evt_ref("pi0pipm", "#pi^{0}#pi^{+}#pi^{-}", pi0, pip, pim)); // pi-, pi0, pi+
  //refs.push_back(evt_ref("pi0epem", "#pi^{0}e^{+}e^{-}", pi0, ep, em)); // e+, e-, pi0: no e+e- events in DPM
  //refs.push_back(evt_ref("pi0jpsi", "#pi^{0}J/#psi", pi0, jpsi)); // pi0, jpsi: no jpsi events in DPM
  refs.push_back(evt_ref("pi02pipm", "#pi^{0}#pi^{0}#pi^{+}#pi^{-}", pi0, pi0, pip, pim)); // pi-, pi0, pi0, pi+
  refs.push_back(evt_ref("pi0pipm2", "#pi^{0}#pi^{+}#pi^{-}#pi^{+}#pi^{-}", pi0, pip, pim, pip, pim)); // pi-, pi-, pi0, pi+, pi+
  //refs.push_back(evt_ref("pi0pi0jpsi", "#pi^{0}#pi^{0}J/#psi", pi0, pi0, jpsi)); // pi0, pi0, jpsi: no jpsi events in DPM

  // input
  string file_name_in = "Background-micro.root";
  TFile *file_in = TFile::Open(file_name_in.c_str());
  TTree *data_in = (TTree*) file_in->Get("data");
  TClonesArray *part_array = new TClonesArray("TParticle");
  data_in->SetBranchAddress("Particles",&part_array);
  int Nevt = data_in->GetEntries();

  string file_name_out = "filt_complete.root";
  std::cout << "Filtered output file = " << file_name_out << std::endl;
  TFile *file_out = TFile::Open(file_name_out.c_str(),"RECREATE");
  std::vector<TTree*> data_out;
  for (unsigned int iref=0; iref < refs.size(); ++iref) {
    //data_out.push_back(new TTree(Form("data%d",iref),Form("DPM Background %d",iref)));
    if (rxn==refs[iref].name)
      data_out.push_back(new TTree("data",refs[iref].title.c_str()));
    else
      data_out.push_back(new TTree(refs[iref].name.c_str(),refs[iref].title.c_str()));
    data_out[iref]->Branch("Particles", &part_array, 32000, 2);
  }

  for (int ievt=0; ievt<Nevt; ievt++) {

    if (ievt%10000==0) std::cout << "event " << ievt << "/" << Nevt << std::endl;

    data_in->GetEntry(ievt);

    // veto events with more than three particles in final state
    int npart = part_array->GetEntriesFast();

    std::vector<int> pdg_codes;
    for (int ipart=0; ipart < npart; ++ipart) {
      const TParticle* part = (TParticle*) part_array->At(ipart);
      if (verb) std::cout << " part " << ipart << " pdg= " << part->GetPDG()->PdgCode() << std::endl;
      pdg_codes.push_back(part->GetPDG()->PdgCode());
    }

    print_interesting_events(ievt, pdg_codes);

    std::sort(pdg_codes.begin(),pdg_codes.end());

    for (unsigned int iref=0; iref<refs.size(); ++iref) {
      if (compare_sorted(pdg_codes, refs[iref])) {
	if (verb) print_log(ievt, pdg_codes);
	data_out[iref]->Fill();
      }
    }

  }

  file_out->cd();
  for (unsigned int iref=0; iref<refs.size(); ++iref) {
    data_out[iref]->Write();
  }
  file_out->Write();

  // summary
  cout << "Total Event Count = " << Nevt << endl;
  for (unsigned int iref=0; iref < refs.size(); ++iref) {
    cout << refs[iref].title << " event type count= " << refs[iref].count << endl;
  }

}

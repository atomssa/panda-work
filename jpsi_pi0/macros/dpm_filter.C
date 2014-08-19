#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include <iostream>
#include <string>
#include <algorithm>

bool compare_unsorted(const vector<int> &evt, const vector<int> &ref) {

  if (evt.size()!=ref.size()) return false;
  vector<int> visited;
  for (unsigned int ievt=0; ievt<evt.size(); ++ievt) {

    // search in reference, and mark index visited if foun
    bool found = false;
    for (unsigned int iref=0; iref<ref.size(); ++iref) {
      if ( evt[ievt] == ref[iref]) {
	found = true;

	// check if this index of ref was already
	bool was_visited = false;
	for (unsigned int ivis=0; ivis<visited.size(); ++ivis) {
	  if (int(iref) == visited[ivis]) {
	    was_visited = true;
	    break;
	  }
	}
	if (was_visited) continue;

	visited.push_back(iref);
	break;
      }
    }
    
    if (!found) return false;

  }
  return true;
}

bool compare_sorted(const vector<int> &evt, const vector<int> &ref) {
  if (evt.size()!=ref.size()) return false;
  for (unsigned int ievt=0; ievt<evt.size(); ++ievt) {
    if (evt[ievt]!=ref[ievt]) return false;
  }
  return true;
}

void print_event(const vector<int> &evt) {
    cout << "---> Event pdg codes : (";
    for (unsigned int ipdg=0; ipdg<evt.size(); ++ipdg) {
      cout << evt[ipdg];
      if (ipdg!=(evt.size()-1)) cout << ",";
    }
    cout << ")" << endl;
}

void print_log(int ievt, int npart, vector<int> &pdg_codes){
  std::cout << "Event(Type3Pi) " << ievt << " Tca.entries= " << npart << std::endl;
  print_event(pdg_codes);	
  cout << endl;
}


void dpm_filter(string file_name_in, string file_name_out) {

  bool verb = false;
  
  // input
  TFile *file_in = TFile::Open(file_name_in.c_str());
  TTree *data_in = (TTree*) file_in->Get("data");
  TClonesArray *part_array = new TClonesArray("TParticle");
  data_in->SetBranchAddress("Particles",&part_array);
  int Nevt = data_in->GetEntries();

  // filtered output
  //size_t s = file_name.find(".root", s);
  //if (s == string::npos) {
  //  std::cout << "You have to supply a root file" << std::endl;
  //  return;
  //}
  //file_name.replace(s,5,"_filt.root");
  std::cout << "Filtered output file = " << file_name_out << std::endl;
  TFile *file_out = TFile::Open(file_name_out.c_str(),"RECREATE");
  TTree *data_out_pipmpi0 = new TTree("data","DPM Background");
  TTree *data_out_jpsipi0 = new TTree("data_jpsi_pi0","data_jpsi_pi0");
  data_out_jpsipi0->Branch("Particles", &part_array, 32000, 0);
  data_out_pipmpi0->Branch("Particles", &part_array, 32000, 0);

  std::vector<int> ref1;
  ref1.push_back(-211);
  ref1.push_back(111);
  ref1.push_back(211);
  int type1_count = 0;
  
  std::vector<int> ref2;
  ref2.push_back(-11);
  ref2.push_back(11);
  ref2.push_back(111);
  int type2_count = 0;

  std::vector<int> ref3;
  ref3.push_back(443);
  ref3.push_back(111);
  int type3_count = 0;
  
  for (int ievt=0; ievt<Nevt; ievt++) {

    if (ievt%10000==0)
      std::cout << "event " << ievt << "/" << Nevt << std::endl;
	
    data_in->GetEntry(ievt);

    // veto events with more than three particles in final state
    int npart = part_array->GetEntriesFast();
    if (npart>3) continue;
    
    std::vector<int> pdg_codes;
    for (int ipart=0; ipart < npart; ++ipart) {
      const TParticle* part = (TParticle*) part_array->At(ipart);
      if (verb) std::cout << " part " << ipart << " pdg= " << part->GetPDG()->PdgCode() << std::endl;
      pdg_codes.push_back(part->GetPDG()->PdgCode());
    }

    std::sort(pdg_codes.begin(),pdg_codes.end());
    
    if (compare_sorted(pdg_codes,ref1)) {
      if (verb) print_log(ievt, npart, pdg_codes);
      data_out_pipmpi0->Fill();
      type1_count++;
    }

    if (compare_sorted(pdg_codes,ref2)) {
      if (verb) print_log(ievt, npart, pdg_codes);
      data_out_jpsipi0->Fill();
      type2_count++;
    }

    if (compare_sorted(pdg_codes,ref3)) {
      if (verb) print_log(ievt, npart, pdg_codes);
      data_out_jpsipi0->Fill();
      type3_count++;
    }



    
  }

  file_out->cd();
  data_out_jpsipi0->Write();
  data_out_pipmpi0->Write();
  file_out->Write();


  // summary
  cout << "Total Event Count = " << Nevt << endl;
  cout << "type1(2Pi) count = " << type1_count << endl;
  cout << "type2(JpsiPi0) count = " << type2_count << endl;
  cout << "type2(JpsiPi0Undecayed) count = " << type3_count << endl;  
  
}

#include "filler.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TClonesArray.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;

static const double rtd= TMath::RadToDeg();
static const double m_2pi = 2.0*TMath::Pi();
static const double m_pi = TMath::Pi();

bool compare_unsorted(const vector<int>&, const vector<int>&);
bool compare_sorted(const vector<int>&, const vector<int>&);

void print_event(const vector<int>&);
void print_log(int, int, vector<int>&);

bool accept_pc(const TLorentzVector &v);
bool accept_pi0(const TLorentzVector &v);
bool accept_elec(const TLorentzVector &v);

int main(const int argc, const char **argv) {


  bool verb = false;

  if (argc<2) {
    std::cout << "Need one argument [file-name]" << std::endl;
    return -1;
  }

  TChain *data_in = new TChain("data","data_in");
  ifstream inf;
  inf.open(argv[1]);
  string root_file;
  while(true) {
    inf >> root_file;
    if (!inf.good()) break;
    data_in->AddFile(root_file.c_str());
  }

  const double m_prot= 0.938;
  const double p_antip = 5.513;
  const double E_antip = TMath::Hypot(m_prot, p_antip);
  const double beta_cm = p_antip/(E_antip + m_prot);
  cout << "betac_cm = " << beta_cm << endl;
  const TVector3 boost_vector(0,0,-beta_cm);
  boost_vector.Print();
  const double m_max = 1.1 * TMath::Sqrt(2*m_prot + 2*p_antip*m_prot);

  TClonesArray *part_array = new TClonesArray("TParticle");
  data_in->SetBranchAddress("Particles",&part_array);

  std::vector<int> ref1;
  ref1.push_back(-211);
  ref1.push_back(111);
  ref1.push_back(211);

  const char *pp[] = {"pim", "pi0", "pip", "pipm"};
  const char *pp_title[] = {"#pi^{+}", "#pi^{0}", "#pi^{-}", "(#pi^{+}#pi^{-})"};
  std::vector<filler*> fillers_lab;

  // (pi0) mom and angles
  const char *lab_frame[] = {"lab","Lab"};
  fillers_lab.push_back(new mom_filler1d(1, pp, 200, 0, 6, lab_frame, pp_title));
  fillers_lab.push_back(new the_filler1d(1, pp, 200, 0, m_pi, lab_frame, pp_title));
  fillers_lab.push_back(new phi_filler1d(1, pp, 200, -m_pi, m_pi, lab_frame, pp_title));

  // (pi+) - (pi-) system mass, mom and angles
  //fillers_lab.push_back(new pair_mass_filler1d(0, 2, pp, 200, 0, 5, Form("#pi^{+}-#pi^{+} pair invariant mass;M^{inv}_{#pi^{+}-#pi^{+}}")));
  fillers_lab.push_back(new pair_mom_filler1d(0, 2, pp, 200, 0, 6, lab_frame, pp_title));
  fillers_lab.push_back(new pair_the_filler1d(0, 2, pp, 200, 0, m_pi, lab_frame, pp_title));
  fillers_lab.push_back(new pair_phi_filler1d(0, 2, pp, 200, -m_pi, m_pi, lab_frame, pp_title));
  fillers_lab.push_back(new pair_oa_filler1d(0, 2, pp, 200, 0, m_pi, lab_frame, pp_title));
  fillers_lab.push_back(new pair_oa_filler1d(1, 3, pp, 200, 0, m_pi, lab_frame, pp_title));

  // frame independent ones
  fillers_lab.push_back(new pair_mass_filler1d(0, 2, pp, 200, 0, m_max, nullptr, pp_title));
  fillers_lab.push_back(new dalitz_filler2d(0, 2, 1, pp, 200, 0, m_max, nullptr, pp_title));

  std::vector<filler*> fillers_cm;
  const char *cm_frame[] = {"cm","CM"};
  fillers_cm.push_back(new mom_filler1d(1, pp, 200, 0, 6, cm_frame, pp_title));
  fillers_cm.push_back(new the_filler1d(1, pp, 200, 0, m_pi, cm_frame, pp_title));
  fillers_cm.push_back(new pair_mom_filler1d(0, 2, pp, 200, 0, 6, cm_frame, pp_title));
  fillers_cm.push_back(new pair_the_filler1d(0, 2, pp, 200, 0, m_pi, cm_frame, pp_title));
  fillers_cm.push_back(new pair_oa_filler1d(0, 2, pp, 200, 0, m_pi, cm_frame, pp_title));
  fillers_cm.push_back(new pair_oa_filler1d(1, 3, pp, 200, 0, m_pi, cm_frame, pp_title));


  const int Nevt = data_in->GetEntries();

  for (int ievt=0; ievt<Nevt; ievt++) {

    if (ievt%10000==0)
      std::cout << "event " << ievt << "/" << Nevt << std::endl;
    data_in->GetEntry(ievt);

    // veto events with more than three particles in final state
    int npart = part_array->GetEntriesFast();
    if (npart!=3) {
      std::cout << "npart= " << npart << std::endl;
    }

    TLorentzVector p4pi0, p4pip, p4pim, p4pipm;
    int cpi0=0,cpip=0,cpim=0;
    for (int ipart=0; ipart < npart; ++ipart) {
      const TParticle* part = (TParticle*) part_array->At(ipart);
      if (verb) {
	cout << "part= " << part << endl;
	std::cout << " part " << ipart << " pdg= " << part->GetPDG()->PdgCode() << std::endl;
	part->Print();
      }
      int pdg = part->GetPDG()->PdgCode();
      // This would insure ordering
      if ( pdg == -211 ) { ((TParticle*)part_array->At(ipart))->Momentum(p4pim); cpim++; }
      if ( pdg == 111  ) { ((TParticle*)part_array->At(ipart))->Momentum(p4pi0); cpi0++; }
      if ( pdg == 211  ) { ((TParticle*)part_array->At(ipart))->Momentum(p4pip); cpip++; }
    }
    p4pipm = p4pip + p4pim;

    // check event type
    if ( cpim!=1 || cpip!=1 || cpi0!=1 ) {
      std::cout << "Event wrong type " << cpim << " pi- " << cpip << " pi+ " << cpi0 << " pi0 found where 1 expected" << std::endl;
      continue;
    }

    // check acceptance
    //if (!accept_pc(p4pim)||!accept_pc(p4pip)||accept_pi0(p4pi0)) continue;

    vector<TLorentzVector> p4s;
    p4s.push_back(p4pim);
    p4s.push_back(p4pi0);
    p4s.push_back(p4pip);
    p4s.push_back(p4pipm);

    for (auto fill: fillers_lab) (*fill)(p4s);
    for (auto fill: fillers_cm) (*fill)(p4s,boost_vector);

  }

  TFile *fout = TFile::Open("out.root","RECREATE");
  fout->cd();
  for (auto fill: fillers_lab) fill->Write();
  for (auto fill: fillers_cm) fill->Write();
  fout->Write();
  fout->Close();

  return 0;

}

bool accept_pc(const TLorentzVector &v) {
  return true;
}

bool accept_pi0(const TLorentzVector &v) {
  return true;
}

bool accept_elec(const TLorentzVector &v) {
  return true;
}


void fill(string file_name_in) {

  bool verb = false;

  // input
  TFile *file_in = TFile::Open(file_name_in.c_str());
  TTree *data_in = (TTree*) file_in->Get("data");
  TClonesArray *part_array = new TClonesArray("TParticle");
  data_in->SetBranchAddress("Particles",&part_array);

  std::vector<int> ref1;
  ref1.push_back(-211);
  ref1.push_back(111);
  ref1.push_back(211);
  int type1_count = 0;

  double m_prot= 0.938;
  double p_antip = 5.513;
  double E_antip = TMath::Hypot(m_prot, p_antip);
  double beta_cm = p_antip/(E_antip + m_prot);
  cout << "betac_cm = " << beta_cm << endl;

  TVector3 boost_vector(0,0,-beta_cm);
  boost_vector.Print();

  string tag = "bg";
  //TH1F* dthe_cm = new TH1F("dthe_cm","dthe_cm",200,179.99,180.01);
  TH1F* dthe_cm = new TH1F(Form("dthe_cm_%s",tag.c_str()),Form("dthe_cm_%s",tag.c_str()),200,0,200);
  TH1F* dthe_lab = new TH1F(Form("dthe_lab_%s",tag.c_str()),Form("dthe_lab_%s",tag.c_str()),200,0,185);
  TH1F* dphi_cm = new TH1F(Form("dphi_cm_%s",tag.c_str()),Form("dphi_cm_%s",tag.c_str()),200,179.99,180.01);

  TH1F* the_pi0_cm = new TH1F(Form("the_pi0_cm_%s",tag.c_str()),Form("%s #theta_{#pi^{0}}^{cm};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pi0_cm->SetLineColor(2);
  the_pi0_cm->SetLineWidth(2);
  TH1F* the_pippim_cm = new TH1F(Form("the_pippim_cm_%s",tag.c_str()),Form("%s #theta_{#pi^{+}#pi^{-}}^{cm};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pippim_cm->SetLineColor(4);
  the_pippim_cm->SetLineWidth(2);
  TH1F* the_pi0_lab = new TH1F(Form("the_pi0_lab_%s",tag.c_str()),Form("%s #theta_{#pi^{0}}^{lab};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pi0_lab->SetLineColor(2);
  the_pi0_lab->SetLineWidth(2);
  TH1F* the_pippim_lab = new TH1F(Form("the_pippim_lab_%s",tag.c_str()),Form("%s #theta_{#pi^{+}#pi^{-}}^{lab};#theta[deg];dN/d#theta",tag.c_str()),200,0,180);
  the_pippim_lab->SetLineColor(4);
  the_pippim_lab->SetLineWidth(2);
  TH1F* phi_pi0_cm = new TH1F(Form("phi_pi0_cm_%s",tag.c_str()),Form("%s #phi_{#pi^{0}}^{cm};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pi0_cm->SetLineColor(2);
  phi_pi0_cm->SetLineWidth(2);
  TH1F* phi_pippim_cm = new TH1F(Form("phi_pippim_cm_%s",tag.c_str()),Form("%s #phi_{#pi^{+}#pi^{-}}^{cm};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pippim_cm->SetLineColor(4);
  phi_pippim_cm->SetLineWidth(2);
  TH1F* phi_pi0_lab = new TH1F(Form("phi_pi0_lab_%s",tag.c_str()),Form("%s #phi_{#pi^{0}}^{lab};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pi0_lab->SetLineColor(2);
  phi_pi0_lab->SetLineWidth(2);
  TH1F* phi_pippim_lab = new TH1F(Form("phi_pippim_lab_%s",tag.c_str()),Form("%s #phi_{#pi^{+}#pi^{-}}^{lab};#phi[deg];dN/d#phi",tag.c_str()),200,0,180);
  phi_pippim_lab->SetLineColor(4);
  phi_pippim_lab->SetLineWidth(2);

  TH1F* minv_pippim_lab = new TH1F(Form("minv_pippim_lab_%s",tag.c_str()),Form("%s minv_pippim_lab; M_{inv} [GeV/c^2]; dN/dM_{inv}",tag.c_str()),200,0,4);
  minv_pippim_lab->SetLineColor(2);
  minv_pippim_lab->SetLineWidth(2);
  TH1F* minv_pippim_cm = new TH1F(Form("minv_pippim_cm_%s",tag.c_str()),Form("%s minv_pippim_cm; M_{inv} [GeV/c^2]; dN/dM_{inv}",tag.c_str()),200,0,4);
  minv_pippim_cm->SetLineColor(4);
  minv_pippim_cm->SetLineWidth(2);

  int Nevt = data_in->GetEntries();
  cout << "nevt= " << Nevt << endl;

  for (int ievt=0; ievt<Nevt; ievt++) {

    if (ievt%10000==0)
      std::cout << "event " << ievt << "/" << Nevt << std::endl;

    data_in->GetEntry(ievt);

    // veto events with more than three particles in final state
    int npart = part_array->GetEntriesFast();
    if (npart!=3) {
      std::cout << "npart= " << npart << std::endl;
    }

    std::vector<int> pdg_codes;
    std::vector<int> indices(3);
    for (int ipart=0; ipart < npart; ++ipart) {
      const TParticle* part = (TParticle*) part_array->At(ipart);
      if (verb) {
	cout << "part= " << part << endl;
	std::cout << " part " << ipart << " pdg= " << part->GetPDG()->PdgCode() << std::endl;
	part->Print();
      }
      int pdg = part->GetPDG()->PdgCode();
      if (pdg==-211) indices[0] = ipart;
      if (pdg==111) indices[1] = ipart;
      if (pdg==211) indices[2] = ipart;
      pdg_codes.push_back(part->GetPDG()->PdgCode());
    }
    std::sort(pdg_codes.begin(),pdg_codes.end());

    if (!compare_sorted(pdg_codes,ref1)) {
      if (verb) print_log(ievt, npart, pdg_codes);
      std::cout << "Event not good" << std::endl;
    }
    type1_count++;

    TLorentzVector p4pi0, p4pip, p4pim, p4pic;
    ((TParticle*)part_array->At(indices[0]))->Momentum(p4pim);
    ((TParticle*)part_array->At(indices[1]))->Momentum(p4pi0);
    ((TParticle*)part_array->At(indices[2]))->Momentum(p4pip);
    p4pic = p4pip + p4pim;
    TLorentzVector cmp4pi0(p4pi0);
    TLorentzVector cmp4pip(p4pip);
    TLorentzVector cmp4pim(p4pim);
    TLorentzVector cmp4pic(p4pic);
    cmp4pi0.Boost(boost_vector);
    cmp4pip.Boost(boost_vector);
    cmp4pim.Boost(boost_vector);
    cmp4pic.Boost(boost_vector);

    p4pic = p4pip + p4pim;
    cmp4pic = cmp4pip + cmp4pim;

    TLorentzVector mom_tot_cm = cmp4pi0 + cmp4pip + cmp4pim;
    TLorentzVector mom_tot = p4pi0 + p4pip + p4pim;

    //cmp4pip.Print();

    //cout << "Mpi0 = " << p4pi0.M() << " Mpip= " << p4pip.M() << " Mpim= " << p4pim.M() << endl;
    //cout << "Minv = " << mom_tot.M() << "  Minv(cm) = " << mom_tot_cm.M() << endl;
    //cout << "Minv(pi+pi-)= " << p4pic.M() << " Minv(cm,pi+pi-)= " << p4pic.M() << endl;
    //cout << "CM Total Mom = "; mom_tot_cm.Print();
    //cout << "LAB Tot Mom= "; mom_tot.Print();
    //
    //cout << "THETA LAB: pi0= " << rtd*p4pi0.Theta() << " pip= " << rtd*p4pip.Theta() << " pim= " << rtd*p4pim.Theta()
    //	 << " pich= " << rtd*p4pic.Theta() << " D0-ch= " << rtd*(p4pi0.Theta()+p4pic.Theta()) <<endl;
    //
    //cout << "THETA  CM: pi0= " << rtd*cmp4pi0.Theta() << " pip= " << rtd*cmp4pip.Theta() << " pim= " << rtd*cmp4pim.Theta()
    //	 << " pich= " << rtd*cmp4pic.Theta() << " D0-ch= " << rtd*(cmp4pi0.Theta()+cmp4pic.Theta()) << endl;
    //
    //cout << "PHI   LAB: pi0= " << rtd*p4pi0.Phi() << " pip= " << rtd*p4pip.Phi() << " pim= " << rtd*p4pim.Phi()
    //	 << " pich= " << rtd*p4pic.Phi() << " D0-ch= " << rtd*(p4pi0.Phi()-p4pic.Phi()) << endl;
    //
    //cout << "PHI    CM: pi0= " << rtd*cmp4pi0.Phi() << " pip= " << rtd*cmp4pip.Phi() << " pim= " << rtd*cmp4pim.Phi()
    //	 << " pich= " << rtd*cmp4pic.Phi() << " D0-ch= " << rtd*(cmp4pi0.Phi()-cmp4pic.Phi()) << endl;

    dthe_cm->Fill( rtd*(cmp4pi0.Theta()+cmp4pic.Theta()) );
    dthe_lab->Fill( rtd*(p4pi0.Theta()+p4pic.Theta()) );
    dphi_cm->Fill( rtd*(cmp4pi0.Phi()-cmp4pic.Phi()) );

    the_pi0_cm->Fill( rtd*cmp4pi0.Theta() );
    the_pippim_cm->Fill( rtd*cmp4pic.Theta() );
    the_pi0_lab->Fill( rtd*p4pi0.Theta() );
    the_pippim_lab->Fill( rtd*p4pic.Theta() );
    phi_pi0_cm->Fill( rtd*cmp4pi0.Phi() );
    phi_pippim_cm->Fill( rtd*cmp4pic.Phi() );
    phi_pi0_lab->Fill( rtd*p4pi0.Phi() );
    phi_pippim_lab->Fill( rtd*p4pic.Phi() );

    minv_pippim_cm->Fill( cmp4pic.M() );
    minv_pippim_lab->Fill( p4pic.M() );

  }
  // summary
  std::cout << "Total Event Count = " << Nevt << std::endl;
  std::cout << "Number of good events = " << type1_count << std::endl;

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->cd();
  dphi_cm->Draw();
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  dthe_cm->Draw();
  //TCanvas *c2 = new TCanvas("c2","c2");
  //dthe_lab->Draw();
  TCanvas *c_minv = new TCanvas("cminv","cminv");
  c_minv->cd();
  minv_pippim_cm->Draw();
  minv_pippim_lab->Draw("same");

  TLegend *tl = new TLegend(0.35,0.65,0.65,0.89);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  TCanvas *c3 = new TCanvas("c3","c3",1400,1000);
  c3->Divide(2,2);

  c3->cd(1);
  the_pi0_cm->SetTitle(Form("%s, %s", the_pi0_cm->GetTitle(), the_pippim_cm->GetTitle() ));
  the_pi0_cm->Draw();
  tl->AddEntry(the_pi0_cm,Form("#pi^{0}, %s",tag.c_str()),"l");
  the_pippim_cm->Draw("same");
  tl->AddEntry(the_pippim_cm,Form("#pi^{+}#pi^{-}, %s",tag.c_str()),"l");
  tl->Draw();

  c3->cd(2);
  gPad->SetLogy();
  the_pippim_lab->SetTitle(Form("%s, %s", the_pi0_lab->GetTitle(), the_pippim_lab->GetTitle() ));
  the_pippim_lab->Draw();
  the_pi0_lab->Draw("same");
  tl->Draw();

  c3->cd(3);
  phi_pi0_cm->SetTitle(Form("%s, %s", phi_pi0_cm->GetTitle(), phi_pippim_cm->GetTitle() ));
  phi_pi0_cm->SetMinimum(0.0);
  phi_pi0_cm->Draw();
  phi_pippim_cm->Draw("same");
  tl->Draw();

  c3->cd(4);
  phi_pi0_lab->SetTitle(Form("%s, %s", phi_pi0_lab->GetTitle(), phi_pippim_lab->GetTitle() ));
  phi_pi0_lab->SetMinimum(0.0);
  phi_pi0_lab->Draw();
  phi_pippim_lab->Draw("same");
  tl->Draw();

  c3->Print(Form("ang_dist_pbarp_pi0pippim_%s.png",tag.c_str()));
  //if (model==0) {
  //  c3->Print("ang_dist_pbarp_pi0pippim_dpm.png");
  //} else {
  //  c3->Print("ang_dist_pbarp_pi0pippim_ftf.png");
  //}

  //TFile *fout = TFile::Open(Form("output_%s.root",fname.c_str()),"RECREATE");
  //fout->cd();
  //dthe_cm->Write();
  //dthe_lab->Write();
  //dphi_cm->Write();
  //the_pi0_cm->Write();
  //the_pippim_cm->Write();
  //the_pi0_lab->Write();
  //the_pippim_lab->Write();
  //phi_pi0_cm->Write();
  //phi_pippim_cm->Write();
  //phi_pi0_lab->Write();
  //phi_pippim_lab->Write();
  //minv_pippim_lab->Write();
  //minv_pippim_cm->Write();







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

bool compare_sorted(const vector<int> &evt, const vector<int> &ref) {
  if (evt.size()!=ref.size()) return false;
  for (unsigned int ievt=0; ievt<evt.size(); ++ievt) {
    if (evt[ievt]!=ref[ievt]) return false;
  }
  return true;
}

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

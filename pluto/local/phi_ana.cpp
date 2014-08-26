#include <iostream>
#include <vector>

// ROOT
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

// Pluto
#include "PParticle.h"
#include "PReaction.h"

using namespace std;

void kinematics(vector<PParticle*>, vector<double> &);
void kinematics(vector<PParticle*>, vector<PParticle*>, vector<double> &);

int main() {

  int verbose = 0;

  TFile *f = new TFile("phi_dalitz.root");
  TTree *d = (TTree*) f->Get("data");

  TCanvas *tc_kin1 = new TCanvas("tc_kin1","tc_kin1",1400,1000);
  tc_kin1->Divide(2,2);
  TH1F* h_oa_kp_km_lab = new TH1F("h_oa_kp_km_lab","K^{+}K^{-} opening angle, Lab",200,0,TMath::Pi());
  TH1F* h_oa_kp_km_cm = new TH1F("h_oa_kp_km_cm","K^{+}K^{-} opening angle, CM",200,0,TMath::Pi());
  TH1F* h_dth_kp_km_cm = new TH1F("h_dth_kp_km_cm","K^{+}K^{-} #Delta#theta, CM",200,-TMath::Pi(),TMath::Pi());
  TH1F* h_phi_kp_km_cm = new TH1F("h_phi_kp_km_cm","K^{+}K^{-} #phi, CM",200,TMath::Pi(),TMath::Pi());
  h_phi_kp_km_cm->SetMinimum(0);

  TCanvas *tc_kin2 = new TCanvas("tc_kin2","tc_kin2",1400,1000);
  tc_kin2->Divide(2,2);
  TH1F* h_mom_kpkm = new TH1F("h_mom_kpkm","p_{tot} of K+ and K-",100,0,4);
  TH1F* h_minv_kp_km1 = new TH1F("h_minv_kp_km1","M_{K^{+}K^{-}}",200,0.9,2.2);
  TH1F* h_minv_kp_km2 = new TH1F("h_minv_kp_km2","M_{K^{+}K^{-}}",200,0.9,2.2);
  TH1F* h_minv_tot = new TH1F("h_minv_tot","M_{K^{+}K^{-}K^{+}K^{-}}",100,3,3.5);

  TClonesArray *part_array = new TClonesArray("PParticle");
  d->SetBranchAddress("Particles",&part_array);

  int nent = d->GetEntries();
  cout << "Nent = " << nent << endl;
  //for (int ient=0; ient<nent; ++ient) {
  for (int ient=0; ient<1; ++ient) {
    d->GetEntry(ient);

    //if (verbose>2||ient%(nent/10)==0)
      cout << "========== new event " << ient << " =============" << endl;
    int npart=part_array->GetEntries();
    if (verbose>2) cout << "npart= " << npart << endl;

    vector<PParticle*> chan1, chan2, chan3;


    TLorentzVector mom4;
    for (int ipart=0; ipart<npart; ++ipart) {

      PParticle *part = (PParticle*) part_array->At(ipart);
      part->SetIndex(ipart);
      if (verbose>3) cout << "part " << ipart << " index = " << part->GetIndex() << " parent index= " << part->GetParentIndex() << endl;
      if ( part->Is("anti_p") || part->Is("phi") )
	chan1.push_back(part);

      if ( part->Is("phi") ) {
	if ( chan2.size()==0 )
	  chan2.push_back(part);
	else
	  chan3.push_back(part);
      }

      if ( part->Is("K+") || part->Is("K-") ) {
	mom4 += part->Vect4();
	h_mom_kpkm->Fill(part->Vect4().Vect().Mag());
	if (chan2.size()>0 && part->GetParentIndex() == chan2[0]->GetIndex() ) chan2.push_back(part);
	else if (chan3.size()>0 && part->GetParentIndex() == chan3[0]->GetIndex() ) chan3.push_back(part);
	else cout << "This must be an ERROR!~~" << endl;
      }

    }
    h_minv_tot->Fill(mom4.M());


    if (verbose>2) {
      cout << "Chan1[0] ids,idx: "; for (unsigned int i=0; i<chan1.size(); ++i) cout << chan1[i]->ID() << "," << chan1[i]->GetIndex() << " "; cout << endl;
      cout << "Chan2[0] ids,idx: "; for (unsigned int i=0; i<chan2.size(); ++i) cout << chan2[i]->ID() << "," << chan2[i]->GetIndex() << " "; cout << endl;
      cout << "Chan3[0] ids,idx: "; for (unsigned int i=0; i<chan3.size(); ++i) cout << chan3[i]->ID() << "," << chan3[i]->GetIndex() << " "; cout << endl;
    }

    int ires = 0;
    vector<double> result;
    kinematics(chan2, result);
    h_oa_kp_km_lab->Fill(result.at(ires++));
    h_minv_kp_km1->Fill(result.at(ires++));
    h_oa_kp_km_cm->Fill(result.at(ires++));
    h_dth_kp_km_cm->Fill(result.at(ires++));
    h_phi_kp_km_cm->Fill(result.at(ires++));


    result.empty();

    kinematics(chan3, result);
    ires = 0;
    h_oa_kp_km_lab->Fill(result.at(ires++));
    h_minv_kp_km2->Fill(result.at(ires++));
    h_oa_kp_km_cm->Fill(result.at(ires++));
    h_dth_kp_km_cm->Fill(result.at(ires++));
    h_phi_kp_km_cm->Fill(result.at(ires++));


    //result.empty();
    //kinematics(chan2,chan3,result);
    //h_minv_tot->Fill(result[0]);

  }

  tc_kin1->cd(1);
  h_oa_kp_km_lab->Draw();
  tc_kin1->cd(2);
  h_oa_kp_km_cm->Draw();
  tc_kin1->cd(3);
  h_dth_kp_km_cm->Draw();
  tc_kin1->cd(4);
  h_phi_kp_km_cm->Draw();
  tc_kin1->Print("kin1.png");

  tc_kin2->cd(1);
  h_mom_kpkm->Draw();
  tc_kin2->cd(2);
  h_minv_kp_km1->Draw();
  tc_kin2->cd(3);
  h_minv_kp_km2->Draw();
  tc_kin2->cd(4);
  h_minv_tot->Draw();
  tc_kin2->Print("kin2.png");

  TFile *fout = TFile::Open("out.root","RECREATE");
  fout->cd();
  tc_kin1->Write();
  tc_kin2->Write();
  fout->Close();

}



void kinematics(vector<PParticle*> chan, vector<double> &res) {

  TLorentzVector p_p = chan[0]->Vect4();
  TLorentzVector p_c1 = chan[1]->Vect4();
  TLorentzVector p_c2 = chan[2]->Vect4();
  res.push_back(p_c1.Angle(p_c2.Vect())); //0

  if (chan[1]->GetParentIndex()!=chan[2]->GetParentIndex())
    cout << "M(K+K-)= " << (p_c1+p_c2).M() << " id1= " << chan[1]->ID() << " id2= " << chan[2]->ID() << " pid1= " << chan[1]->GetParentIndex() << " pid2= " << chan[2]->GetParentIndex() << endl;

  res.push_back( (p_c1+p_c2).M() ); // 1

  p_p = p_p;

  cout << "(boostv)p_p = "; p_p.Print();
  cout << "(boostv)p_p = "; p_p.Vect().Print();
  cout << "(boostv)|p_p| = " << p_p.Vect().Mag() << endl;
  cout << "(before)p_c1= "; p_c1.Print();
  cout << "(before)p_c2= "; p_c2.Print();

  p_c1.Boost(p_p.Vect());
  p_c2.Boost(p_p.Vect());

  cout << "( after)p_c1= "; p_c1.Print();
  cout << "( after)p_c2= "; p_c2.Print();

  res.push_back(p_c1.Angle(p_c2.Vect()));  //2
  //res.push_back(p_c1.Theta() - p_c2.Theta()); //3
  res.push_back(p_c1.Theta()); //3
  //res.push_back(p_c1.Phi() - p_c2.Phi() ); //4
  res.push_back(p_c1.Phi()); //4

  cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
}

void kinematics(vector<PParticle*> chan1, vector<PParticle*> chan2, vector<double> &res) {

  TLorentzVector p_c1_c1 = chan1[1]->Vect4();
  TLorentzVector p_c1_c2 = chan1[2]->Vect4();

  TLorentzVector p_c2_c1 = chan1[1]->Vect4();
  TLorentzVector p_c2_c2 = chan1[2]->Vect4();

  res.push_back( (p_c1_c2+p_c1_c2+p_c1_c2+p_c1_c2).M() );


}

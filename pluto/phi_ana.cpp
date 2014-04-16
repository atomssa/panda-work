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
  
  TFile *f = new TFile("ntuples/phi_dalitz.root");
  TTree *d = (TTree*) f->Get("data");

  TCanvas *tc_kin0 = new TCanvas("tc_kin0","tc_kin0",1400,1000);
  tc_kin0->Divide(2,2);
  TH1F* h_oa_phi_lab = new TH1F("h_oa_phi_lab","#phi-#phi opening angle, Lab",200,0,TMath::Pi());
  TH1F* h_oa_phi_cm = new TH1F("h_oa_phi_cm","#phi-#phi opening angle, CM",100,0.995*TMath::Pi(),1.005*TMath::Pi());
  TH1F* h_costh_phi_cm = new TH1F("h_costh_phi_cm","cos(#theta_{#phi}), CM",200,-1.05,1.05);  
  TH1F* h_phi_phi_cm = new TH1F("h_phi_phi_cm","#phi-#phi #phi, CM",200,TMath::Pi(),TMath::Pi());
  TH1F* h_minv_phi = new TH1F("h_minv_phi","M_{#phi-#phi}",200,0.9,2.2);

  TCanvas *tc_kin1 = new TCanvas("tc_kin1","tc_kin1",1400,1000);
  tc_kin1->Divide(2,2);
  TH1F* h_oa_kpkm_lab = new TH1F("h_oa_kpkm_lab","K^{+}K^{-} opening angle, Lab",200,0,TMath::Pi());
  TH1F* h_oa_kpkm_cm = new TH1F("h_oa_kpkm_cm","K^{+}K^{-} opening angle, CM",100,0.995*TMath::Pi(),1.005*TMath::Pi());
  TH1F* h_costh_kpkm_cm = new TH1F("h_costh_kpkm_cm","cos(#theta_{K+}), CM",200,-1.05,1.05);  
  TH1F* h_phi_kpkm_cm = new TH1F("h_phi_kpkm_cm","K^{+}K^{-} #phi, CM",200,TMath::Pi(),TMath::Pi());
  h_phi_kpkm_cm->SetMinimum(0);
  
  TCanvas *tc_kin2 = new TCanvas("tc_kin2","tc_kin2",1400,1000);
  tc_kin2->Divide(2,2);
  TH1F* h_mom_kpkm = new TH1F("h_mom_kpkm","p of K+ and K-",100,0,4);
  TH1F* h_minv_kpkm1 = new TH1F("h_minv_kpkm1","M_{K^{+}K^{-}}",200,0.9,2.2);
  TH1F* h_minv_kpkm2 = new TH1F("h_minv_kpkm2","M_{K^{+}K^{-}}",200,0.9,2.2);
  TH1F* h_minv_tot = new TH1F("h_minv_tot","M_{K^{+}K^{-}K^{+}K^{-}}",100,3,3.5);
  
  TClonesArray *part_array = new TClonesArray("PParticle");
  d->SetBranchAddress("Particles",&part_array);
  
  int nent = d->GetEntries();
  cout << "Nent = " << nent << endl;
  for (int ient=0; ient<nent; ++ient) {
    //for (int ient=0; ient<1; ++ient) {    
    d->GetEntry(ient);

    if (verbose>2||ient%(nent/10)==0)
      cout << "========== new event " << ient << " =============" << endl;
    int npart=part_array->GetEntries();
    if (verbose>2) cout << "npart= " << npart << endl;

    vector<PParticle*> chan1, chan2, chan3;

    TLorentzVector mom4;
    for (int ipart=0; ipart<npart; ++ipart) {

      PParticle *part = (PParticle*) part_array->At(ipart);
      part->SetIndex(ipart);

      if (verbose>3) cout << "part " << ipart << " index = " << part->GetIndex() << " parent index= " << part->GetParentIndex() << endl;
      if ( part->Is("anti_p") || part->Is("phi") || part->ID()==14015 )
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
      cout << "Chan1[0] ids,idx: "; for (unsigned int i=0; i<chan1.size(); ++i) cout << chan1[i]->ID() << "," << chan1[i]->GetIndex() << "," << chan1[i]->GetParentIndex() << "  "; cout << endl;
      cout << "Chan2[0] ids,idx: "; for (unsigned int i=0; i<chan2.size(); ++i) cout << chan2[i]->ID() << "," << chan2[i]->GetIndex() << "," << chan2[i]->GetParentIndex() << "  "; cout << endl;
      cout << "Chan3[0] ids,idx: "; for (unsigned int i=0; i<chan3.size(); ++i) cout << chan3[i]->ID() << "," << chan3[i]->GetIndex() << "," << chan3[i]->GetParentIndex() << "  "; cout << endl;
    }

    int ires = 0;

    vector<double> result1;
    result1.empty();
    kinematics(chan1,result1);
    h_oa_phi_lab->Fill(result1.at(ires++));
    h_minv_phi->Fill(result1.at(ires++));
    h_oa_phi_cm->Fill(result1.at(ires++));
    h_costh_phi_cm->Fill(result1.at(ires++));
    h_phi_phi_cm->Fill(result1.at(ires++));
    
    vector<double> result2;
    result2.empty();
    kinematics(chan2,result2);
    ires=0;
    h_oa_kpkm_lab->Fill(result2.at(ires++));
    h_minv_kpkm1->Fill(result2.at(ires++));
    h_oa_kpkm_cm->Fill(result2.at(ires++));
    h_costh_kpkm_cm->Fill(result2.at(ires++));
    h_phi_kpkm_cm->Fill(result2.at(ires++));
    
    
    vector<double> result3;    
    result3.empty();
    kinematics(chan3,result3);
    ires = 0;
    h_oa_kpkm_lab->Fill(result3.at(ires++));
    h_minv_kpkm2->Fill(result3.at(ires++));
    h_oa_kpkm_cm->Fill(result3.at(ires++));
    h_costh_kpkm_cm->Fill(result3.at(ires++));
    h_phi_kpkm_cm->Fill(result3.at(ires++));


    //result.empty();
    //kinematics(chan2,chan3,result);
    //h_minv_tot->Fill(result[0]);
    
  }


  tc_kin0->cd(1);
  h_oa_phi_lab->Draw();
  tc_kin0->cd(2);
  h_oa_phi_cm->Draw();
  tc_kin0->cd(3);
  h_costh_phi_cm->Draw();
  tc_kin0->cd(4);
  h_phi_phi_cm->Draw();  
  tc_kin0->Print("kin0.png");

  tc_kin1->cd(1);
  h_oa_kpkm_lab->Draw();
  tc_kin1->cd(2);
  h_oa_kpkm_cm->Draw();
  tc_kin1->cd(3);
  h_costh_kpkm_cm->Draw();
  tc_kin1->cd(4);
  h_phi_kpkm_cm->Draw();  
  tc_kin1->Print("kin1.png");

  tc_kin2->cd(1);
  h_mom_kpkm->Draw();
  tc_kin2->cd(2);
  h_minv_kpkm1->Draw();
  tc_kin2->cd(3);
  h_minv_kpkm2->Draw();
  tc_kin2->cd(4);
  h_minv_tot->Draw();
  tc_kin2->Print("kin2.png");

  TFile *fout = TFile::Open("out.root","RECREATE");
  fout->cd();
  tc_kin0->Write();
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
    cout << "M(K+K-)= " << (p_c1+p_c2).M() << " id1= " << chan[1]->ID() << " id2= " << chan[2]->ID() 
	 << " pid1= " << chan[1]->GetParentIndex() << " pid2= " << chan[2]->GetParentIndex() << endl;

  res.push_back( (p_c1+p_c2).M() ); // 1

  p_c1.Boost(-(p_p.BoostVector()));
  p_c2.Boost(-(p_p.BoostVector()));

  res.push_back( p_c1.Angle(p_c2.Vect()) );  //2
  //res.push_back( p_c1.Theta() - p_c2.Theta() ); //3
  res.push_back(TMath::Cos(p_c1.Theta())); //3  
  //res.push_back( p_c1.Phi() - p_c2.Phi() ); //4
  res.push_back(p_c1.Phi()); //4

}

void kinematics(vector<PParticle*> chan1, vector<PParticle*> chan2, vector<double> &res) {

  TLorentzVector p_c1_c1 = chan1[1]->Vect4();
  TLorentzVector p_c1_c2 = chan1[2]->Vect4();
  
  TLorentzVector p_c2_c1 = chan1[1]->Vect4();
  TLorentzVector p_c2_c2 = chan1[2]->Vect4();

  res.push_back( (p_c1_c2+p_c1_c2+p_c1_c2+p_c1_c2).M() );
  
  
}

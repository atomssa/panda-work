//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: EvtPpbarJpsiPi0.cpp,v 1.2 2008/01/08 17:15:11 steinke Exp $
//
// Description:
//            Generator of pbar p -> jpsi + pi0
//
// Author List:
// 
//
//------------------------------------------------------------------------
#include "EvtGenModels/EvtPpbarJpsiPi0.hh"

#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtStdHep.hh"

// Temporary to store all events regardless of acceptance
//#include "TFile.h"
//#include "TTree.h"
//#include "TClonesArray.h"
//#include "TParticle.h"
//#include "TParticlePDG.h"
//#include "TLorentzVector.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TClass.h"
#include "TParticlePDG.h"

//using std::fstream;
using std::string;
using namespace std;

EvtPpbarJpsiPi0::~EvtPpbarJpsiPi0()
{
  fFile->cd();
  fTree->Write();
  fFile->Write();
  fFile->Close();
}

std::string EvtPpbarJpsiPi0::getName()
{
  return "PpbarJpsiPi0";
}

EvtDecayBase* EvtPpbarJpsiPi0::clone()
{
  return new EvtPpbarJpsiPi0;
}

void EvtPpbarJpsiPi0::initProbMax()
{
  setProbMax(500);
}

void EvtPpbarJpsiPi0::init()
{
  checkNDaug(2);
  checkNArg(1);
  _s = getArg(0); 

  // Temporary to store all events regardless of acceptance   
  fFile = new TFile("Signal-micro.root","RECREATE","ROOT_Tree"); 
  cout << "init_root_tree() : TParticle class: " << TClass::GetClass("TParticle") << endl;
  fEvt = new TClonesArray("TParticle",100);
  fEvtTlv = new TClonesArray("TLorentzVector",100);
  fTree = new TTree("data","TDA Signal");
  fTree->Branch("prob",&fProb,"prob/F");
  fTree->Branch("Particles",&fEvt, 32000, 2);
  fTree->Branch("PartTlvs",&fEvtTlv, 32000, 2);
}

void EvtPpbarJpsiPi0::decay(EvtParticle* p)
{
  p->initializePhaseSpace(getNDaug(),getDaugs());
   
  static const EvtId pi0ID = EvtPDL::getId("pi0");
  EvtParticle* pi0 =  p->getDaug(0);
  
  if(pi0->getId() != pi0ID )
    pi0= p->getDaug(1);
  if(pi0->getId() != pi0ID )
    cout << "EvtPpbarJpsiPi0::decay():\n wrong id of produced particles!"<<endl;
  
  EvtVector4R p4 = pi0->getP4Lab();
  double _Mp = 0.938;
  double _Mpi0 = 0.135;
  double _Mjpsi2 = 9.59;

  double pi = 3.14159;
  double fpi = 0.093;// GeV
  double gpiNN = 13;
  double M0 = 8776.88; //GeV
  double fphi = 0.413; //GeV
  double fN = 0.005;//GeV
  double alpha_s = 0.25;
  double M = 3.;//GeV
  double C = pow(4.*pi*alpha_s,3)*(fN*fN*fphi/fpi)*(10./81.);

  double _Ebeam = (_s -2.*pow(_Mp,2))/(2.*_Mp);//5.59(s=12.25)   9.72(s=20)
  double _Pbeam = sqrt(pow(_Ebeam,2)-pow(_Mp,2));//5.51(s=12.25)    9.68(s=20)
  double _t = pow((_Ebeam-p4.get(0)),2)-pow((0-p4.get(1)),2)-pow((0-p4.get(2)),2)-pow((_Pbeam-p4.get(3)),2);
  //double _qsit = (_Mjpsi2 - _t - pow(_Mp,2))/(2.*_s-_Mjpsi2+_t-3.*pow(_Mp,2));
  double _qsi = M*M/(2*_s-M*M);
  double _deltaT2 = ((1-_qsi)/(1+_qsi))*(_t - 2.*_qsi*(pow(_Mp,2)/(1+_qsi)-pow(_Mpi0,2)/(1-_qsi)));
  double t_max = 2*_qsi*(_Mp*_Mp*(_qsi-1)+_Mpi0*_Mpi0*(_qsi+1))/(_qsi*_qsi-1);
  double x = _t;
  
  double prob;
  if (x<-0.5 || x >= t_max)
    prob =0;
  else
    {
      double I = M0*fpi*gpiNN*_Mp*(1-_qsi)/(x-_Mp*_Mp)/(1+_qsi);
      double Iprim = M0*fpi*gpiNN*_Mp/(x-_Mp*_Mp);
      double MT2 = 0.25*C*C*(2*(_qsi+1)/_qsi/pow(M,8))*(I*I-_deltaT2*Iprim*Iprim/pow(_Mp,2));
      double lamda2 = pow(_s,2)-4*_s*pow(_Mp,2);
      prob = MT2/(16*pi*lamda2)*0.38941*pow(10,9);
    }

  setProb(prob);

  // Temporary to store all events regardless of acceptance 
  bool verb_flag = true;
  fProb = prob;
  fEvtStdHep.init();
  p->makeStdHep(fEvtStdHep);
  fEvt->Clear();
  int tca_idx = 0;
  for(Int_t iPart=0; iPart<fEvtStdHep.getNPart(); iPart++){
    // add track
    //const Int_t nFD = fEvtStdHep.getFirstDaughter(iPart);
    //const Int_t nLD = fEvtStdHep.getLastDaughter(iPart);
    //if( nFD==-1 && nLD==-1 )
    //{
	const int fId = fEvtStdHep.getStdHepID(iPart);
	TLorentzVector _p4, _v4;
	set_p4_v4(fEvtStdHep, iPart, _p4, _v4);
	// onto the free store for the TCA
	//TParticle *tmp = new((*fEvt)[tca_idx++]) TParticle(fId,1,0,0,0,0,_p4,_v4);
	TLorentzVector *tmp = new((*fEvtTlv)[tca_idx++]) TLorentzVector();
	tmp->SetPxPyPzE(_p4.Px(),_p4.Py(),_p4.Pz(),_p4.E());
	if(verb_flag) {
	  cout << "- I -: New particle " << iPart << " Id= " << fId << " Prob= " << fProb << " -> ";
	  tmp->Print();
	}
	//}
  }
  fTree->Fill();
  if (verb_flag) cout <<"==== EvtPpbarJpsiPi0 ==="<<endl;
  
  return;
}

//_________ Simple helper to convert 4-mom and 4-vtx from EvtVector4R to TLorentzVector __________________________
void EvtPpbarJpsiPi0::set_p4_v4(EvtStdHep &part, const int &ipart, TLorentzVector &_p4, TLorentzVector &_v4) {

  const int fId = part.getStdHepID(ipart);
  const EvtVector4R v4 = part.getX4(ipart);
  const EvtVector4R p4 = part.getP4(ipart);

  const Double_t mm2cm = 0.1;
  _p4.SetPxPyPzE(p4.get(1), p4.get(2), p4.get(3), p4.get(0));
  _v4.SetPxPyPzE(v4.get(1)*mm2cm, v4.get(2)*mm2cm, v4.get(3)*mm2cm, v4.get(0));
  //_p4.Print();
}


































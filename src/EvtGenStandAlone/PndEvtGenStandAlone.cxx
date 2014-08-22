// -------------------------------------------------------------------------
// -----                PndEvtGenStandAlone source file                  -----
// -----                                                          -----
// -------------------------------------------------------------------------


// FAIR headers
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

#include "FairTask.h"

#include <iostream>
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
#include "PndEvtGenStandAlone.h"

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"

/// use of external generators
#if EVTGEN_EXTERNAL
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#endif

#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

using std::endl;
using std::ofstream;
using std::cout;

//define class for generating random nubers
class EvtRootRandomEngine:public EvtRandomEngine{
public:
  EvtRootRandomEngine(int s=0) {seed=s;}
  double random();
  int seed;
};

double EvtRootRandomEngine::random(){
  static TRandom3 randengine(seed);
  return randengine.Rndm();
}


// -----   Default constructor   ------------------------------------------
PndEvtGenStandAlone::PndEvtGenStandAlone() {
  SetName("PndEvtGenStandAlone");
  fVerb=0;
}
// ------------------------------------------------------------------------

// -----   Standard constructor   -----------------------------------------
PndEvtGenStandAlone::PndEvtGenStandAlone(TString particle,TString decfile,Double_t Mom, Long_t Seed,TString defaultDECAY,TString defaultPDL) {

  PndEvtGenStandAlone();

  cout << "<I> PndEvtGenStandAlone"<<endl;
  cout << "<I> Particle: "<<particle<<endl;
  cout << "<I> decfile: "<<decfile<<endl;
  if(Mom>0) cout << "<I> pbar-Momentum: "<<Mom<<endl;
  if(Mom==0) cout << "<I> Momentum: "<<Mom<<endl;
  if(Mom<0) cout << "<I> CMS energy: "<<Mom<<endl;
  cout << "<I> Rnd Seed: "<<Seed<<endl;
  TString work = getenv("VMCWORKDIR");
  if (defaultDECAY=="") defaultDECAY = work + "/pgenerators/EvtGen/EvtGen/Private/DECAY.DEC";
  if (defaultPDL=="") defaultPDL = work + "/pgenerators/EvtGen/EvtGen/Private/evt.pdl";

  //Initialize the generator - read in the decay table and particle properties
  EvtRandomEngine* myRandomEngine=0;

  // Make sure that the seed is always set if default value (-1) is used, JGM, August 2011
  if (Seed<0) 
    {
      Seed = gRandom->GetSeed();
      cout << "<I> Rnd Seed changed to " << Seed << endl;
    } 
  myRandomEngine=new EvtRootRandomEngine(Seed);

  // Set up the default external generator list: Photos, Pythia and/or Tauola ONLY if available
#if EVTGEN_EXTERNAL
  EvtExternalGenList genList;
  EvtAbsRadCorr* radCorrEngine = genList.getPhotosModel();
  std::list<EvtDecayBase*> extraModels = genList.getListOfModels();

  // Create the EvtGen generator object
  fGenerator=new EvtGen(defaultDECAY,defaultPDL,myRandomEngine,
			 radCorrEngine, &extraModels);
#else
  //If you don't want to use external generators, use the following:
  fGenerator=new EvtGen(defaultDECAY,defaultPDL,myRandomEngine);
#endif

  //If I wanted a user decay file, I would read it in now.
  if(decfile!="") fGenerator->readUDecay(decfile.Data());

  PART=EvtPDL::getId(std::string(particle.Data()));

  if ( (particle.Contains("pbarp") || particle.Contains("pbard")) && Mom==0)
    {
      cerr <<"\033[5m\033[31m -E  ******  FATAL ERROR: <particle> is '" << particle.Data() << "'; MUST give pbar momentum or cms energy!\033[0m"<<endl;
      exit(0);
    }
 
  double val=-3.0969;
  fMomentum = 0.0;
  fEnergy = 0.0;
  double mp=0.93827;
  double md=1.875613;

  if ( (particle.Contains("pbarp") || particle.Contains("pbard")) && Mom!=0){
    val=Mom;
  }else{
    if(PART.getId()==-1){
      cerr << "Particle \""<<particle<<"\" is unknown!!!"<<endl<<"Check your Macro for spelling mistake."<<endl;
      exit(0);
    }
    val=-EvtPDL::getMass(PART);
  }
  
  // val is the momentum of the pbar beam
  if (val>0){  
    fMomentum = val;
    if ( particle.Contains("pbarpSystem") ) fEnergy = mp+sqrt(fMomentum*fMomentum+mp*mp);
    if ( particle.Contains("pbardSystem") ) fEnergy = md+sqrt(fMomentum*fMomentum+mp*mp);
  }
  else  //val is -E_cm
    {
      val=-val;
      fEnergy = val*val/(2*mp);
      fMomentum = sqrt(fEnergy*fEnergy-val*val);
    }
  
  cout <<"\n############# Generating with following conditions:\n\n";
  cout <<"incident 4-mom : ("<<fEnergy<<", 0, 0, "<<fMomentum<<"), m = "<<sqrt(fEnergy*fEnergy-fMomentum*fMomentum)<<endl;
  cout <<"\n######################\n\n"<<endl;

}
// ------------------------------------------------------------------------

// -----   Destructor   ---------------------------------------------------
PndEvtGenStandAlone::~PndEvtGenStandAlone() {
}


void PndEvtGenStandAlone::init_root_tree() {

  iEvt = 0;
  //   Root initialization 
  fFile = new TFile("Signal-nano.root","RECREATE","ROOT_Tree"); 

  fNpart=0;
  cout << "init_root_tree() : TParticle class: " << TClass::GetClass("TParticle") << endl;
  fEvt = new TClonesArray("TParticle",100);

  fTree = new TTree("data","TDA Signal");
  fTree->Branch("Npart",&fNpart,"Npart/I");
  fTree->Branch("Particles",&fEvt, 32000, 0);

}

// ------------------------------------------------------------------------
InitStatus PndEvtGenStandAlone::Init() {
  Initialize();
  return kSUCCESS;
}

// ------------------------------------------------------------------------
void PndEvtGenStandAlone::Initialize() {
  init_root_tree();
}

// -----  Execute Event Loop ----------------
void PndEvtGenStandAlone::Exec() {
  generate_event(fEvt);
  iEvt++;
}

// -----  Execute Event Loop ----------------
void PndEvtGenStandAlone::Exec(TClonesArray *evt) {
  generate_event(evt);
  iEvt++;
}

// -----  Execute Event Loop ----------------
void PndEvtGenStandAlone::Exec(Option_t* opt) {
  Exec();
}

// -----  Do finish (automatically callecd by FairTask)
void PndEvtGenStandAlone::Finish() {
  Finalize();
}

//----- Do finish (if not called as FairTask)
void PndEvtGenStandAlone::Finalize() {
  close_root_file();
}

// -----  Save tree to output file --------------
void PndEvtGenStandAlone::close_root_file() {
  fFile->cd();
  fTree->Write();
  fFile->Write();  
  fFile->Close();
}

//// -----   Public method ReadEvent   --------------------------------------
//Bool_t PndEvtGenStandAlone::generate_event() {
//  
//  if (iEvt%1000==0) cout << "PndEvtGenStandAlone::generate_event generating event " << iEvt << endl;
//
//  // Loop to create nEvents, starting from an Upsilon(4S)
//  // Set up the parent particle
//  EvtParticle *parent;
//  EvtVector4R pInit(fEnergy,  0.0000, -0.0000,  fMomentum);
//  parent=EvtParticleFactory::particleFactory(PART,pInit);
//  parent->setDiagonalSpinDensity();  
//
//  fGenerator->generateDecay(parent);
//
//  fEvtStdHep.init();
//  parent->makeStdHep(fEvtStdHep);
//
//  if (verb_flag()) print_detail(parent);
//    
//  // Write the output
//  fNpart = fEvtStdHep.getNPart();
//
//  fEvt->Clear();
//  for(Int_t iPart=0; iPart<fNpart; iPart++){
//
//    // add track
//    const Int_t nFD = fEvtStdHep.getFirstDaughter(iPart);
//    const Int_t nLD = fEvtStdHep.getLastDaughter(iPart);
//
//    if( nFD==-1 && nLD==-1 )
//      {
//
//	const int fId = fEvtStdHep.getStdHepID(iPart);
//	TLorentzVector p4, v4;
//	set_p4_v4(fEvtStdHep, iPart, p4, v4);
//
//	// onto the free store for the TCA
//	TParticle *tmp = new((*fEvt)[iPart]) TParticle(fId,1,0,0,0,0,p4,v4);
//	if(verb_flag()) {
//	  cout << "- I -: New particle " << iPart << " -> " << endl;
//	  tmp->Print();
//	}
//	
//      }
//
//  }
//
//  if (!compare()) print_detail(parent);
//  
//  if (verb_flag()) cout <<"==== compare end ==="<<endl;
//
//  fTree->Fill();
//
//  parent->deleteTree();  
//
//  return kTRUE;
//
//}

// -----   Public method ReadEvent   --------------------------------------
Bool_t PndEvtGenStandAlone::generate_event(TClonesArray*evt) {
  
  if (iEvt%1000==0) cout << "PndEvtGenStandAlone::generate_event generating event " << iEvt << endl;

  cout << "==================================================== From code:" << endl;

  // Loop to create nEvents, starting from an Upsilon(4S)
  // Set up the parent particle
  EvtParticle *parent;
  EvtVector4R pInit(fEnergy,  0.0000, -0.0000,  fMomentum);
  parent=EvtParticleFactory::particleFactory(PART,pInit);
  parent->setDiagonalSpinDensity();  

  fGenerator->generateDecay(parent);

  fEvtStdHep.init();
  parent->makeStdHep(fEvtStdHep);

  if (verb_flag()) print_detail(parent);
    
  // Write the output
  fNpart = fEvtStdHep.getNPart();

  int tca_idx = 0;
  for(Int_t iPart=0; iPart<fNpart; iPart++){

    // add track
    const Int_t nFD = fEvtStdHep.getFirstDaughter(iPart);
    const Int_t nLD = fEvtStdHep.getLastDaughter(iPart);

    if( nFD==-1 && nLD==-1 )
      {

	const int fId = fEvtStdHep.getStdHepID(iPart);
	TLorentzVector p4, v4;
	set_p4_v4(fEvtStdHep, iPart, p4, v4);

	// onto the free store for the TCA
	TParticle *tmp = new((*evt)[tca_idx++]) TParticle(fId,1,0,0,0,0,p4,v4);
	//if(verb_flag()) {
	cout << "- I -: New particle " << iPart << " -> " << endl;
	tmp->Print();
	//}

      }
  }

  evt->Print();

  if (!compare()) print_detail(parent);
  
  if (verb_flag()) cout <<"==== compare end ==="<<endl;

  fTree->Fill();

  parent->deleteTree();  

  return kTRUE;

}

Bool_t PndEvtGenStandAlone::compare() {
  if (ref_topo.size() == 0 ) return true; // someone forgot to set the reference topo
  const int size = fEvt->GetEntriesFast();
  if (size!=ref_topo.size()) return false;
  std::vector<Int_t> evt_topo(size);
  for (unsigned int ipart=0; ipart<size; ++ipart) {
    evt_topo.push_back( ((TParticle*) fEvt->At(ipart))->GetPDG()->PdgCode() );
  }
  std::sort(evt_topo.begin(),evt_topo.end());
  for (unsigned int ipart=0; ipart<evt_topo.size(); ++ipart) {
    if (evt_topo[ipart]!=ref_topo[ipart]) return false;
  }
  return true;
}

inline Bool_t PndEvtGenStandAlone::verb_flag() {
  return fVerb>1 || ( fVerb==1 && (iEvt<10||iEvt%100==0) );
}

//_________ Simple helper to convert 4-mom and 4-vtx from EvtVector4R to TLorentzVector __________________________
void PndEvtGenStandAlone::set_p4_v4(EvtStdHep &part, const int &ipart, TLorentzVector &_p4, TLorentzVector &_v4) {

  const int fId = part.getStdHepID(ipart);
  const EvtVector4R v4 = part.getX4(ipart);
  const EvtVector4R p4 = part.getP4(ipart);

  const Double_t mm2cm = 0.1;
  _p4.SetPxPyPzE(p4.get(1), p4.get(2), p4.get(3), p4.get(0));
  _v4.SetPxPyPzE(v4.get(1)*mm2cm, v4.get(2)*mm2cm, v4.get(3)*mm2cm, v4.get(0));
  //_p4.Print();
}

//_____ Function to print event detail every now and then _________
void PndEvtGenStandAlone::print_detail(EvtParticle *parent) {

  cout << "PndEvtGenStandAlone::ReadEvent "<<iEvt <<" "<<fEnergy<<" "<<fMomentum << endl;
  parent->printParticle();
  // Write out the results
  EvtHepMCEvent theEvent;
  theEvent.constructEvent(parent);
  HepMC::GenEvent* genEvent = theEvent.getEvent();
  genEvent->print(std::cout);
  cout <<"==== now compare ==="<<endl;

}

// ------------------------------------------------------------------------
ClassImp(PndEvtGenStandAlone)


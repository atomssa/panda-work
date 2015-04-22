void anatup( Double_t pBeam=3.3077729, Int_t Did=19)
// Analyse ee events, and write tuples
{

  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();

  FILE *fp;
  Int_t NTmax=100000;
  //   Int_t NFmax=10;   // max number of files
  //   Int_t NTmax=100;   // max number per file
  Int_t NEVcount=0;   // max number per file
  //   Int_t NTmax=10000;
  //   Int_t NTmax=200;
  Int_t ip=0;
  Int_t iPart=3;

  Double_t mProt= 0.938272;
  Double_t mElec= 0.000511;
  Double_t mMuon= 0.105658;
  Double_t mPion= 0.139570;
  Double_t mZero= 0.134977;
  Double_t mZeroFit= 0.13;
  //  Double_t mZeroCut= 0.06;
  Double_t mZeroCut= 0.02;
  //  Double_t pBeam=4.00000;
  Double_t eBeam=sqrt(pBeam*pBeam+mProt*mProt);
  Double_t eSystem=eBeam+mProt;
  Double_t mElec2 = mElec*mElec;
  Double_t mZero2= mZero*mZero;
  Double_t radeg=180./3.1415926535;

  // electroncuts
  //   Double_t EmcCUT=0.8;
  //   Double_t SttCUT=0.5;
  Double_t EmcCUT=0.5;
  Double_t SttCUT=0.5;
  Double_t DiscCUT=0.5;
  Double_t DrcCUT=0.5;
  Double_t MuonCUT=0.5;
  Double_t MvdCUT=0.5;
  Double_t TOTCUT=0.50;
  Double_t TOT90CUT=0.90;
  Double_t TOT95CUT=0.95;
  Double_t TOT98CUT=0.98;
  Double_t TOT99CUT=0.99;

  // Set up the Lorentzvectors of the system
  TLorentzVector target(0.0, 0.0, 0.0, mProt);
  Double_t eBeam=sqrt(pBeam*pBeam+mProt*mProt);
  TLorentzVector beam(0.0, 0.0, pBeam, eBeam);
  TLorentzVector W = beam + target;
  TVector3 bSigma(0,0,-W.Pz()/W.E());


  // Construct filenames

  TString Directory[]={"rfiles3/","rootfiles/"      // vandewie
		       ,"ee01/","ee02/","ee03/","ee04/","ee05/"    // gosia + photos
		       ,"ee06/","ee07/","ee08/","ee09/","ee10/"
		       ,"ee11/","ee12/","ee13/","ee14/","ee15/"
		       ,"ee16/","ee17/"
		       ,"eeno01/","eeno02/","eeno03/","eeno04/"    // gosia + nophotos
  };

  TString name = "_complete.root";
  TString inRecoFile   = Directory[Did]+"reco"+name;
  TString inDigiFile   = Directory[Did]+"digi"+name;
  TString inSimFile    = Directory[Did]+"sim"+name;
  TString inPidFile    = Directory[Did]+"pid"+name;

  TString outAnaFile   = Directory[Did]+"ana"+name;
  TFile *out = TFile::Open(outAnaFile,"RECREATE");

  TNtuple* NTev = new TNtuple("NTev","NTev",
			      "p1MC:th1MC:ph1MC:p2MC:th2MC:ph2MC:costhMC:Q2:p1:th1:ph1:EMC1:NX1:pEMC1:pSTT1:pDIS1:pDRC1:pMVD1:pCOM1:p2:th2:ph2:EMC2:NX2:pEMC2:pSTT2:pDIS2:pDRC2:pMVD2:pCOM2:bestTH:bestPH:bestCOST"
			      ,68000);
  Float_t atuple[33];

  cout << "filename:" << inPidFile << endl;
  TFile *inFile = TFile::Open(inPidFile,"READ");
  TTree *lhe=(TTree *) inFile->Get("cbmsim") ;
  TFile *outFile = TFile::Open(outAnaFile,"new");

  // adding other files as friends
  lhe->AddFriend("cbmsim",inSimFile);
  lhe->AddFriend("cbmsim",inDigiFile);

  PndEmcMapper::Init(6);
  // get the data  (correspondage with Inspect in TBrowser)


  TClonesArray* cCand_array=new TClonesArray("PndPidCandidate");
  lhe->SetBranchAddress("PidChargedCand", &cCand_array);

  TClonesArray* nCand_array=new TClonesArray("PndPidCandidate");
  lhe->SetBranchAddress("PidNeutralCand", &nCand_array);

  TClonesArray* mc_array=new TClonesArray("PndMCTrack");
  lhe->SetBranchAddress("MCTrack", &mc_array);

  TClonesArray* stthit_array=new TClonesArray("PndSttHit");
  lhe->SetBranchAddress("STTHit", &stthit_array);

  TClonesArray* sttpoint_array=new TClonesArray("PndSttPoint");
  lhe->SetBranchAddress("STTPoint", &sttpoint_array);


  TClonesArray* cluster_array=new TClonesArray("PndEmcCluster");
  lhe->SetBranchAddress("EmcCluster",&cluster_array);

  TClonesArray* digi_array=new TClonesArray("PndEmcSharedDigi");
  lhe->SetBranchAddress("EmcSharedDigi",&digi_array);

  TClonesArray* bump_array=new TClonesArray("PndEmcBump");
  lhe->SetBranchAddress("EmcBump",&bump_array);

  TClonesArray* drc_array=new TClonesArray("PndPidProbability");
  lhe->SetBranchAddress("PidAlgoDrc", &drc_array);

  TClonesArray* disc_array=new TClonesArray("PndPidProbability");
  lhe->SetBranchAddress("PidAlgoDisc", &disc_array);

  TClonesArray* mvd_array=new TClonesArray("PndPidProbability");
  lhe->SetBranchAddress("PidAlgoMvd", &mvd_array);

  TClonesArray* stt_array=new TClonesArray("PndPidProbability");
  lhe->SetBranchAddress("PidAlgoStt", &stt_array);

  TClonesArray* emcb_array=new TClonesArray("PndPidProbability");
  lhe->SetBranchAddress("PidAlgoEmcBayes", &emcb_array);


  // canvas and stuff
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetPalette(1);

  int off=32;
  int start=250;


  cMA = new TCanvas("cMA","cMA",200,0, 1000, 1000);
  cMA->Divide(3,3);

  //  cMCelec1 = new TCanvas("cMCelec1","cMCelec1",250,0, 1200, 1000);
  //  cMCelec1->Divide(3,2);

  // histos
  TH1F *h_costheta_mc = new TH1F("h_costheta_mc","cos_theta_ep",20,-1,1);
  TH1F *h_costheta = new TH1F("h_costheta","cos_theta_ep",20,-1,1);
  TH1F *h_costheta_sel = new TH1F("h_costheta_sel","cos_theta_ep",20,-1,1);
  TH1F *h_efficiency = new TH1F("h_efficiency","efficiency",20,-1,1);
  TH1F *h_theta = new TH1F("h_theta","theta_ep",400,-400,400);
  TH1F *h_dtheta = new TH1F("h_dtheta","dtheta_ep",100,170,190);
  TH1F *h_dphi = new TH1F("h_dphi","dphi_ep",100,170,190);

  TH1D* hcutE   = new TH1D("hcutE"  ,"hcutE",50,-eSystem,eSystem);
  TH1D* hcutx   = new TH1D("hcutx"  ,"hcutx",50,-1,1);
  TH1D* hcuty   = new TH1D("hcuty"  ,"hcuty",50,-1,1);
  TH1D* hcutz   = new TH1D("hcutz"  ,"hcutz",50,-1,1);
  TH1D* hMCz    = new TH1D("hMCz"  ,"hMCz",50,0,eSystem);
  TH1D* hpt2    = new TH1D("hpt2"  ,"hpt2",50,0,4);
  TH1D* hpair   = new TH1D("hpair"  ,"hpair",11,-0.5,10.5);

  cout << " finished histos " << endl;

  // process the data

  // Lorentz vectors MV
  TLorentzVector mcTrack[4];
  TLorentzVector QQ_MC;

  // loop over  events
  Double_t MCEnergy, MCTheta, MCPhi;
  Double_t MC1Energy, MC1Theta, MC1Phi;
  Double_t MC2Energy, MC2Theta, MC2Phi;


  NTevents=lhe->GetEntriesFast();
  cout << "NTevents: " << NTevents << endl;
  if(NTevents>NTmax) NTevents=NTmax;
  for (Int_t j=0; j< NTevents ; j++) {

    lhe->GetEntry(j);   // kinematics
    NEVcount++;
    if(j%1000 == 0) cout << "event: " << j << endl;

    //       cout << "processing event: " << j<< endl ;
    Int_t nmc      = mc_array->GetEntriesFast();
    Int_t ncCand=cCand_array->GetEntriesFast();
    Int_t nnCand=nCand_array->GetEntriesFast();
    if(j<5) {
      cout << " nMC: " << mc_array->GetEntriesFast() ;
      cout << " ncCand: " << cCand_array->GetEntriesFast() << endl;
      cout << " nnCand: " << nCand_array->GetEntriesFast() << endl;
    }

    // Loop over electron tracks, store in mcTrack[0,1], pizero_MC, QQ_MC
    Float_t mc_pp = 0, mc_E = 0, mc_mass = 0;
    Float_t mc_px = 0, mc_py = 0, mc_pz = 0;
    Int_t mc_pdg;
    Int_t nepi=0;
    Int_t mc0=0, mc1=1;
    if(Did>1) {mc0=1; mc1=2;}
    // positron
    PndMCTrack *mctrack = (PndMCTrack*)mc_array->At(mc0);
    //       Int_t MotherID = mctrack->GetMotherID();
    //       mc_pdg = (int) (mctrack->GetPdgCode());

    mc_pp = mctrack->GetMomentum().Mag();
    mc_px = mctrack->GetMomentum().Px();
    mc_py = mctrack->GetMomentum().Py();
    mc_pz = mctrack->GetMomentum().Pz();
    mc_E=TMath::Sqrt(mc_pp*mc_pp+mElec2);
    mcTrack[0].SetPxPyPzE(mc_px,mc_py,mc_pz,mc_E);
    // electron
    PndMCTrack *mctrack = (PndMCTrack*)mc_array->At(mc1);
    //       Int_t MotherID = mctrack->GetMotherID();
    //       mc_pdg = (int) (mctrack->GetPdgCode());

    mc_pp = mctrack->GetMomentum().Mag();
    mc_px = mctrack->GetMomentum().Px();
    mc_py = mctrack->GetMomentum().Py();
    mc_pz = mctrack->GetMomentum().Pz();
    mc_E=TMath::Sqrt(mc_pp*mc_pp+mElec2);
    mcTrack[1].SetPxPyPzE(mc_px,mc_py,mc_pz,mc_E);

    // elec1
    MC1Energy = mcTrack[0].E();
    MC1Theta  = radeg*(mcTrack[0].Theta());
    MC1Phi    = radeg*(mcTrack[0].Phi());
    //       hmc1E->Fill(MC1Energy);
    //       hmc1TH->Fill(MC1Theta);
    //       hmc1PH->Fill(MC1Phi);
    // elec2
    MC2Energy = mcTrack[1].E();
    MC2Theta  = radeg*(mcTrack[1].Theta());
    MC2Phi    = radeg*(mcTrack[1].Phi());
    //       hmc1E->Fill(MC2Energy);
    //       hmc1TH->Fill(MC2Theta);
    //       hmc1PH->Fill(MC2Phi);

    // quadrivecteur
    double elecMC=mcTrack[0].Angle(mcTrack[1].Vect());
    //       double Q2=4*mcTrack[0].E()*mcTrack[1].E()*(sin(elecMC/2))*(sin(elecMC/2));
    QQ_MC=mcTrack[0]+mcTrack[1];
    double Q2=QQ_MC.M2();
    //       hmcQ2->Fill(Q2);
    TLorentzVector *elecMC0 = new TLorentzVector(mcTrack[0].Px(),
						 mcTrack[0].Py(),
						 mcTrack[0].Pz(),
						 mcTrack[0].E());
    TLorentzVector *elecMC1 = new TLorentzVector(mcTrack[1].Px(),
						 mcTrack[1].Py(),
						 mcTrack[1].Pz(),
						 mcTrack[1].E());

    // boost
    elecMC0->Boost(bSigma);
    elecMC1->Boost(bSigma);
    Double_t THMC0=radeg*(elecMC0->Theta());
    Double_t THMC1=radeg*(elecMC1->Theta());
    Double_t THtot=THMC0+THMC1;
    Double_t COSTMC=TMath::Cos(elecMC0->Theta());
    h_costheta_mc-> Fill(COSTMC);
    // print some event data
    if(j<5) {
      cout << "event:" << j << endl;
      cout << "MC-eepi:" << MC1Energy << " " << MC2Energy << endl;
      cout << "theta:" << MC1Theta << " " << MC2Theta  << endl;
      cout << "phi:" << MC1Phi << " " << MC2Phi << endl;
      cout << "THMC0:" << THMC0 << "THMC1:" << THMC1 << endl;
      cout << "THtot:" << THtot << " Q2:" << Q2 << endl;
    }

    // fill tuple
    atuple[0]=MC1Energy;
    atuple[1]=MC1Theta;
    atuple[2]=MC1Phi;
    atuple[3]=MC2Energy;
    atuple[4]=MC2Theta;
    atuple[5]=MC2Phi;
    atuple[6]=COSTMC;
    atuple[7]=Q2;


    // print some event data
    if(j<5) {
      cout << "event:" << j << endl;
      cout << "MC-eepi:" << MC1Energy << " " << MC2Energy << " " << MCEnergy << endl;
      cout << "theta:" << MC1Theta << " " << MC2Theta << " " << MCTheta << endl;
      cout << "phi:" << MC1Phi << " " << MC2Phi << " " << MCPhi << endl;
      cout << "COSTMC:" << COSTMC << " Q2:" << Q2 << endl;
    }

    // Analysis starts here
    // loop over charged candidate tracks
    Float_t cc_pp = 0, cc_E = 0, cc_mass = 0, cc_TH = 0;
    Float_t cc_px = 0, cc_py = 0, cc_pz = 0;

    TLorentzVector reTrack[4], QQ_RE;

    Int_t nelec_pair = 0;
    Double_t bestTH=-999;
    Double_t bestPH=-999;
    Double_t bestCOST=-999;
    Int_t ix=-1, iy=-1;

    for (Int_t nc1 = 0; nc1 < ncCand; nc1++) {
      PndPidCandidate *pc1 = (PndPidCandidate*)cCand_array->At(nc1);
      Int_t Charge1 = pc1->GetCharge();
      if (Charge1<0) continue;
      cc_pp = pc1->GetMomentum().Mag();
      //         if (cc_pp>eSystem) continue;
      cc_px = pc1->GetMomentum().Px();
      cc_py = pc1->GetMomentum().Py();
      cc_pz = pc1->GetMomentum().Pz();
      cc_E=TMath::Sqrt(cc_pp*cc_pp+mElec2);
      reTrack[0].SetPxPyPzE(cc_px,cc_py,cc_pz,cc_E);
      for (Int_t nc2 = 0; nc2 < ncCand; nc2++) {
	PndPidCandidate *pc2 = (PndPidCandidate*)cCand_array->At(nc2);
	Int_t Charge2 = pc2->GetCharge();
	if (Charge2>0) continue;
	cc_pp = pc2->GetMomentum().Mag();
	//           if (cc_pp>eSystem) continue;
	cc_px = pc2->GetMomentum().Px();
	cc_py = pc2->GetMomentum().Py();
	cc_pz = pc2->GetMomentum().Pz();
	cc_E=TMath::Sqrt(cc_pp*cc_pp+mElec2);
	reTrack[1].SetPxPyPzE(cc_px,cc_py,cc_pz,cc_E);
	nelec_pair++;
	// selection on best kinematics
	// elec1
	RE1Energy = reTrack[0].E();
	RE1Theta  = radeg*(reTrack[0].Theta());
	RE1Phi    = radeg*(reTrack[0].Phi());
	// elec2
	RE2Energy = reTrack[1].E();
	RE2Theta  = radeg*(reTrack[1].Theta());
	RE2Phi    = radeg*(reTrack[1].Phi());
	// quadrivecteur
	double elecRE=reTrack[0].Angle(reTrack[1].Vect());
	//       double Q2=4*mcTrack[0].E()*mcTrack[1].E()*(sin(elecMC/2))*(sin(elecMC/2));
	QQ_RE=reTrack[0]+reTrack[1];
	double Q2=QQ_RE.M2();
	// Center of mass angle of electrons
	// boost vector
	TLorentzVector *elecRE0 = new TLorentzVector(reTrack[0].Px(),
						     reTrack[0].Py(),
						     reTrack[0].Pz(),
						     reTrack[0].E());
	TLorentzVector *elecRE1 = new TLorentzVector(reTrack[1].Px(),
						     reTrack[1].Py(),
						     reTrack[1].Pz(),
						     reTrack[1].E());
	elecRE0->Boost(bSigma);
	elecRE1->Boost(bSigma);
	Double_t THRE0=radeg*(elecRE0->Theta());
	Double_t COSTRE=TMath::Cos(elecRE0->Theta());
	Double_t THRE1=radeg*(elecRE1->Theta());
	Double_t PHRE0=radeg*(elecRE0->Phi());
	Double_t PHRE1=radeg*(elecRE1->Phi());
	Double_t THtot=THRE0+THRE1;
	Double_t PHtot=PHRE0-PHRE1;
	if(PHtot< -90) PHtot+360;
	if(abs(THtot-180) < abs(bestTH-180)) {
	  bestTH=THtot;
	  bestPH=PHtot;
	  bestCOST=COSTRE;
	  ix=nc1;
	  iy=nc2;
	}
	// print some event data
	if(j<5) {
	  cout << "event:" << j << endl;
	  cout << "RE-eepi:" << RE1Energy << " " << RE2Energy << endl;
	  cout << "theta:" << RE1Theta << " " << RE2Theta  << endl;
	  cout << "phi:" << RE1Phi << " " << RE2Phi << endl;
	  cout << "THRE0:" << THRE0 << "THRE1:" << THRE1 << endl;
	  cout << "THtot:" << THtot << " Q2:" << Q2 << endl;
	  cout << "PHtot:" << PHtot << endl;
	}
      }       // nc2
    }       // nc1
    h_costheta-> Fill(COSTRE);
    h_dtheta->Fill(bestTH);
    h_dphi->Fill(bestPH);
    hpair->Fill((double) nelec_pair);
    // Rebelotte with the best solution
    if(nelec_pair<1) continue;
    //
    PndPidCandidate *pc1 = (PndPidCandidate*)cCand_array->At(ix);
    PndPidCandidate *pc2 = (PndPidCandidate*)cCand_array->At(iy);
    // positron
    cc_pp = pc1->GetMomentum().Mag();
    cc_px = pc1->GetMomentum().Px();
    cc_py = pc1->GetMomentum().Py();
    cc_pz = pc1->GetMomentum().Pz();
    cc_E=TMath::Sqrt(cc_pp*cc_pp+mElec2);
    reTrack[0].SetPxPyPzE(cc_px,cc_py,cc_pz,cc_E);
    // electron
    cc_pp = pc2->GetMomentum().Mag();
    cc_px = pc2->GetMomentum().Px();
    cc_py = pc2->GetMomentum().Py();
    cc_pz = pc2->GetMomentum().Pz();
    cc_E=TMath::Sqrt(cc_pp*cc_pp+mElec2);
    reTrack[1].SetPxPyPzE(cc_px,cc_py,cc_pz,cc_E);

    // elec1
    // detector data
    PndPidProbability *drc_ele = (PndPidProbability *)drc_array->At(ix);
    PndPidProbability *disc_ele = (PndPidProbability *)disc_array->At(ix);
    PndPidProbability *mvd_ele = (PndPidProbability *)mvd_array->At(ix);
    PndPidProbability *stt_ele = (PndPidProbability *)stt_array->At(ix);
    PndPidProbability *emcb_ele = (PndPidProbability *)emcb_array->At(ix);
    Double_t k_drc_e = drc_ele->GetElectronPidProb();
    Double_t k_disc_e = disc_ele->GetElectronPidProb();
    Double_t k_mvd_e = mvd_ele->GetElectronPidProb();
    Double_t k_stt_e = stt_ele->GetElectronPidProb();
    Double_t k_emcb_e = emcb_ele->GetElectronPidProb();
    Double_t xx_e = (k_drc_e/(1-k_drc_e))*(k_disc_e/(1-k_disc_e))
      *(k_mvd_e/(1-k_mvd_e))*(k_stt_e/(1-k_stt_e))
      *(k_emcb_e/(1-k_emcb_e));
    Double_t k_comb_e = xx_e/(xx_e+1);
    atuple[ 8]=reTrack[0].P();
    atuple[ 9]=reTrack[0].Theta();
    atuple[10]=reTrack[0].Phi();
    atuple[11]=pc1->GetEmcRawEnergy();
    atuple[12]=pc1->GetEmcNumberOfCrystals();
    atuple[13]=k_emcb_e;
    atuple[14]=k_stt_e;
    atuple[15]=k_disc_e;
    atuple[16]=k_drc_e;
    atuple[17]=k_mvd_e;
    atuple[18]=k_comb_e;

    // elec2
    PndPidProbability *drc_posi = (PndPidProbability *)drc_array->At(iy);
    PndPidProbability *disc_posi = (PndPidProbability *)disc_array->At(iy);
    PndPidProbability *mvd_posi = (PndPidProbability *)mvd_array->At(iy);
    PndPidProbability *stt_posi = (PndPidProbability *)stt_array->At(iy);
    PndPidProbability *emcb_posi = (PndPidProbability *)emcb_array->At(iy);
    Double_t k_drc_p = drc_posi->GetElectronPidProb();
    Double_t k_disc_p = disc_posi->GetElectronPidProb();
    Double_t k_mvd_p = mvd_posi->GetElectronPidProb();
    Double_t k_stt_p = stt_posi->GetElectronPidProb();
    Double_t k_emcb_p = emcb_posi->GetElectronPidProb();
    Double_t xx_p = (k_drc_p/(1-k_drc_p))*(k_disc_p/(1-k_disc_p))
      *(k_mvd_p/(1-k_mvd_p))*(k_stt_p/(1-k_stt_p))
      *(k_emcb_p/(1-k_emcb_p));
    Double_t k_comb_p = xx_p/(xx_p+1);

    atuple[19]=reTrack[1].P();
    atuple[20]=reTrack[1].Theta();
    atuple[21]=reTrack[1].Phi();
    atuple[22]=pc2->GetEmcRawEnergy();
    atuple[23]=pc2->GetEmcNumberOfCrystals();
    atuple[24]=k_emcb_p;
    atuple[25]=k_stt_p;
    atuple[26]=k_disc_p;
    atuple[27]=k_drc_p;
    atuple[28]=k_mvd_p;
    atuple[29]=k_comb_p;
    atuple[30]=bestTH;
    atuple[31]=bestPH;
    atuple[32]=bestCOST;

    NTev->Fill(atuple);
    if(j<5) {
      cout << "k_comb_e:" << k_comb_e<< " k_comb_p:" << k_comb_p << endl;
    }
    // standard cuts for histos
    Bool_t com_ele_1 = drc_ele->GetElectronPidProb() > 0.05 &&
      disc_ele->GetElectronPidProb() > 0.05 &&
      mvd_ele->GetElectronPidProb() > 0.05 &&
      stt_ele->GetElectronPidProb() > 0.05 &&
      emcb_ele->GetElectronPidProb() > 0.05;
    Bool_t com_ele_2 = k_comb_e > 0.9;
    Bool_t com_ele_3 = pc1->GetEmcNumberOfCrystals() > 5;
    Bool_t com_ele_4 = bestTH >=178. && bestTH <= 182.;
    Bool_t com_ele_5 = bestPH >=178. && bestPH <= 182.;
    Bool_t com_posi_1 = drc_posi->GetElectronPidProb() > 0.05 &&
      disc_posi->GetElectronPidProb() > 0.05 &&
      mvd_posi->GetElectronPidProb() > 0.05 &&
      stt_posi->GetElectronPidProb() > 0.05 &&
      emcb_posi->GetElectronPidProb() > 0.05;
    Bool_t com_posi_2 =  k_comb_p > 0.9;
    Bool_t com_posi_3 = pc2->GetEmcNumberOfCrystals() > 5;
    if(j<5) {
      cout << com_ele_1 << com_ele_2 << com_ele_3 << com_ele_4
	   << com_ele_4 << com_posi_1 << com_posi_2 << com_posi_3 << endl;
    }
    if (com_ele_1&&com_ele_2&&com_ele_3&&com_ele_4
        &&com_ele_5&&com_posi_1&&com_posi_2&&com_posi_3) {
      h_costheta_sel->Fill(bestCOST);

    }
  }  // loop j over events

  cout << " NEVcount: " << NEVcount << endl;
  h_efficiency->Divide(h_costheta_sel,h_costheta_mc,1,1,"B");

  // output

  /*
    cMCelec1->cd(1); hcutE->Draw();
    cMCelec1->cd(2); hcutx->Draw();
    cMCelec1->cd(3); hcuty->Draw();
    cMCelec1->cd(4); hcutz->Draw();
    //   cMCelec1->cd(5); hpt2->Draw();
    cMCelec1->cd(0);
  */

  cMA->cd(1);gPad->SetLogy(); h_costheta_mc->Draw();
  cMA->cd(2);gPad->SetLogy(); h_costheta->Draw();
  cMA->cd(3);gPad->SetLogy(); h_costheta_sel->Draw();
  cMA->cd(4);gPad->SetLogy(); h_efficiency->Draw();
  cMA->cd(5);gPad->SetLogy(); h_dtheta->Draw();
  cMA->cd(6);gPad->SetLogy(); h_dphi->Draw();
  cMA->cd(7); hpair->Draw();
  cMA->cd(0);

  out->cd();

  hpair->Write();
  h_costheta_mc->Write();
  h_costheta->Write();
  h_costheta_sel->Write();
  h_efficiency->Write();
  h_dtheta->Write();
  h_dphi->Write();


  NTev->Write();

  out->Save();




  cout << " Yahoo! " << endl;

}

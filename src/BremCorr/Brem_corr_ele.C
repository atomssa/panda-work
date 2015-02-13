void Brem_corr_ele() {

  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();

  TString Directory[1]={"/vol0/oct_2013/fichierroot/electron10000/barrel/1GeV/"};
  //TString Directory[1]={"/vol0/simulations_stt/macros_new/fichierroot_new/e100000/m2tfpa/"};

  Int_t Did=0;

  TString name = "_complete.root";
  TString inRecoFile   = (Directory[Did])+"reco"+name;
  TString inDigiFile   = Directory[Did]+"digi"+name;
  TString inSimFile    = Directory[Did]+"sim"+name;
  TString inPidSTTFile = Directory[Did]+"pid"+name;
  TString outAnaFile   = Directory[Did]+"brem_corr_ele_new"+name;

  TStopwatch timer;
  timer.Start();

  TFile *inFile = TFile::Open(inPidSTTFile,"READ");
  TTree *lhe=(TTree *) inFile->Get("cbmsim") ;
  TFile *outFile = TFile::Open(outAnaFile,"new");

  // adding other files as friends
  lhe->AddFriend("cbmsim",inSimFile);
  lhe->AddFriend("cbmsim",inDigiFile);
  lhe->AddFriend("cbmsim",inRecoFile);

  PndEmcMapper::Init(1);

  TClonesArray* mc_array=new TClonesArray("PndMCTrack");
  lhe->SetBranchAddress("MCTrack", &mc_array);

  TClonesArray* ChargedCand_array=new TClonesArray("PndPidCandidate");
  lhe->SetBranchAddress("PidChargedCand", &ChargedCand_array);

  TClonesArray* NeutralCand_array=new TClonesArray("PndPidCandidate");
  lhe->SetBranchAddress("PidNeutralCand", &NeutralCand_array);

  TClonesArray* emc_array=new TClonesArray("PndEmcBump");
  lhe->SetBranchAddress("EmcBump", &emc_array);

  TClonesArray* emc_array1=new TClonesArray("PndEmcCluster");
  lhe->SetBranchAddress("EmcCluster", &emc_array1);

  TClonesArray* emc_digi_array=new TClonesArray("PndEmcDigi");
  lhe->SetBranchAddress("EmcDigi", &emc_digi_array);

  TH1F *h_phi_bump = new TH1F("h_phi_bump","phi bumps",160,-180,180);

  TH1F *h_resol_corr_sep = new TH1F("h_resol_corr_sep","Momentum resolution with correction(method1+2)",500,-0.2,1);
  TH1F *h_resol_corr_all = new TH1F("h_resol_corr_all","Momentum resolution with correction(method1)",500,-0.2,1);

  TH1F *h_resol = new TH1F("h_resol","Momentum resolution w/o correction",500,-0.2,1);
  TH1F *h_photon_energy = new TH1F("h_photon_energy","E_ph_MC - E_ph_found",400,-2,2);
  TH1F *h_case = new TH1F("h_case","The photon case  ",5,0,5);
  TH1F *h_Eph_case1 = new TH1F("h_Eph_case1","First case photons energy",10000,0,10);
  TH1F *h_Eph_case2 = new TH1F("h_Eph_case2","Second case photons energy",10000,0,10);

  Int_t NTevents=lhe->GetEntriesFast();
  cout << "Number of events: " << NTevents << endl;

  for (Int_t j= 0; j<NTevents; j++){
    lhe->GetEntry(j);

    Float_t McMomOfEle = 0, McThetaOfEle = 0, McPhiOfEle = 0;
    Float_t RecMomOfEle = 0, RecThetaOfEle = 0, RecPhiOfEle = 0;
    Float_t PhotonTotEnergySep = 0., PhotonTotEnergyMerg = 0.;
    Float_t RecMomOfEleCor = 0, RecThetaOfEleCor = 0, RecPhiOfEleCor = 0;
    Float_t PhotonTotMCEnergy = 0, PhotonTotEnergy = 0;

    Int_t IndiceEleTrack = 0;

    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< endl;
    cout << "Event No."<< j << endl;

    PndMCTrack *EleMcTrack = (PndMCTrack *) mc_array->At(0);
    McMomOfEle = EleMcTrack->GetMomentum().Mag();
    McThetaOfEle = EleMcTrack->GetMomentum().Theta()*TMath::RadToDeg();
    McPhiOfEle = EleMcTrack->GetMomentum().Phi()*TMath::RadToDeg();

    Int_t NumOfMCTrack = mc_array->GetEntriesFast();

    for (Int_t k = 0; k < NumOfMCTrack; k++)
      {
	PndMCTrack *McTrack = (PndMCTrack *) mc_array->At(k);
	Int_t pdgcode = McTrack->GetPdgCode();
	Int_t MotherID = McTrack->GetMotherID();
	if (MotherID == 0 && pdgcode == 22)
	  {
            Double_t Sx = McTrack->GetStartVertex().x();
            Double_t Sy = McTrack->GetStartVertex().y();
            Double_t Sz = McTrack->GetStartVertex().z();
            Double_t Sr = sqrt(pow(Sx,2)+pow(Sy,2));
            if (Sr < 42. && Sz<195. )
	      {
                Double_t PhotonMCEnergy = McTrack->Get4Momentum().E();
                PhotonTotMCEnergy += PhotonMCEnergy;
	      }
	  }
      }
    //cout << McMomOfEle << endl;
    Int_t NumOfChargedCand = ChargedCand_array->GetEntriesFast();
    if(NumOfChargedCand == 0)
      {
	cout << "no primary electron track found!!" << endl;
	continue;
      }

    // find closest reconstructed track to Mc by momentum
    for (Int_t i_ChargeCand=0; i_ChargeCand<NumOfChargedCand; i_ChargeCand++)
      {
	PndPidCandidate *EleRecTrack = (PndPidCandidate*) ChargedCand_array->At(i_ChargeCand);
	if (EleRecTrack->GetMcIndex()!=0) continue;
	PndMCTrack* CorrespondMcTrack = (PndMCTrack*)mc_array->At(EleRecTrack->GetMcIndex());
	if (CorrespondMcTrack->GetMotherID()!=-1) continue;
	if (  fabs(RecMomOfEle-McMomOfEle) > fabs(EleRecTrack->GetMomentum().Mag()-McMomOfEle) )
	  {
	    RecMomOfEle = EleRecTrack->GetMomentum().Mag();
	    RecThetaOfEle = EleRecTrack->GetMomentum().Theta()*TMath::RadToDeg();
	    RecPhiOfEle = EleRecTrack->GetMomentum().Phi()*TMath::RadToDeg();
	    IndiceEleTrack = i_ChargeCand;
	  }
      } // chargedcandidate loop


    RecMomOfEleCor = RecMomOfEle;

    h_resol->Fill((McMomOfEle-RecMomOfEle)/McMomOfEle);

    Int_t NumOfNeutralCand = NeutralCand_array->GetEntriesFast();

    for(Int_t i_NeutralCand = 0; i_NeutralCand<NumOfNeutralCand; i_NeutralCand++)
      {
	Float_t PhotonEnergySep = 0, PhotonThetaSep = 0, PhotonPhiSep = 0;
	PndPidCandidate *PhotonCand = (PndPidCandidate *) NeutralCand_array->At(i_NeutralCand);
	PhotonEnergySep = PhotonCand->GetEmcCalEnergy();
	PndEmcBump *PhotonBump = (PndEmcBump *) emc_array->At(PhotonCand->GetEmcIndex());
	PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
	PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();


	Float_t Pt = RecMomOfEle*TMath::Sin(RecThetaOfEle/TMath::RadToDeg());
	Float_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
	Float_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(RecThetaOfEle/57.3)*57.3;
	Float_t RealDeltaPhi = PhotonPhiSep - RecPhiOfEle;
	Float_t RealDeltaTheta = PhotonThetaSep - RecThetaOfEle;
	if (RecThetaOfEle <= 23.)
	  {
	    Float_t PhiCutUp = DeltaPhiForward, PhiCutDown = -1.;
	    Float_t ThetaCutUp = 2., ThetaCutDown = -2.;
	  }
	else if (RecThetaOfEle > 23.)
	  {
	    Float_t PhiCutUp = DeltaPhiBarrel, PhiCutDown = -1.;
	    Float_t ThetaCutUp = 2., ThetaCutDown = -2.;
	  }
	Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
	Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
	if (PhiCut && ThetaCut) PhotonTotEnergySep += PhotonEnergySep;

      }//loop neutralcand

    if(PhotonTotEnergySep > RecMomOfEle/100.)
      {
	RecMomOfEleCor += PhotonTotEnergySep;
	cout << "Seperated photon found: " << PhotonTotEnergySep <<" GeV" << endl;
	h_case->Fill(1);
	h_Eph_case1->Fill(PhotonTotEnergySep);
      }

    h_resol_corr_sep->Fill((McMomOfEle-RecMomOfEleCor)/McMomOfEle);

    PndPidCandidate *SelectedEleCand = (PndPidCandidate*) ChargedCand_array->At(IndiceEleTrack);
    if (SelectedEleCand->GetEmcIndex() < 0) continue;
    PndEmcBump *EleBump = (PndEmcBump *) emc_array->At(SelectedEleCand->GetEmcIndex());
    Int_t EleRefCluster = EleBump->GetClusterIndex();
    if (EleRefCluster < 0) continue;
    PndEmcCluster *EleCluster = (PndEmcCluster *) emc_array1->At(EleRefCluster);
    Int_t EleDigiSize = EleCluster->DigiList().size();
    std::vector<Int_t> EleDigiList = EleCluster->DigiList();

    for (Int_t i_digi = 0; i_digi<EleDigiSize; i_digi++)
      {

        PndEmcDigi *EleEmcDigi = (PndEmcDigi *) emc_digi_array->At(EleDigiList[i_digi]);
        Double_t EleEmcDigiPhi = EleEmcDigi->GetPhi()*TMath::RadToDeg();
        Double_t EleEmcDigiEnergy = EleEmcDigi->GetEnergy();
        h_phi_bump->Fill(EleEmcDigiPhi,EleEmcDigiEnergy);
      }

    Int_t TotNumOfHitPhi = h_phi_bump->GetNbinsX();

    Double_t DepoEnergyList[100];
    Int_t DepoENoOfBinList[100];
    Int_t GapSizeList[100];
    Int_t DepoEListIndex = 1;

    DepoEnergyList[0] = 0;
    DepoENoOfBinList[0] = -10;
    GapSizeList[0] = -1;
    for (Int_t i_phi = 1; i_phi <= TotNumOfHitPhi; i_phi++)
      {
	Double_t BinValue = h_phi_bump->GetBinContent(i_phi);
	if (BinValue != 0)
          {
	    DepoEnergyList[DepoEListIndex] = BinValue;
	    DepoENoOfBinList[DepoEListIndex] = i_phi;
	    GapSizeList[DepoEListIndex] = DepoENoOfBinList[DepoEListIndex] - DepoENoOfBinList[DepoEListIndex-1];
	    DepoEListIndex++;
          }
      }
    DepoEnergyList[DepoEListIndex] = 0;
    DepoENoOfBinList[DepoEListIndex] = -1;
    GapSizeList[DepoEListIndex] = -1;

    Int_t StartIndice = 0;
    for (Int_t i = 1;i < DepoEListIndex;i++)
      {
	if (GapSizeList[i] > GapSizeList[StartIndice]) StartIndice = i;
      }
    if(StartIndice > 1)
      {
	Int_t i_corr = 0;
	Double_t DepoEnergyListCorr[100];
	DepoEnergyListCorr[0] = DepoEnergyList[0];
	for (Int_t i2 = StartIndice; i2 < DepoEListIndex;i2++)
	  {
	    i_corr++;
	    DepoEnergyListCorr[i_corr] = DepoEnergyList[i2];
	  }
	for(Int_t i2 = 1;i2 < StartIndice;i2++)
	  {
	    i_corr++;
	    DepoEnergyListCorr[i_corr] = DepoEnergyList[i2];
	  }
	DepoEnergyListCorr[DepoEListIndex] = DepoEnergyList[DepoEListIndex];
      } //StartIndice
    else if (StartIndice = 1)
      {
	for (Int_t i2 = 0; i2 <= DepoEListIndex ; i2++) DepoEnergyListCorr[i2] = DepoEnergyList[i2];
      }

    Int_t Case[100];
    Double_t PhiBump[100], Poid[100];
    Int_t PhiBumpIndex = 0, PoidIndex = 0;
    Case[0] = -3;
    Poid[PoidIndex] = 0;
    Int_t IndiceVally = 0;
    for (Int_t n_sel = 1; n_sel < DepoEListIndex; n_sel++)
      {
	if(DepoEnergyListCorr[n_sel-1] < DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] < DepoEnergyListCorr[n_sel+1] ) Case[n_sel] = 1;
	else if(DepoEnergyListCorr[n_sel-1] < DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] > DepoEnergyListCorr[n_sel+1] )
          {
	    PoidIndex++;
	    Case[n_sel] = 0;
	    Poid[PoidIndex] = DepoEnergyListCorr[n_sel];
          }
	else if(DepoEnergyListCorr[n_sel-1] > DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] > DepoEnergyListCorr[n_sel+1] ) Case[n_sel] = -1;
	else if(DepoEnergyListCorr[n_sel-1] > DepoEnergyListCorr[n_sel] && DepoEnergyListCorr[n_sel] < DepoEnergyListCorr[n_sel+1] ) Case[n_sel] = -2;
      }
    Poid[PoidIndex+1] = 0;

    Int_t iPoid = 0;

    for (Int_t n_sel = 1; n_sel < DepoEListIndex; n_sel++)
      {
	if (Case[n_sel] == -2 || n_sel == DepoEListIndex-1)
          {
	    iPoid++;
	    PhiBump[PhiBumpIndex] = DepoEnergyListCorr[IndiceVally]*(Poid[iPoid]/(Poid[iPoid]+Poid[iPoid-1]));
	    for(Int_t p = IndiceVally;p < n_sel;p++) PhiBump[PhiBumpIndex] += DepoEnergyListCorr[p];
	    PhiBump[PhiBumpIndex] += DepoEnergyListCorr[n_sel]*(Poid[iPoid]/(Poid[iPoid]+Poid[iPoid+1]));
	    PhiBumpIndex++;
	    IndiceVally = n_sel;
          }
      }
    PhiBumpIndex-=1;

    Double_t EnergyCut = 0.15/TMath::Sin(RecThetaOfEle*TMath::DegToRad());
    Double_t EleEnergy = 0;
    for (Int_t i_phibump = PhiBumpIndex; i_phibump >= 0; i_phibump--)
      {

	if(PhiBump[i_phibump] > EnergyCut)
	  {
	    for(Int_t r = 0; r<i_phibump; r++) PhotonTotEnergyMerg += PhiBump[r];
	    for(r = i_phibump; r <= PhiBumpIndex;r++) EleEnergy += PhiBump[r];
	    i_phibump = -1;
	  }
      }

    h_phi_bump->Reset();

    if(PhotonTotEnergyMerg > RecMomOfEle/100.)
      {
	RecMomOfEleCor += PhotonTotEnergyMerg;
	cout << "Merged photon found: " << PhotonTotEnergyMerg <<" GeV" << endl;
	h_case->Fill(2);
	h_Eph_case2->Fill(PhotonTotEnergyMerg);

      }

    h_resol_corr_all->Fill((McMomOfEle-RecMomOfEleCor)/McMomOfEle);
    PhotonTotEnergy = PhotonTotEnergyMerg + PhotonTotEnergySep;
    cout << PhotonTotMCEnergy << "  " << PhotonTotEnergy  << endl;
    if (PhotonTotMCEnergy != 0 || PhotonTotEnergy != 0) h_photon_energy->Fill(PhotonTotMCEnergy-PhotonTotEnergy);
    if(PhotonTotEnergyMerg < RecMomOfEle/100. && PhotonTotEnergySep < RecMomOfEle/100.) h_case->Fill(0);
  }//loop NTevents

  outFile->cd();
  h_resol->Write();
  h_resol_corr_sep->Write();
  h_resol_corr_all->Write();
  h_photon_energy->Write();
  h_case->Write();
  h_Eph_case1->Write();
  h_Eph_case2->Write();
  outFile->Save();

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

  cout << " Test passed" << endl;
  cout << " All ok " << endl;

}

#include "TROOT.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "TStopwatch.h"
#include <iostream>

//class PndPidCandidate;
//class PndMCTrack;
//class PndEmcDigi;
//class PndEmcBump;
//class PndEmcCluster;

//#include "PndPidCandidate.h"
//#include "PndMCTrack.h"
//#include "PndEmcDigi.h"
//#include "PndEmcBump.h"
//#include "PndEmcCluster.h"

using namespace std;

void Brem_corr_ele(TString dir="./")
{

  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  //gSystem->Load("libPndData");
  //gSystem->Load("libPid");

  TStopwatch timer;
  timer.Start();

  //TChain *lhe = new TChain("cbmsim","cbmsim");
  //for (int i =0; i<10; ++i) {
  //  if (i==34||i==44) continue;
  //  TString fName= Form("/vol0/panda/work/jpsi_pi0/grid.out/esim_trunk/runall.%d.out/pid_complete.root",i);
  //  cout << "Adding " << fName << " ...  " << endl;
  //  lhe->AddFile(fName);
  //}
  //TFile *outFile = TFile::Open("brem_ala_binsong.root","RECREATE");

  TString inFileName = dir+"/pid_complete.root";
  TString outFileName   = dir+"/brem_ala_binsong.root";
  cout << "inFile= " << inFileName << endl;
  cout << "outFile= " << outFileName << endl;
  TFile *inFile = TFile::Open(inFileName);
  TTree *lhe=(TTree *) inFile->Get("cbmsim") ;
  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  if (1) {
    lhe->AddFriend("cbmsim",dir+"/sim_complete.root");
    lhe->AddFriend("cbmsim",dir+"/digi_complete.root");
    lhe->AddFriend("cbmsim",dir+"/reco_complete.root");
  }

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
  TH1F *h_case = new TH1F("h_case","The photon case  ",5,0,5);
  TH1F *h_Eph_case1 = new TH1F("h_Eph_case1","First case photons energy",10000,0,10);
  TH1F *h_Eph_case2 = new TH1F("h_Eph_case2","Second case photons energy",10000,0,10);

  Int_t NTevents=lhe->GetEntries();
  cout << "Number of events: " << NTevents << endl;
  int verb = 0;

  for (Int_t j= 0; j<NTevents; j++){

    lhe->GetEntry(j);

    if (verb||j%100==0) cout << "<<<< Event No. " << j << " <<<<<" << endl;
    int NumOfChargedCand = ChargedCand_array->GetEntriesFast();
    if (verb) cout << "NumChCand= " << NumOfChargedCand << endl;
    for (Int_t i_ChargeCand=0; i_ChargeCand<NumOfChargedCand; i_ChargeCand++) {

      Double_t PhotonTotEnergySep = 0., PhotonTotEnergyMerg = 0.;
      //Double_t RecMomOfEleCor = 0, RecThetaOfEleCor = 0, RecPhiOfEleCor = 0;

      PndPidCandidate *EleRecTrack = (PndPidCandidate*) ChargedCand_array->At(i_ChargeCand);
      if (verb) cout << "EleRecTrackMcIndex= " << EleRecTrack->GetMcIndex() << endl;
      if (EleRecTrack->GetMcIndex()>4||EleRecTrack->GetMcIndex()<0) continue; // five primary tracks per event

      PndMCTrack* CorrespondMcTrack = (PndMCTrack*)mc_array->At(EleRecTrack->GetMcIndex());
      //if (CorrespondMcTrack->GetMotherID()!=-1) continue;
      Double_t McMomOfEle = CorrespondMcTrack->GetMomentum().Mag();
      Double_t McThetaOfEle = CorrespondMcTrack->GetMomentum().Theta()*TMath::RadToDeg();
      Double_t McPhiOfEle = CorrespondMcTrack->GetMomentum().Phi()*TMath::RadToDeg();

      Double_t RecMomOfEle = EleRecTrack->GetMomentum().Mag();
      Double_t RecThetaOfEle = EleRecTrack->GetMomentum().Theta()*TMath::RadToDeg();
      Double_t RecPhiOfEle = EleRecTrack->GetMomentum().Phi()*TMath::RadToDeg();
      Double_t RecMomOfEleCor = RecMomOfEle;

      Double_t mom_min = 0.5;
      Double_t mom_max = 1.0;
      Double_t th_min = 0;
      Double_t th_max = 180;
      //if ((McMomOfEle<mom_min)||(McMomOfEle>mom_max)) continue;
      //if ((McThetaOfEle<th_min)||(McThetaOfEle>th_max)) continue;
      //if (verb) {
      //	cout << "========== Track " << i_ChargeCand << " ===============" << endl;
      //	cout << "McMom= " << McMomOfEle << " McThe= " << McThetaOfEle << endl;
      //}

      //if (  fabs(RecMomOfEle-McMomOfEle) > fabs(EleRecTrack->GetMomentum().Mag()-McMomOfEle) ) { /*??*/ }
      h_resol->Fill((McMomOfEle-RecMomOfEle)/McMomOfEle);

      Int_t NumOfNeutralCand = NeutralCand_array->GetEntriesFast();

      for(Int_t i_NeutralCand = 0; i_NeutralCand<NumOfNeutralCand; i_NeutralCand++)
	{
	  Double_t PhotonEnergySep = 0, PhotonThetaSep = 0, PhotonPhiSep = 0;
	  PndPidCandidate *PhotonCand = (PndPidCandidate *) NeutralCand_array->At(i_NeutralCand);
	  PhotonEnergySep = PhotonCand->GetEmcCalEnergy();
	  if (PhotonCand->GetEmcIndex()<0 || PhotonCand->GetEmcIndex()>emc_array->GetEntriesFast()) continue;
	  PndEmcBump *PhotonBump = (PndEmcBump *) emc_array->At(PhotonCand->GetEmcIndex());
	  PhotonThetaSep = PhotonBump->position().Theta()*TMath::RadToDeg();
	  PhotonPhiSep = PhotonBump->position().Phi()*TMath::RadToDeg();

	  Double_t Pt = RecMomOfEle*TMath::Sin(RecThetaOfEle/TMath::RadToDeg());
	  Double_t DeltaPhiBarrel = TMath::ASin(0.12/Pt)*2.*TMath::RadToDeg();
	  Double_t DeltaPhiForward = (0.6*2.0/Pt)*TMath::Tan(RecThetaOfEle/57.3)*57.3;
	  Double_t RealDeltaPhi = PhotonPhiSep - RecPhiOfEle;
	  Double_t RealDeltaTheta = PhotonThetaSep - RecThetaOfEle;

	  Double_t PhiCutUp=0, PhiCutDown=0;
	  Double_t ThetaCutUp=0, ThetaCutDown=0;
	  if (RecThetaOfEle <= 23.)
	    {
	      PhiCutUp = DeltaPhiForward; PhiCutDown = -1.;
	      ThetaCutUp = 2.; ThetaCutDown = -2.;
	    }
	  else if (RecThetaOfEle > 23.)
	    {
	      PhiCutUp = DeltaPhiBarrel; PhiCutDown = -1.;
	      ThetaCutUp = 2.; ThetaCutDown = -2.;
	    }
	  Bool_t PhiCut = RealDeltaPhi <= PhiCutUp && RealDeltaPhi >= PhiCutDown;
	  Bool_t ThetaCut = RealDeltaTheta <= ThetaCutUp && RealDeltaTheta >= ThetaCutDown;
	  if (PhiCut && ThetaCut) PhotonTotEnergySep += PhotonEnergySep;

	}//loop neutralcand

      if(PhotonTotEnergySep > RecMomOfEle/100.)
	{
	  RecMomOfEleCor += PhotonTotEnergySep;
	  if (verb) cout << "Seperated photon found: " << PhotonTotEnergySep <<" GeV" << endl;
	  h_case->Fill(1);
	  h_Eph_case1->Fill(PhotonTotEnergySep);
	}

      h_resol_corr_sep->Fill((McMomOfEle-RecMomOfEleCor)/McMomOfEle);

      if (EleRecTrack->GetEmcIndex()<0 || EleRecTrack->GetEmcIndex()> emc_array->GetEntriesFast()) continue;
      PndEmcBump *EleBump = (PndEmcBump *) emc_array->At(EleRecTrack->GetEmcIndex());
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
      Double_t DepoEnergyListCorr[100];
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
	      for(Int_t r = i_phibump; r <= PhiBumpIndex;r++) EleEnergy += PhiBump[r];
	      i_phibump = -1;
	    }
	}

      h_phi_bump->Reset();

      if(PhotonTotEnergyMerg > RecMomOfEle/100.)
	{
	  RecMomOfEleCor += PhotonTotEnergyMerg;
	  if (verb) cout << "Merged photon found: " << PhotonTotEnergyMerg <<" GeV" << endl;
	  h_case->Fill(2);
	  h_Eph_case2->Fill(PhotonTotEnergyMerg);
	}

      Double_t PhotonTotEnergy = 0;
      h_resol_corr_all->Fill((McMomOfEle-RecMomOfEleCor)/McMomOfEle);
      PhotonTotEnergy = PhotonTotEnergyMerg + PhotonTotEnergySep;
      if (verb) cout << "PhotonTotEnergy= " << PhotonTotEnergy  << endl;
      if(PhotonTotEnergyMerg < RecMomOfEle/100. && PhotonTotEnergySep < RecMomOfEle/100.) h_case->Fill(0);

    }

  }//loop NTevents

  outFile->cd();
  h_resol->Write();
  h_resol_corr_sep->Write();
  h_resol_corr_all->Write();
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


enum DetectorId {
/** kDCH must be the 1st id, and kHYP must be the last one. Please put new detectors in between!! **/
    kDCH,kDRC,kDSK,kEMC,kGEM,kLUMI,kMDT,kMVD,kRPC,kSTT,kTPC,kTOF,kFTS,kHYPG,kHYP};

void track_check(Int_t nEntries = 0)
{
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  Text_t  tOutputFile[50];
  
  TString inPidFile = "evt_pid.root", inSimFile = "evt_points.root";
 
  TFile *inFile = TFile::Open(inSimFile,"READ");
  
  TTree *tree=(TTree *) inFile->Get("cbmsim") ;
  tree->AddFriend("cbmsim",inPidFile);
    
  TClonesArray* cand_array=new TClonesArray("PndPidCandidate");
  tree->SetBranchAddress("PidChargedCand", &cand_array);

  TClonesArray* mc_array=new TClonesArray("PndMCTrack");
  tree->SetBranchAddress("MCTrack", &mc_array);

  TFile *out = TFile::Open("output.root","RECREATE");
  TNtuple *nt = new TNtuple("nt","nt","evt:id:mc_p:mc_theta:mc_phi:mc_pid:mc_stt:mc_tpc:mc_mvd:mc_gem:p:theta:phi:mult:stt:tpc:mvd");
  TNtuple *ntEvt = new TNtuple("ntEvt","ntEvt","evt:acc_stt:acc_tpc:eff:effct:cand:mc_D1:mc_D2:reco_D1:reco_D2:mc_all:reco_all");
  Float_t mass_k = 0.493677; 
  Float_t mass_pi =0.13957018;

  if (nEntries==0) nEntries =  tree->GetEntriesFast();
  for (Int_t j=0; j< nEntries; j++){
    tree->GetEntry(j);
    //if (cand_array->GetEntriesFast()==0) continue;
    if ((j%100)==0)    cout << "processing event " << j << "\n";
    
    Float_t mc_mom = 0, mc_theta = 0, mc_phi = 0;
    Float_t rec_mom = -1, rec_theta = -1, rec_phi = 0;
    Int_t stt_mccount = 0, tpc_mccount = 0,  reco_stt = 0, reco_tpc = 0,reco_mvd = 0, reco_gem = 0, reco_count = 0, reco_ctcount = 0;
    TLorentzVector mc_k[6], reco_k[6], mc_D1, reco_D1, mc_D2, reco_D2, mc_all, reco_all;
    for (Int_t mc = 0; mc < mc_array->GetEntriesFast(); mc++)
      {
	PndMCTrack *mctrack = (PndMCTrack*)mc_array->At(mc);
	if (mctrack->GetMotherID()!=-1) continue;
	if (mctrack->GetNPoints(kSTT))  stt_mccount++;
	if (mctrack->GetNPoints(kTPC))  tpc_mccount++;
	mc_k[mc] = mctrack->Get4Momentum();
	mc_mom = mctrack->GetMomentum().Mag();
	mc_theta = mctrack->GetMomentum().Theta()*TMath::RadToDeg();
	mc_phi = mctrack->GetMomentum().Phi()*TMath::RadToDeg();
	Int_t mc_pid = mctrack->GetPdgCode();
	Int_t cand_mult = 0;
	for (Int_t pp=0; pp<cand_array->GetEntriesFast(); pp++)
	  {
	    PndPidCandidate *pidCand = (PndPidCandidate*)cand_array->At(pp);
	    if (pidCand->GetMcIndex()!=mc) continue;
	    if ( (cand_mult==0) || ((cand_mult>0) && (fabs(rec_mom-mc_mom)> fabs(pidCand->GetMomentum().Mag()-mc_mom))) )
	      {
		rec_mom = pidCand->GetMomentum().Mag();
		rec_theta = pidCand->GetMomentum().Theta();
		rec_phi = pidCand->GetMomentum().Phi();
		if (mc==0||mc==3)
		  reco_k[mc].SetXYZM(pidCand->GetMomentum().X(),pidCand->GetMomentum().Y(),pidCand->GetMomentum().Z(),mass_k);
		else
		  reco_k[mc].SetXYZM(pidCand->GetMomentum().X(),pidCand->GetMomentum().Y(),pidCand->GetMomentum().Z(),mass_pi);
		reco_stt = pidCand->GetSttHits();
	       	reco_tpc = pidCand->GetTpcHits();	
		reco_mvd = pidCand->GetMvdHits();       
	      }	
	    cand_mult++;
	    
	  } // end of candidate loop
	
	if (cand_mult>0) reco_count++;
	if ((cand_mult>0) && ((reco_stt>0) || (reco_tpc>0)))  reco_ctcount++;
	Float_t ntuple_nt[] = {
	  j,mc, mc_mom,mc_theta,mc_phi, mc_pid,
	  mctrack->GetNPoints(kSTT), mctrack->GetNPoints(kTPC), mctrack->GetNPoints(kMVD), mctrack->GetNPoints(kGEM), 
	  rec_mom, rec_theta*TMath::RadToDeg(), rec_phi*TMath::RadToDeg(), cand_mult, reco_stt, reco_tpc, reco_mvd
	};
	nt->Fill(ntuple_nt);
	
      } // end of MC loop
    mc_D1 = mc_k[0] + mc_k[1] + mc_k[2];
    mc_D2 = mc_k[3] + mc_k[4] + mc_k[5];
    reco_D1 = reco_k[0] + reco_k[1] + reco_k[2];
    reco_D2 = reco_k[3] + reco_k[4] + reco_k[5];
    mc_all = mc_D1 + mc_D2;
    reco_all = reco_D1 + reco_D2;
    
    Float_t ntuple_evt[] = {
      j, stt_mccount, tpc_mccount, reco_count,  reco_ctcount,cand_array->GetEntriesFast(),
      mc_D1.M(), mc_D2.M(), reco_D1.M(), reco_D2.M(), mc_all.M(), reco_all.M()
    };
    ntEvt->Fill(ntuple_evt);
    
    
  } // end of event loop
  
  nt->Draw("mc_theta>>hMcTheta(75,0,150)");
  nt->Draw("mc_theta>>hMcSttTheta(75,0,150)","mc_stt>0");
  nt->Draw("mc_theta>>hMcSttMvdTheta(75,0,150)","(mc_stt>0)||(mc_mvd>0)"); 
  nt->Draw("mc_theta>>hMcTpcTheta(75,0,150)","mc_tpc>0");
  nt->Draw("mc_theta>>hMcTpcMvdTheta(75,0,150)","(mc_tpc>0)||(mc_tpc>0)");
  nt->Draw("mc_theta>>hMcMvdTheta(75,0,150)","mc_mvd>0");
  nt->Draw("mc_theta>>hMcGemTheta(75,0,150)","mc_gem>0");
  nt->Draw("mc_theta>>hRecoTheta(75,0,150)","mult>0"); 
  nt->Draw("mc_theta>>hRecoCtTheta(75,0,150)","mult>0&&(stt>0||tpc>0)");
  nt->Draw("mc_theta>>hRecoEffTheta(75,0,150)","mult>0"); 
  nt->Draw("mc_theta>>hRecoEffCtTheta(75,0,150)","mult>0&&(stt>0||tpc>0)");
  hRecoEffTheta->Divide(hMcTheta);
  hRecoEffCtTheta->Divide(hMcTheta);
  nt->Draw("mc_theta>>hMcAccSttTheta(75,0,150)","mc_stt>0");
  hMcAccSttTheta->Divide(hMcTheta);
  nt->Draw("mc_theta>>hMcAccTpcTheta(75,0,150)","mc_tpc>0");
  hMcAccTpcTheta->Divide(hMcTheta);
  
 
  nt->Draw("mc_phi>>hMcPhi(50,-200,200)");
  nt->Draw("mc_phi>>hMcSttPhi(50,-200,200)","mc_stt>0");
  nt->Draw("mc_phi>>hMcSttMvdPhi(50,-200,200)","(mc_stt>0)||(mc_mvd>0)"); 
  nt->Draw("mc_phi>>hMcTpcPhi(50,-200,200)","mc_tpc>0");
  nt->Draw("mc_phi>>hMcTpcMvdPhi(50,-200,200)","(mc_tpc>0)||(mc_mvd>0)");
  nt->Draw("mc_phi>>hMcMvdPhi(50,-200,200)","mc_mvd>0");
  nt->Draw("mc_phi>>hMcGemPhi(50,-200,200)","mc_gem>0");
  nt->Draw("mc_phi>>hRecoPhi(50,-200,200)","mult>0"); 
  nt->Draw("mc_phi>>hRecoCtPhi(50,-200,200)","mult>0&&(stt>0||tpc>0)");
  nt->Draw("mc_phi>>hRecoEffPhi(50,-200,200)","mult>0");
  nt->Draw("mc_phi>>hRecoEffCtPhi(50,-200,200)","mult>0&&(stt>0||tpc>0)");
  hRecoEffPhi->Divide(hMcPhi);
  hRecoEffCtPhi->Divide(hMcPhi);
  nt->Draw("mc_phi>>hMcAccSttPhi(50,-200,200)","mc_stt>0");
  hMcAccSttPhi->Divide(hMcPhi);
  nt->Draw("mc_phi>>hMcAccTpcPhi(50,-200,200)","mc_tpc>0");
  hMcAccTpcPhi->Divide(hMcPhi);

  nt->Draw("mc_p>>hMcP(100,0,5)");
  nt->Draw("mc_p>>hMcSttP(100,0,5)","mc_stt>0"); 
  nt->Draw("mc_p>>hMcSttMvdP(100,0,5)","(mc_stt>0)||(mc_mvd>0)");
  nt->Draw("mc_p>>hMcTpcP(100,0,5)","mc_tpc>0"); 
  nt->Draw("mc_p>>hMcTpcMvdP(100,0,5)","(mc_tpc>0)||(mc_mvd>0)");
  nt->Draw("mc_p>>hMcMvdP(100,0,5)","mc_mvd>0");
  nt->Draw("mc_p>>hMcGemP(100,0,5)","mc_gem>0");
  nt->Draw("mc_p>>hRecoP(100,0,5)","mult>0");
  nt->Draw("mc_p>>hRecoCtP(100,0,5)","mult>0&&(stt>0||tpc>0)");
  nt->Draw("mc_p>>hRecoEffP(100,0,5)","mult>0"); 
  nt->Draw("mc_p>>hRecoEffCtP(100,0,5)","mult>0&&(stt>0||tpc>0)");
  hRecoEffP->Divide(hMcP);
  hRecoEffCtP->Divide(hMcP);
  nt->Draw("mc_p>>hMcAccSttP(100,0,5)","mc_stt>0");
  hMcAccSttP->Divide(hMcP);
  nt->Draw("mc_p>>hMcAccTpcP(100,0,5)","mc_tpc>0");
  hMcAccTpcP->Divide(hMcP);

  nt->Draw("(mc_p-p):mc_p>>hresp_p(100,0,5,100,-0.5,0.5)","mult>0");
  nt->Draw("(mc_p-p):mc_theta>>hresp_theta(100,0,100,100,-0.5,0.5)","mult>0","colz");
  nt->Draw("(mc_theta-theta):mc_theta>>hrestheta_theta(100,0,100,100,-2,2)","mult>0");
  nt->Draw("(mc_theta-theta):mc_p>>hrestheta_p(100,0,5,100,-2,2)","mult>0");
  nt->Draw("(mc_phi-phi):mc_phi>>hresphi_phi(100,-200,200,100,-2,2)","mult>0");

  ntEvt->Draw("reco_D1>>hD(100,1.5,2.2)","eff==6");
  ntEvt->Draw("reco_D2>>hD2(100,1.5,2.2)","eff==6");
  hD->Add(hD2);
  ntEvt->Draw("reco_all>>hAll(100,3.2,4.4)","eff==6");
  out->cd();
  
  nt->Write();  
  ntEvt->Write();

  hMcTheta->Write(); hMcSttTheta->Write(); hMcSttMvdTheta->Write(); hMcMvdTheta->Write(); hMcGemTheta->Write(); hRecoTheta->Write();
  hMcTpcTheta->Write(); hMcTpcMvdTheta->Write(); hRecoEffTheta->Write(); hMcAccSttTheta->Write(); hMcAccTpcTheta->Write(); 
  hMcPhi->Write();   hMcSttPhi->Write();   hMcSttMvdPhi->Write();   hMcMvdPhi->Write();   hMcGemPhi->Write();   hRecoPhi->Write();
  hMcTpcPhi->Write();   hMcTpcMvdPhi->Write();  hRecoEffPhi->Write(); hMcAccSttPhi->Write(); hMcAccTpcPhi->Write(); 
  hMcP->Write();     hMcSttP->Write();     hMcSttMvdP->Write();     hMcMvdP->Write();     hMcGemP->Write();     hRecoP->Write();
  hMcTpcP->Write();     hMcTpcMvdP->Write(); hRecoEffP->Write();hMcAccSttP->Write(); hMcAccTpcP->Write(); 
  hRecoCtTheta->Write(); 
  hRecoCtPhi->Write(); 
  hRecoCtP->Write(); 
  hRecoEffCtTheta->Write(); 
  hRecoEffCtPhi->Write(); 
  hRecoEffCtP->Write(); 
  

  
  hresp_p->Write(); 
  hresp_theta->Write();
  hrestheta_theta->Write(); 
  hrestheta_p->Write();
  hresphi_phi->Write();
  hD->Write();
  hAll->Write();
  
  out->Save();
  
}	

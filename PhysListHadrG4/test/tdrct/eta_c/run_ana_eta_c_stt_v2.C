// Use PndAnalysis intead of PndEventReader
class TCandList;
class TCandidate;
class TFitParams;

void run_ana_eta_c_stt_v2(int nevts=0, bool usePID=true)
{
	TString simFile = "evt_points_stt.root";
	TString recoFile  = "evt_reco_stt.root";
	TString pidFile  = "evt_pid_stt.root";
	TString parFile = "evt_params_stt.root";
	TString outFile = "etac_histo_stt.root";

	enum DetectorId {kDCH,kDRC,kDSK,kEMC,kGEM,kLUMI,kMDT,
	kMVD,kRPC,kSTT,kTPC,kTOF,kFTS,kHYPG,kHYP};
	gStyle->SetOptFit(1011);

	TStopwatch timer;
	timer.Start();
	
	gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");

	//////////////// Histogram initialization ///////////////////////////////
	TH1F *h_etac_nocut=new TH1F("h_etac_nocut","#eta_{c}m(#phi,#phi), (no cuts);E, GeV",100,2.5,3.5);
	TH1F *h_mphi_nocuts=new TH1F("h_mphi_nocuts","#phi: m(K+ K-) (no cuts);E, GeV",100,0.95,1.5);

	TH1F *h_etac_pid=new TH1F("h_etac_pid","#eta_{c}m(#phi,#phi), (MC PID);E, GeV",100,2.5,3.5);
	TH1F *h_mphi_pid=new TH1F("h_mphi_pid","#phi: m(K+ K-) (MC PID);E, GeV",100,0.95,1.5);
	
	TH1F *h_etac_4c=new TH1F("h_etac_4c","#eta_{c}m(#phi,#phi), 4C-fit;E, GeV",100,2.5,3.5);
	TH1F *h_etac_4c_refit=new TH1F("h_etac_4c_refit","#eta_{c}m(#phi,#phi), 4C-fit;E, GeV",100,2.5,3.5);
	TH1F *h_mphi_4c=new TH1F("h_mphi_4c","#phi: m(K+ K-) (4C-fit);E, GeV",100,0.95,1.5);
	
	TH1F *h_etac_vtx=new TH1F("h_etac_vtx","#eta_{c}m(#phi,#phi), Vertex fit;E, GeV",100,2.5,3.5);
	TH1F *h_mphi_vtx=new TH1F("h_mphi_vtx","#phi: m(K+ K-) (Vertex fit);E, GeV",100,0.95,1.5);
			
	TH1F *h_etac_vtx_2=new TH1F("h_etac_vtx_2","#eta_{c}m(#phi,#phi), Vertex fit;E, GeV",100,2.5,3.5);
	TH1F *h_mphi_vtx_2=new TH1F("h_mphi_vtx_2","#phi: m(K+ K-) (Vertex fit);E, GeV",100,0.95,1.5);

	TH1F *h_etac_phimass_4c=new TH1F("h_etac_phimass_4c","#eta_{c}m(#phi,#phi), (cut on #phi mass);E,
	GeV",100,2.8,3.2);
	TH1F *h_mphi_final_4c=new TH1F("h_mphi_final_4c","#phi: m(K+ K-);E, GeV",100,0.95,1.1);

	TH1F *h_etac_phimass_vtx=new TH1F("h_etac_phimass_vtx","#eta_{c}m(#phi,#phi), (cut on #phi mass);E,
	GeV",100,2.8,3.2);
	TH1F *h_mphi_final_vtx=new TH1F("h_mphi_final_vtx","#phi: m(K+ K-);E, GeV",100,0.95,1.1);

	TH1F *h_etac_phimass_vtx_2=new TH1F("h_etac_phimass_vtx_2","#eta_{c}m(#phi,#phi), (cut on #phi mass);E,
	GeV",100,2.8,3.2);
	TH1F *h_mphi_final_vtx_2=new TH1F("h_mphi_final_vtx_2","#phi: m(K+ K-);E, GeV",100,0.95,1.1);

	TH1F *h_etac_phimassfit=new TH1F("h_etac_phimassfit","#eta_{c}m(#phi,#phi);E, GeV",100,2.8,3.2);
	TH1F *h_mphi_final_massfit=new TH1F("h_mphi_final_massfit","#phi: m(K+ K-);E, GeV",100,0.95,1.1);

	TH1F *nc=new TH1F("nc","n charged",20,0,20);

	TH1F *h_chi2_4c=new TH1F("h_chi2_4c","#chi^{2} 4C-fit;#chi^{2}/N_{df}",100,0,100);
	TH1F *h_chi2b_4c=new TH1F("h_chi2b_4c","#chi^{2} 4C-fit;#chi^{2}/N_{df}",100,0,100);
	TH1F *h_chi2_vtx=new TH1F("h_chi2_vtx","#chi^{2} vertex;#chi^{2}/N_{df}",100,0,100);
	TH1F *h_chi2b_vtx=new TH1F("h_chi2b_vtx","#chi^{2} vertex;#chi^{2}/N_{df}",100,0,100);
	
	TH1F *h_chi2_prefit=new TH1F("h_chi2_prefit","#chi^{2} prefit;#chi^{2}",100,0,100);
	
	TH1F *h_chi2_mass=new TH1F("h_chi2_mass","#chi^{2} Mass constraint fit;#chi^{2}",100,0,100);
	
	TH2F *hvpos = new TH2F("hvpos","(x,y) projection of fitted decay
	vertex",100,-5,5,100,-5,5);
	TH1F *hvzpos = new TH1F("hvzpos","z position of fitted decay vertex",100,-10,10);
	TH1F *hvtxresX = new TH1F("hvtxresX","X resolution of fitted decay
	vertex",100,-0.1,0.1);
	TH1F *hvtxresY = new TH1F("hvtxresY","Y resolution of fitted decay
	vertex",100,-0.1,0.1);
	TH1F *hvtxresZ = new TH1F("hvtxresZ","Z resolution of fitted decay
	vertex",100,-0.1,0.1);
	
	TH2F *h_theta_p = new TH2F("h_theta_p","Theta vs p",100,0,180,100,0.,3);
	TH2F *h_theta_p_nr = new TH2F("h_theta_p_nr","Theta vs p",100,0,180,100,0.,3); // not reconstructed
	
	TH2F *h_dp_p = new TH2F("h_dp_p","Delta p vs p",100,-3.,3,100,0,3);
	TH2F *h_dp_theta = new TH2F("h_dp_theta","Delta p vs theta",100,-3.,3,100,0,180);
	TH1F *h_dp=new TH1F("h_dp","Delta p",100,-3.,3.);
	TH1F *h_dp_low=new TH1F("h_dp_low","Delta p",100,-3.,3);
	TH1F *h_dp_high=new TH1F("h_dp_high","Delta p",100,-3.,3);

	// Number of events in file and number of reconstructed eta_c to store in root file
	TH1F *n_events=new TH1F("n_events","total number of events",1,0,1);
	TH1F *n_etac_4c=new TH1F("n_etac_4c","number of reconstructed eta_c (4C-fit)",1,0,1);
	TH1F *n_etac_vtx=new TH1F("n_etac_vtx","number of reconstructed eta_c (Vertex fit)",1,0,1);
	TH1F *n_etac_vtx_2=new TH1F("n_etac_vtx_2","number of reconstructed eta_c (Vertex fit)",1,0,1);
	
	/////////////////////// Analysis /////////////////////////////////
	
	// Set up the Analysis
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	//FairLogger::GetLogger()->SetLogScreenLevel("DEBUG4");
	//FairLogger::GetLogger()->SetLogVerbosityLevel("MEDIUM");
	FairRunAna* fRun = new FairRunAna();
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	
	fRun->SetInputFile(simFile);
	fRun->AddFriend(recoFile);
	fRun->AddFriend(pidFile);
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(parFile.Data());
	rtdb->setFirstInput(parIO);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(outFile.Data());
	fRun->Init();  

	TFile *out = TFile::Open(outFile,"RECREATE");
	
	PndAnalysis* theAnalysis = new PndAnalysis();

	// the candidates lists we need
	TCandList p1, p2, phi1, phi1_pid, phi1_massfit, phi1_massfit_cut;
	TCandList etac, etac_pid,  etac_nocut, etac_massfit, etac_massfit_cut;
	TCandList etac_vtx;
	TCandList mctracks;

	int evts = theAnalysis->GetEntries();
	if (nevts>0 && nevts<evts) evts=nevts; 
	n_events->SetBinContent(1,evts);
	
	TPidMassSelector *phiMassSel=new TPidMassSelector("phi",1.02,0.2);
	TPidPlusSelector *kplusSel=new TPidPlusSelector("kplus");
	TPidMinusSelector *kminusSel=new TPidMinusSelector("kminus");

	int n_reco_4c=0, n_reco_vtx=0, n_reco_vtx_2=0;

	TLorentzVector ini(0,0,3.6772,4.7333);

	int i=0,j=0, k=0, l=0;

	Double_t mc_mom, mc_theta, rec_mom, rec_theta, delta_p;
	int ievt;
	while (theAnalysis->GetEvent() && ievt++<evts)
	{    
		etac_vtx.Cleanup();
		phi1_massfit.Cleanup();
		phi1_massfit_cut.Cleanup();
		if ((ievt%100)==0)  cout<<"evt " << ievt << "\n";

		theAnalysis->FillList(mctracks,"McTruth");    
		theAnalysis->FillList(p1,"Charged");
		theAnalysis->FillList(p2,"Charged");
		
		p1.Select(kplusSel);
		p2.Select(kminusSel);
		
		int n_removed=0;
		int ii=0;
		for (Int_t l=0;l<p1.GetLength();l++){
			ii=l-n_removed;
			if((p1[ii].GetMicroCandidate().GetSttHits())==0){
				p1.Remove(p1[ii]);
				n_removed++;
			}
		}
		
		int n_removed=0;
		int ii=0;
		for (Int_t l=0;l<p2.GetLength();l++){
			ii=l-n_removed;
			if((p2[ii].GetMicroCandidate().GetSttHits())==0){
				p2.Remove(p2[ii]);
				n_removed++;
			}
		}

		int nchrg=p1.GetLength()+p2.GetLength();
		nc->Fill(nchrg);

		for (j=0;j<p1.GetLength();++j) {
			p1[j].SetMass(TRho::Instance()->GetPDG()->GetParticle(321)->Mass());
		}
		for (j=0;j<p2.GetLength();++j) {
			p2[j].SetMass(TRho::Instance()->GetPDG()->GetParticle(321)->Mass());
		}
		
		phi1.Combine(p1,p2);
		
		for (j=0;j<phi1.GetLength();++j) h_mphi_nocuts->Fill(phi1[j].M());
		
		phi1.Select(phiMassSel);
		etac_nocut.Combine(phi1,phi1);
		
		for (l=0;l<etac_nocut.GetLength();++l) {
			h_etac_nocut->Fill(etac_nocut[l].M());
		}

		TVector3 mcVertex, mcD1Vertex, mcD2vertex;
		FairMCEventHeader* evthead =
		(FairMCEventHeader*)(FairRootManager::Instance()->GetObject("EvtHeader"));
		evthead->GetVertex(mcVertex);
		
		// MC PID
		// Leave only kaons in particle lists
		if (usePID)
		{
		int n_removed=0;
		int ii=0;
		for (l=0;l<p1.GetLength();++l) {
			ii=l-n_removed;
			Int_t mcIndex = p1[ii].GetMicroCandidate().GetMcIndex();
			if (mcIndex>-1){
				if ((mctracks[mcIndex].PdgCode()!=321))
				{
					p1.Remove(p1[ii]);
					n_removed++;
				}
			}
			else
			{
				std::cout<<"stt h: " << p1[ii].GetMicroCandidate().GetSttHits() << std::endl;
				std::cout<<"Kaon list 1, element "<<l<<" has no assosiated mcTRack"<<std::endl;
				p1.Remove(p1[ii]);
				n_removed++;
			}
		}
		n_removed=0;
		ii=0;
		for (l=0;l<p2.GetLength();++l) {
			ii=l-n_removed;
			Int_t mcIndex = p2[ii].GetMicroCandidate().GetMcIndex();
			if (mcIndex>-1){
				if ((mctracks[mcIndex].PdgCode()!=-321))
				{
					p2.Remove(p2[ii]);
					n_removed++;
				} 
			}
			else
			{ 
				std::cout<<"stt h: " << p2[ii].GetMicroCandidate().GetSttHits() << std::endl;
				std::cout<<"Kaon list 2, element "<<l<<" has no assosiated mcTRack"<<std::endl;
				p2.Remove(p2[ii]);
				n_removed++;
			}
		}
		}

		phi1_pid.Combine(p1,p2);
		
		for (j=0;j<phi1_pid.GetLength();++j) h_mphi_pid->Fill(phi1_pid[j].M());
		phi1_pid.Select(phiMassSel);
		etac.Combine(phi1_pid,phi1_pid); 
		
		for (l=0;l<etac.GetLength();++l) {
			h_etac_pid->Fill(etac[l].M());
		}

		////////////// 4C-fit ///////////////
		int best_i=0;
		double best_chi2=10000;
		TCandidate *ccfit = new TCandidate();
		double m_phi1, m_phi2;
		for (l=0;l<etac.GetLength();++l) {
			Pnd4CFitter fitter(etac[l],ini);
			fitter.FitConserveMasses();
			double chi2=fitter.GetChi2();
			if (chi2<best_chi2)
			{
				best_chi2=chi2;
				best_i = l;
				ccfit = (TCandidate*)fitter.FittedCand(etac[l]);
				TCandidate *phi1best = fitter.FittedCand(*(etac[l].Daughter(0)));
				TCandidate *phi2best = fitter.FittedCand(*(etac[l].Daughter(1)));
				TCandidate *k1best = fitter.FittedCand(*(etac[l].Daughter(0)->Daughter(0)));
				TCandidate *k2best = fitter.FittedCand(*(etac[l].Daughter(0)->Daughter(1)));
				TCandidate *k3best = fitter.FittedCand(*(etac[l].Daughter(1)->Daughter(0)));
				TCandidate *k4best = fitter.FittedCand(*(etac[l].Daughter(1)->Daughter(1)));
				TLorentzVector tlvk1 = k1best->P4();
				TLorentzVector tlvk2 = k2best->P4();
				TLorentzVector tlvk3 = k3best->P4();
				TLorentzVector tlvk4 = k4best->P4();
				m_phi1= (tlvk1+tlvk2).M();
				m_phi2= (tlvk3+tlvk4).M();
			}
			h_chi2_4c->Fill(chi2/9); // Ndf=3N-3=9
		}
		
		if(/*(best_chi2<270)&&*/(etac.GetLength()!=0))
		{
			h_chi2b_4c->Fill(best_chi2/9); // Ndf=3N-3=9
			h_etac_4c->Fill(etac[best_i].M());
			h_etac_4c_refit->Fill(ccfit->M());
			h_mphi_4c->Fill(m_phi1);
			h_mphi_4c->Fill(m_phi2);
			h_mphi_final_4c->Fill(m_phi1);
			h_mphi_final_4c->Fill(m_phi2);
			if (((m_phi1>1.02-0.02)&&(m_phi1<1.02+0.02))
				&& ((m_phi2>1.02-0.02)&&(m_phi2<1.02+0.02)))
			{
				h_etac_phimass_4c->Fill(etac[best_i].M());
				if ((etac[best_i].M()>2.9)&&(etac[best_i].M()<3.06))
					n_reco_4c++;
			}
		}
		
		////////////// Vertex fit /////////////
		TCandidate *k1, *k2, *k3, *k4, *phi1tmp, *phi2tmp, *etac_tmp;
		
		//Combine 4 kaons directly to candidates
		for (j=0;j<etac.GetLength();++j)
		{
			phi1tmp=etac[j].Daughter(0);
			phi2tmp=etac[j].Daughter(1);

			k1=phi1tmp->Daughter(0);
			k2=phi1tmp->Daughter(1);
			k3=phi2tmp->Daughter(0);
			k4=phi2tmp->Daughter(1);
					
			etac_tmp=k1->Combine(*k2,*k3,*k4);
			etac_vtx.Add(*etac_tmp);
		}

		for (l=0;l<etac_vtx.GetLength();++l) {
			h_etac_pid->Fill(etac_vtx[l].M());
		}
		
		int best_i=0;
		double best_chi2=10000;
		TCandidate *etacfit_best=0;
		TCandidate *k1fit_best=0, *k2fit_best=0, *k3fit_best=0, *k4fit_best=0;
		TCandidate *phi1fit_best, *phi2fit_best;
		TVector3 bestPos;
		Float_t etacvtx_mass;
		for (j=0;j<etac_vtx.GetLength();++j)
		{
			PndKinVtxFitter vtxfitter(etac_vtx[j]);        // instantiate a vertex fitter
			vtxfitter.Fit();                          // do the vertex fit

			TCandidate *etacfit=vtxfitter.FittedCand(etac_vtx[j]);  // request the fitted EtaC candidate
			TVector3 etacVtx=etacfit->Pos();                    // and the decay vertex position
			double chi2_vtx=vtxfitter.GlobalChi2();
			h_chi2_vtx->Fill(chi2_vtx/5); // Number degree of freedom 2N-3=5
			// plot mass and vtx x,y projection after fit
			hvpos->Fill(etacVtx.X(),etacVtx.Y());
			hvzpos->Fill(etacVtx.Z());
			if(chi2_vtx<best_chi2)
			{
				best_chi2=chi2_vtx;
				best_i=l;
				etacfit_best=etacfit;
				k1fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(0)));
				k2fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(1)));
				k3fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(2)));
				k4fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(3)));
				etacvtx_mass = etacfit_best->M();
				bestPos = etacfit->Pos(); 
			}
		}
		
		if(/*(best_chi2<150)&&*/(etac.GetLength()!=0)&&(k1fit_best!=0)&&(k3fit_best!=0))
		{
			//h_etac_vtx->Fill(etacfit_best->M());
			h_chi2b_vtx->Fill(best_chi2/5);
			h_etac_vtx->Fill(etacvtx_mass);
			phi1fit_best=k1fit_best->Combine(*k2fit_best);
			phi2fit_best=k3fit_best->Combine(*k4fit_best);
			double m_phi1=phi1fit_best->M();
			double m_phi2=phi2fit_best->M();
			h_mphi_vtx->Fill(m_phi1);
			h_mphi_vtx->Fill(m_phi2);
			h_mphi_final_vtx->Fill(m_phi1);
			h_mphi_final_vtx->Fill(m_phi2);
			hvtxresX->Fill(mcVertex.X()-bestPos.X());
			hvtxresY->Fill(mcVertex.Y()-bestPos.Y());
			hvtxresZ->Fill(mcVertex.Z()-bestPos.Z());

			if (((m_phi1>1.02-0.02)&&(m_phi1<1.02+0.02))&& 
				((m_phi2>1.02-0.02)&&(m_phi2<1.02+0.02)))
			{
				h_etac_phimass_vtx->Fill(etacvtx_mass);
				if ((etacvtx_mass>2.9)&&(etacvtx_mass<3.06))
					n_reco_vtx++;
			}

		}
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////// Prefit selection ////////////////////////////////////
		int best_i=0;
		double best_chi2_prefit=1e6;
		TCandidate *etacprefit_best=0;
		TCandidate *k1prefit_best=0, *k2prefit_best=0, *k3prefit_best=0, *k4prefit_best=0;
		TCandidate *phi1prefit_best, *phi2prefit_best;
		Double_t mphi1, mphi2, metac;
		Double_t sigma_etac=0.0331, sigma_phi=0.00392;
		Double_t m_etac_pdg=2.98, m_phi_pdg=1.02;
		
		for (j=0;j<etac_vtx.GetLength();++j)
		{
			metac=etac_vtx[j].M();
			k1=etac_vtx[j].Daughter(0);
			k2=etac_vtx[j].Daughter(1);
			k3=etac_vtx[j].Daughter(2);
			k4=etac_vtx[j].Daughter(3);
			phi1prefit_best=k1->Combine(*k2);
			phi2prefit_best=k3->Combine(*k4);
			m_phi1=phi1prefit_best->M();
			m_phi2=phi2prefit_best->M();
		
			double chi2_prefit=pow(metac-m_etac_pdg,2)/pow(sigma_etac,2)+
			pow(m_phi1-m_phi_pdg,2)/pow(sigma_phi,2)+pow(m_phi1-m_phi_pdg,2)/pow(sigma_phi,2);
			h_chi2_prefit->Fill(chi2_prefit);

			if(chi2_prefit<best_chi2_prefit)
			{
				best_chi2_prefit=chi2_prefit;
				best_i=j;
				etacprefit_best=etac_vtx[j];
				etacvtx_mass=etacprefit_best->M();
			}
		}
		
		if (etacprefit_best!=0)
		{
			PndKinVtxFitter vtxfitter(*etacprefit_best);        // instantiate a vertex fitter
			vtxfitter.Fit();                          // do the vertex fit
			TCandidate *etacfit=vtxfitter.FittedCand(*etacprefit_best);
			k1prefit_best=vtxfitter.FittedCand(*(etacfit->Daughter(0)));
			k2prefit_best=vtxfitter.FittedCand(*(etacfit->Daughter(1)));
			k3prefit_best=vtxfitter.FittedCand(*(etacfit->Daughter(2)));
			k4prefit_best=vtxfitter.FittedCand(*(etacfit->Daughter(3)));
			phi1prefit_best=k1prefit_best->Combine(*k2prefit_best);
			phi2prefit_best=k3prefit_best->Combine(*k4prefit_best);
			double m_phi1=phi1prefit_best->M();
			double m_phi2=phi2prefit_best->M();
			h_mphi_vtx_2->Fill(m_phi1);
			h_mphi_vtx_2->Fill(m_phi2);
			h_mphi_final_vtx_2->Fill(m_phi1);
			h_mphi_final_vtx_2->Fill(m_phi2);
			if (((m_phi1>1.02-0.02)&&(m_phi1<1.02+0.02))&& 
				((m_phi2>1.02-0.02)&&(m_phi2<1.02+0.02)))
			{
				h_etac_phimass_vtx_2->Fill(etacvtx_mass);
				if ((etacvtx_mass>2.9)&&(etacvtx_mass<3.06))
					n_reco_vtx_2++;
			}
		}
	}
	
	std::cout<<"Number of reconstructed eta_c (4C) = "<<n_reco_4c<<std::endl;
	std::cout<<"Number of reconstructed eta_c (Vertex fit)= "<<n_reco_vtx<<std::endl;
	std::cout<<"Number of reconstructed eta_c (Vertex fit, preselection)= "<<n_reco_vtx_2<<std::endl;
	n_etac_4c->SetBinContent(1,n_reco_4c);
	n_etac_vtx->SetBinContent(1,n_reco_vtx);
	n_etac_vtx_2->SetBinContent(1,n_reco_vtx_2);

	out->cd();
	n_etac_4c->Write();
	n_etac_vtx->Write();
	n_etac_vtx_2->Write();
	n_events->Write();
	h_etac_nocut->Write();
	h_etac_pid->Write();
	h_etac_vtx->Write();
	h_etac_vtx_2->Write();
	h_etac_4c->Write();
	h_etac_4c_refit->Write();
	h_etac_phimass_4c->Write();
	h_etac_phimass_vtx->Write();
	h_etac_phimass_vtx_2->Write();
	h_etac_phimassfit->Write();
	
	h_mphi_nocuts->Write();
	h_mphi_pid->Write();
	h_mphi_vtx->Write();
	h_mphi_vtx_2->Write();
	h_mphi_4c->Write();
	h_mphi_final_4c->Write();
	h_mphi_final_vtx->Write();
	h_mphi_final_vtx_2->Write();
	h_mphi_final_massfit->Write();

	nc->Write();

	h_chi2_4c->Write();
	h_chi2b_4c->Write();
	h_chi2_vtx->Write();
	h_chi2b_vtx->Write();
	h_chi2_mass->Write();
	h_chi2_prefit->Write();
	hvzpos->Write();
	hvpos->Write();
	
	hvtxresX->Write();
	hvtxresY->Write();
	hvtxresZ->Write();

	out->Save();

	timer.Stop();
	Double_t rtime = timer.RealTime();
	Double_t ctime = timer.CpuTime();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

}

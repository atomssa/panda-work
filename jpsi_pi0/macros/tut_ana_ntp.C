class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;
class RhoTuple;

void tut_ana_ntp(int nevts=0)
{
	// *** some variables
	int i=0,j=0, k=0, l=0;
	gStyle->SetOptFit(1011);
	
	// *** the output file for FairRunAna
	TString OutFile="output.root";  
					
	// *** the files coming from the simulation
	TString inPidFile  = "pid_complete.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile  = "simparams.root";
	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";	
	
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* fRun = new FairRunAna();
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	fRun->SetInputFile(inPidFile);
	
	// *** setup parameter database 	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(inParFile);
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(OutFile);
	fRun->Init(); 
	
	// *** create an output file for all histograms
	TFile *out = TFile::Open("output_ana.root","RECREATE");
	
	// *** create ntuples for J/psi and psi(2S)
	RhoTuple *njpsi = new RhoTuple("njpsi","J/psi Analysis");
	RhoTuple *npsip = new RhoTuple("npsip","psi' Analysis");
	
	// *** create some columns which might not be filled sometimes
	njpsi->Column("tjpsim",		0.0f, -999.9f);
	njpsi->Column("tjpsip",		0.0f, -999.9f);
	njpsi->Column("tjpsitht",	0.0f, -999.9f);
	
	npsip->Column("tpsim",		0.0f, -999.9f);
	npsip->Column("tpsip",		0.0f, -999.9f);
	npsip->Column("tpsitht",	0.0f, -999.9f);

	//
	// Now the analysis stuff comes...
	//
	
	
	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** RhoCandLists for the analysis
	RhoCandList muplus, muminus, piplus, piminus, jpsi, psi2s;
	
	// *** Mass selector for the jpsi cands
	double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi
	RhoMassParticleSelector *jpsiMassSel=new RhoMassParticleSelector("jpsi",m0_jpsi,1.0);
	
	// *** Pid Selection Algorithms
	TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";
	
	// *** the lorentz vector of the initial psi(2S)
	TLorentzVector ini(0, 0, 6.231552, 7.240065);
	
	// ***
	// the event loop
	// ***
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
				
		// *** Select with no PID info ('All'); type and mass are set 		
		theAnalysis->FillList(muplus,  "MuonAllPlus",	pidSelection);
		theAnalysis->FillList(muminus, "MuonAllMinus",	pidSelection);
		theAnalysis->FillList(piplus,  "PionAllPlus",	pidSelection);
		theAnalysis->FillList(piminus, "PionAllMinus",	pidSelection);
		
		// *** combinatorics for J/psi -> mu+ mu-
		jpsi.Combine(muplus, muminus);
		jpsi.SetType(443);
		
		// *** some mass pre selection
		//jpsi.Select(jpsiMassSel);
		
		// ***
		// *** do all kind of analysis and store in N-tuple
		// ***
		for (j=0;j<jpsi.GetLength();++j) 
		{
			// get daughters
			RhoCandidate *mup = jpsi[j]->Daughter(0);
			RhoCandidate *mum = jpsi[j]->Daughter(1);
			PndPidCandidate *mup_rec = (PndPidCandidate*)mup->GetRecoCandidate();
			PndPidCandidate *mum_rec = (PndPidCandidate*)mum->GetRecoCandidate();
			
			// get truth information 
			bool mct = theAnalysis->McTruthMatch(jpsi[j]);
			RhoCandidate *true_jpsi = jpsi[j]->GetMcTruth();
			
			// perform vertex fitter
			PndKinVtxFitter vtxfitter(jpsi[j]);	// instantiate a vertex fitter
			vtxfitter.Fit();
			
			RhoCandidate *fitvtx_jpsi = jpsi[j]->GetFit();
			double chi2_vtx = vtxfitter.GetChi2();	// access chi2 of fit
			double prob_vtx = vtxfitter.GetProb();	// access probability of fit
			TVector3 vtxpos(-999.,-999.,-999.);
			if (fitvtx_jpsi) vtxpos = fitvtx_jpsi->Daughter(0)->Pos();
			
			// perform mass fit
			PndKinFitter mfitter(jpsi[j]);		// instantiate the PndKinFitter in psi(2S)
			mfitter.AddMassConstraint(m0_jpsi);	// add the mass constraint
			mfitter.Fit();						// do fit
			
			RhoCandidate *fitmass_jpsi = jpsi[j]->GetFit();
			double chi2_mass = mfitter.GetChi2();	// get chi2 of fit
			double prob_mass = mfitter.GetProb();	// access probability of fit
			
			// *** now write ntuple information
			
			// *** general event info
			njpsi->Column("ev",			(Float_t) i,							-999.9f);
			njpsi->Column("cand",		(Float_t) j,							-999.9f);
			
			// *** basic J/psi info
			njpsi->Column("jpsim",		(Float_t) jpsi[j]->M(),					-999.9f); 
			njpsi->Column("jpsip",		(Float_t) jpsi[j]->P(),					-999.9f); 
			njpsi->Column("jpsipt",		(Float_t) jpsi[j]->P3().Pt(),			-999.9f); 
			njpsi->Column("jpsitht",	(Float_t) jpsi[j]->P3().Theta(),		-999.9f); 
			njpsi->Column("jpsimissm",	(Float_t) (ini-(jpsi[j]->P4())).M(),	-999.9f);
			
			// *** MC truth info
			njpsi->Column("mct",		(Float_t) mct,							-999.9f);
			if (true_jpsi)
			{
				njpsi->Column("tjpsim",	(Float_t) true_jpsi->M(),				-999.9f);
				njpsi->Column("tjpsip",	(Float_t) true_jpsi->M(),				-999.9f);
				njpsi->Column("tjpsitht",(Float_t) true_jpsi->P3().Theta(),		-999.9f);
			}

			// *** fitting info
			njpsi->Column("jpsimvtx",	(Float_t) fitvtx_jpsi->M(),				-999.9f); 
			njpsi->Column("chi2vtx",	(Float_t) chi2_vtx,						-999.9f); 
			njpsi->Column("probvtx",	(Float_t) prob_vtx,						-999.9f); 
			njpsi->Column("vtxx",		vtxpos.X(),								-999.9f);
			njpsi->Column("vtxy",		vtxpos.Y(),								-999.9f);
			njpsi->Column("vtxz",		vtxpos.Z(),								-999.9f);
			
			njpsi->Column("jpsimmass",	(Float_t) fitmass_jpsi->M(),			-999.9f); 
			njpsi->Column("chi2mass",	(Float_t) chi2_mass,					-999.9f); 
			njpsi->Column("probmass",	(Float_t) prob_mass,					-999.9f); 
			
			// *** kinematic info of daughters
			njpsi->Column("mupp",		(Float_t) mup->P(),						-999.9f);
			njpsi->Column("muppt",		(Float_t) mup->P3().Pt(),				-999.9f);
			njpsi->Column("muptht",		(Float_t) mup->P3().Theta(),			-999.9f);
			
			njpsi->Column("mump",		(Float_t) mum->P(),						-999.9f);
			njpsi->Column("mumpt",		(Float_t) mum->P3().Pt(),				-999.9f);
			njpsi->Column("mumtht",		(Float_t) mum->P3().Theta(),			-999.9f);
			
			// *** PID info of daughters
			njpsi->Column("muppid",		(Float_t) mup->GetPidInfo(1),			-999.9f);
			njpsi->Column("mumpid",		(Float_t) mum->GetPidInfo(1),			-999.9f);
			
			// *** and finally FILL Ntuple
			njpsi->DumpData();
			
		}
		
		// *** combinatorics for psi(2S) -> J/psi pi+ pi-
		psi2s.Combine(jpsi, piplus, piminus);
		psi2s.SetType(30443);
		
		for (j=0;j<psi2s.GetLength();++j) 
		{
			// get daughters
			RhoCandidate *jp =  psi2s[j]->Daughter(0);
			RhoCandidate *pip = psi2s[j]->Daughter(1);
			RhoCandidate *pim = psi2s[j]->Daughter(2);
			
			PndPidCandidate *pip_rec = (PndPidCandidate*)pip->GetRecoCandidate();
			PndPidCandidate *pim_rec = (PndPidCandidate*)pim->GetRecoCandidate();
			
			// get truth information 
			bool mct = theAnalysis->McTruthMatch(psi2s[j]);
			RhoCandidate *true_psi = psi2s[j]->GetMcTruth();
			
			// do 4C fit
			PndKinFitter fitter(psi2s[j]);	// instantiate the kin fitter in psi(2S)
			fitter.Add4MomConstraint(ini);	// set 4 constraint
			fitter.Fit();		            // do fit
			RhoCandidate *fit4c_jpsi = psi2s[j]->Daughter(0)->GetFit();	// get fitted J/psi
			
			double chi2_4c = fitter.GetChi2();	// get chi2 of fit
			double prob_4c = fitter.GetProb();	// access probability of fit
			
			// *** general event info
			npsip->Column("ev",		(Float_t) i,							-999.9f);
			npsip->Column("cand",	(Float_t) j,							-999.9f);
			
			// *** basic psi(2s) info
			npsip->Column("psim",	(Float_t) psi2s[j]->M(),				-999.9f); 
			npsip->Column("psip",	(Float_t) psi2s[j]->P(),				-999.9f); 
			npsip->Column("psipt",	(Float_t) psi2s[j]->P3().Pt(),			-999.9f); 
			npsip->Column("psitht",	(Float_t) psi2s[j]->P3().Theta(),		-999.9f); 
			
			// *** basic J/psi info
			npsip->Column("jpsim",	(Float_t) psi2s[j]->M(),				-999.9f); 
			npsip->Column("jpsip",	(Float_t) psi2s[j]->P(),				-999.9f); 
			npsip->Column("jpsipt",	(Float_t) psi2s[j]->P3().Pt(),			-999.9f); 
			npsip->Column("jpsitht",(Float_t) psi2s[j]->P3().Theta(),		-999.9f); 
			
			npsip->Column("jpsim4c",(Float_t) fit4c_jpsi->M(),				-999.9f);
			
			// *** MC truth info
			npsip->Column("mct",	(Float_t) mct,							-999.9f);
			if (true_psi)
			{
				npsip->Column("tpsim",	(Float_t) true_psi->M(),			-999.9f);
				npsip->Column("tpsip",	(Float_t) true_psi->M(),			-999.9f);
				npsip->Column("tpsitht",(Float_t) true_psi->P3().Theta(),	-999.9f);
			}
			
			// *** kinematic info of daughters
			npsip->Column("pipp",	(Float_t) pip->P(),						-999.9f);
			npsip->Column("pippt",	(Float_t) pip->P3().Pt(),				-999.9f);
			npsip->Column("piptht",	(Float_t) pip->P3().Theta(),			-999.9f);
			
			npsip->Column("pimp",	(Float_t) pim->P(),						-999.9f);
			npsip->Column("pimpt",	(Float_t) pim->P3().Pt(),				-999.9f);
			npsip->Column("pimtht",	(Float_t) pim->P3().Theta(),			-999.9f);
			
			// *** PID info of daughters
			npsip->Column("pippid",	(Float_t) pip->GetPidInfo(2),			-999.9f);
			npsip->Column("pimpid",	(Float_t) pim->GetPidInfo(2),			-999.9f);
			
			// *** and finally FILL Ntuple
			npsip->DumpData();
		}		
		
	}
	
	// *** write out all the histos
	out->cd();

	njpsi->GetInternalTree()->Write();
	npsip->GetInternalTree()->Write();
		
	out->Save();
	
}

class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

void ana_jpsi(int nevts=0)
{
  // *** some variables
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

  // *** create some histograms
  TH1F *hjpsim_all = new TH1F("hjpsim_all","J/#psi mass (all)",200,0,5.0);
  TH1F *hjpsim_all_b = new TH1F("hjpsim_all_b","J/#psi mass (all, brem)",200,0,5.0);
  TH2F *hjpsim_all_corr = new TH2F("hjpsim_all_corr","J/#psi mass (brem vs no)",200,0,5.0,200,0,5.0);
  TH1F *hjpsim_ftm = new TH1F("hjpsim_ftm","J/#psi mass (full truth match)",200,0,5.0);
  TH1F *hjpsim_nm = new TH1F("hjpsim_nm","J/#psi mass (no truth match)",200,0,5.0);
  TH1F *hjpsim_diff = new TH1F("hjpsim_diff","J/#psi mass diff to truth",100,-2,2);
  TH2F *hm_vs_ephi = new TH2F("hm_vs_ephi","J/#psi mass (brem) vs. phi of electron",100,-180,180,100,0,5.0);
  TH2F *hm_vs_pphi = new TH2F("hm_vs_pphi","J/#psi mass (brem) vs. phi of positron",100,-180,180,100,0,5.0);
  TH2F *hm_vs_ethe = new TH2F("hm_vs_ethe","J/#psi mass (brem) vs. the of electron",100,0,180,100,0,5.0);
  TH2F *hm_vs_pthe = new TH2F("hm_vs_pthe","J/#psi mass (brem) vs. the of positron",100,0,180,100,0,5.0);
  TH2F *bad_ep_phi_vs_th = new TH2F("bad_ep_phi_vs_the","phi vs theta of dauth elecs for bad jpsi mass",100,0,180,100,-180,180);
  TH2F *good_ep_phi_vs_th = new TH2F("good_ep_phi_vs_the","phi vs theta of dauth elecs for good jpsi mass",100,0,180,100,-180,180);
  TH2F *all_ep_phi_vs_th = new TH2F("all_ep_phi_vs_the","phi vs theta of dauth elecs for all jpsi mass",100,0,180,100,-180,180);

  //
  // Now the analysis stuff comes...
  //
  // *** the data reader object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts= theAnalysis->GetEntries();

  // *** RhoCandLists for the analysis
  RhoCandList eb, pb, jpsib;
  RhoCandList e, p, jpsi;

  // *** Mass selector for the jpsi cands
  double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi
  RhoMassParticleSelector *jpsiMassSel=new RhoMassParticleSelector("jpsi",m0_jpsi,1.0);

  // *** the lorentz vector of the initial psi(2S)
  TLorentzVector ini(0, 0, 6.231552, 7.240065);

  // ***
  // the event loop
  // ***
  int i = 0;
  while (theAnalysis->GetEvent() && i++<nevts)
    {
      if ((i%100)==0) cout<<"evt " << i << endl;

      // *** Select with no PID info ('All'); type and mass are set
      theAnalysis->FillList(e, "ElectronAllMinus");
      theAnalysis->FillList(p, "ElectronAllPlus");
      theAnalysis->FillList(eb, "BremElectronAllMinus");
      theAnalysis->FillList(pb, "BremElectronAllPlus");

      // *** combinatorics for J/psi -> mu+ mu-
      jpsi.Combine(e, p);
      jpsib.Combine(eb, pb);

      // ***
      // *** do the TRUTH MATCH for jpsi
      // ***
      //jpsi.SetType(443);
      //jpsib.SetType(443);

      for (int j=0;j<jpsi.GetLength();++j)
	{
	  hjpsim_all->Fill( jpsi[j]->M() );
	  hjpsim_all_b->Fill( jpsib[j]->M() );
	  hjpsim_all_corr->Fill(jpsi[j]->M(), jpsib[j]->M());

	  int ch0 = jpsib[j]->Daughter(0)->Charge();
	  double eth = jpsib[j]->Daughter(ch0<0?0:1)->P3().Theta()*TMath::RadToDeg();
	  double eph = jpsib[j]->Daughter(ch0<0?0:1)->P3().Phi()*TMath::RadToDeg();
	  double pth = jpsib[j]->Daughter(ch0<0?1:0)->P3().Theta()*TMath::RadToDeg();
	  double pph = jpsib[j]->Daughter(ch0<0?1:0)->P3().Phi()*TMath::RadToDeg();

	  hm_vs_ephi->Fill(eph, jpsib[j]->M());
	  hm_vs_pphi->Fill(pph, jpsib[j]->M());
	  hm_vs_ethe->Fill(eth, jpsib[j]->M());
	  hm_vs_pthe->Fill(pth, jpsib[j]->M());
	  if (jpsib[j]->M()>3.5) {
	    bad_ep_phi_vs_th->Fill(eth,eph);
	    bad_ep_phi_vs_th->Fill(pth,pph);
	  } else {
	    good_ep_phi_vs_th->Fill(eth,eph);
	    good_ep_phi_vs_th->Fill(pth,pph);
	  }
	  all_ep_phi_vs_th->Fill(eth,eph);
	  all_ep_phi_vs_th->Fill(pth,pph);
	  //if (theAnalysis->McTruthMatch(jpsi[j]))
	  //  {
	  //    hjpsim_ftm->Fill( jpsi[j]->M() );
	  //    hjpsim_diff->Fill( jpsi[j]->GetMcTruth()->M() - jpsi[j]->M() );
	  //  }
	  //else
	  //  hjpsim_nm->Fill( jpsi[j]->M() );
	}
    }

  // *** write out all the histos
  out->cd();

  hjpsim_all->Write();
  hjpsim_all_b->Write();
  hjpsim_all_corr->Write();
  //hjpsim_ftm->Write();
  //hjpsim_diff->Write();
  //hjpsim_nm->Write();

  hm_vs_ephi->Write();
  hm_vs_pphi->Write();
  hm_vs_ethe->Write();
  hm_vs_pthe->Write();
  bad_ep_phi_vs_th->Write();
  good_ep_phi_vs_th->Write();
  all_ep_phi_vs_th->Write();

  out->Save();

}

class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

// *** routine to only keep PID matched candidates in list
int SelectTruePid(PndAnalysis *ana, RhoCandList &l)
{
  int removed = 0;

  for (int ii=l.GetLength()-1;ii>=0;--ii)
    {
      if ( !(ana->McTruthMatch(l[ii])) )
	{
	  l.Remove(l[ii]);
	  removed++;
	}
    }

  return removed;
}


void ana(int nevts=0)
{
  // *** some variables
  //int i=0,j=0, k=0, l=0;
  gStyle->SetOptFit(1011);

  // *** the output file for FairRunAna
  TString OutFile="output.root";

  // *** the files coming from the simulation
  TString inPidFile  = "output/pid/pbar_p_jpsi_pi0_100.root";    // this file contains the PndPidCandidates and McTruth
  TString inParFile  = "output/par/pbar_p_jpsi_pi0_100.root";

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
  TH1F *hjpsim_all = new TH1F("hjpsim_all","J/#psi mass (all)",200,0,4.5);
  TH1F *hpi0m_all  = new TH1F("hpi0m_all","#pi^{0} mass (all)",200,0,0.2);
  TH1F *hjpsipi0m_all  = new TH1F("hjpsipi0m_all","J/#psi-#pi^{0} mass (all)",200,0,5);

  TH1F *hjpsim_ftm = new TH1F("hjpsim_ftm","J/#psi mass (full truth match)",200,0,4.5);
  TH1F *hpi0m_ftm  = new TH1F("hpi0m_ftm","#pi^{0} mass (full truth match)",200,0,0.2);
  TH1F *hjpsipi0m_ftm  = new TH1F("hjpsipi0m_ftm","J/#psi-#pi^{0} mass (full truth match)",200,0,5);

  TH1F *hjpsim_nm = new TH1F("hjpsim_nm","J/#psi mass (no truth match)",200,0,4.5);
  TH1F *hpi0m_nm  = new TH1F("hpi0m_nm","#pi^{0} mass (no truth match)",200,0,0.2);
  TH1F *hjpsipi0m_nm  = new TH1F("hjpsipi0m_nm","J/#psi-#pi^{0} mass (no truth match)",200,0,5);

  TH1F *hjpsim_diff = new TH1F("hjpsim_diff","J/#psi mass diff to truth",100,-2,2);
  TH1F *hpi0m_diff  = new TH1F("hpi0m_diff","#pi^{0} mass diff to truth",100,-2,2);
  TH1F *hjpsipi0m_diff  = new TH1F("hjpsipi0m_diff","J/#psi-#pi^{0} mass diff to truth",100,-2,2);

  //TH1F *hjpsim_lpid = new TH1F("hjpsim_lpid","J/#psi mass (loose pid)",200,0,4.5);
  //TH1F *hpi0m_lpid  = new TH1F("hpi0m_lpid","#pi^{0} mass (loose pid)",200,0,5);
  //
  //TH1F *hjpsim_tpid = new TH1F("hjpsim_tpid","J/#psi mass (tight pid)",200,0,4.5);
  //TH1F *hpi0m_tpid  = new TH1F("hpi0m_tpid","#pi^{0} mass (tight pid)",200,0,5);
  //
  //TH1F *hjpsim_trpid = new TH1F("hjpsim_trpid","J/#psi mass (true pid)",200,0,4.5);
  //TH1F *hpi0m_trpid  = new TH1F("hpi0m_trpid","#pi^{0} mass (true pid)",200,0,5);

  // *** the data reader object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts= theAnalysis->GetEntries();

  // *** RhoCandLists for the analysis
  //RhoCandList muplus, muminus, piplus, piminus, jpsi, psi2s;
  RhoCandList ep, em, g1, g2, jpsi, pi0, tot;

  // *** Mass selector for the jpsi cands
  double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi
  RhoMassParticleSelector *jpsiMassSel=new RhoMassParticleSelector("jpsi",m0_jpsi,1.0);

  // *** the lorentz vector of the initial jpsi-pi0 system
  TLorentzVector ini(0, 0, 6.231552, 7.240065);

  // ***
  // the event loop
  // ***
  int i= 0;
  while (theAnalysis->GetEvent() && i++<nevts) {

    if ((i%100)==0) cout<<"evt " << i << endl;

    // *** Select with no PID info ('All'); type and mass are set
    theAnalysis->FillList(ep,  "ElectronAllPlus");
    theAnalysis->FillList(em, "ElectronAllMinus");
    theAnalysis->FillList(g1,  "Neutral");
    theAnalysis->FillList(g2, "Neutral");

    // *** combinatorics for J/psi -> mu+ mu-
    jpsi.Combine(ep, em);
    pi0.Combine(g1,g2);
    tot.Combine(pi0, jpsi);

    // ***
    // *** do the TRUTH MATCH for jpsi
    // ***
    jpsi.SetType(443);
    pi0.SetType(111);

    for (int j=0;j<jpsi.GetLength();++j) {
      hjpsim_all->Fill( jpsi[j]->M() );
      if (theAnalysis->McTruthMatch(jpsi[j])) {
	hjpsim_ftm->Fill( jpsi[j]->M() );
	hjpsim_diff->Fill( jpsi[j]->GetMcTruth()->M() - jpsi[j]->M() );
      } else {
	hjpsim_nm->Fill( jpsi[j]->M() );
      }
    }

    for (int j=0; j<pi0.GetLength(); ++j) {
      hpi0m_all->Fill( pi0[j]->M() );
      if ( theAnalysis->McTruthMatch(pi0[j]) ) {
	hpi0m_ftm->Fill( pi0[j]->M() );
	hpi0m_diff->Fill( pi0[j]->GetMcTruth()->M() - pi0[j]->M() );
      } else {
	hpi0m_nm->Fill( pi0[j]->M() );
      }
    }

    //for (int j=0; j<tot.GetLength(); ++j) {
    //  hpi0jpsim_all->Fill( tot[j]->M() );
    //  //if ( theAnalysis->McTruthMatch(pi0[j]) ) {
    //  //	hpi0jpsim_ftm->Fill( pi0[j]->M() );
    //  //	hpi0jpsim_diff->Fill( pi0[j]->GetMcTruth()->M() - pi0[j]->M() );
    //  //} else {
    //  //	hpi0jpsim_nm->Fill( pi0[j]->M() );
    //  //}
    //}




  }

  out->cd();

  hjpsim_all->Write();
  hpi0m_all->Write();
  hjpsipi0m_all->Write();

  hjpsim_ftm->Write();
  hpi0m_ftm->Write();
  hjpsipi0m_ftm->Write();

  hjpsim_nm->Write();
  hpi0m_nm->Write();
  hjpsipi0m_nm->Write();

  hjpsim_diff->Write();
  hpi0m_diff->Write();
  hjpsipi0m_diff->Write();

  out->Write();
  out->Close();

}

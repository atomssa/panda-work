class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

void theta(int nevts=0)
{
  // *** some variables
  gStyle->SetOptFit(1011);

  // *** the output file for FairRunAna
  TString OutFile="output.root";

  TString dir="/projet/panda/Ermias/tda/elec_flat/runall.50/";
  TString inPidFile  = dir+"pid_complete.root";    // this file contains the PndPidCandidates and McTruth
  TString inParFile  = dir+"simparams.root";

  //// *** the files coming from the simulation
  //TString inPidFile  = TString(gSystem->Getenv("DIR"))+"pid_complete.root";    // this file contains the PndPidCandidates and McTruth
  //TString inParFile  = TString(gSystem->Getenv("DIR"))+"simparams.root";

  cout << "inPid= " << inPidFile << endl;
  cout << "inPar= " << inParFile << endl;

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

  TFile *out = TFile::Open("output_ana.root","RECREATE");

  TH1F *eth_all = new TH1F("eth_all","#theta dist. of e^{-}",100,0,90);
  TH1F *eth_eid = new TH1F("eth_eid","#theta dist. of e^{-}",100,0,90);
  eth_eid->SetLineColor(2);

  TH1F *emom_all = new TH1F("emom_all","mom dist. of e^{-}",100,0,5);
  TH1F *emom_eid = new TH1F("emom_eid","mom dist. of e^{-}",100,0,5);
  emom_eid->SetLineColor(2);

  TH1F *pth_all = new TH1F("pth_all","#theta dist. of e^{+}",100,0,90);
  TH1F *pth_eid = new TH1F("pth_eid","#theta dist. of e^{+}",100,0,90);
  pth_eid->SetLineColor(2);

  TH1F *pmom_all = new TH1F("pmom_all","mom dist. of e^{+}",100,0,5);
  TH1F *pmom_eid = new TH1F("pmom_eid","mom dist. of e^{+}",100,0,5);
  pmom_eid->SetLineColor(2);

  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts= theAnalysis->GetEntries();

  RhoCandList e, p;
  RhoCandList e_eid, p_eid;

  int i = 0;
  while (theAnalysis->GetEvent() && i++<nevts) {
    if ((i%100)==0) cout<<"evt " << i << endl;

    // *** Select with no PID info ('All'); type and mass are set
    theAnalysis->FillList(e, "ElectronAllMinus");
    theAnalysis->FillList(p, "ElectronAllPlus");

    theAnalysis->FillList(e_eid, "ElectronVeryTightMinus", "PidAlgoEmcBayes");
    theAnalysis->FillList(p_eid, "ElectronVeryTightPlus", "PidAlgoEmcBayes");

    for (int ii = 0; ii < e.GetLength(); ++ii) {
      eth_all->Fill(e[ii]->P3().Theta()*TMath::RadToDeg());
      emom_all->Fill(e[ii]->P3().Mag());
    }
    for (int ii = 0; ii < p.GetLength(); ++ii) {
      pth_all->Fill(p[ii]->P3().Theta()*TMath::RadToDeg());
      pmom_all->Fill(p[ii]->P3().Mag());
    }

    for (int ii = 0; ii < e_eid.GetLength(); ++ii) {
      eth_eid->Fill(e_eid[ii]->P3().Theta()*TMath::RadToDeg());
      emom_eid->Fill(e_eid[ii]->P3().Mag());
    }
    for (int ii = 0; ii < p_eid.GetLength(); ++ii) {
      pth_eid->Fill(p_eid[ii]->P3().Theta()*TMath::RadToDeg());
      pmom_eid->Fill(p_eid[ii]->P3().Mag());
    }
  }

  TEfficiency *eeff = new TEfficiency(*eth_eid, *eth_all);
  eeff->SetMarkerStyle(20);
  eeff->SetMarkerSize(1);
  TEfficiency *peff = new TEfficiency(*pth_eid, *pth_all);
  peff->SetMarkerStyle(20);
  peff->SetMarkerSize(1);

  TEfficiency *eeff_mom = new TEfficiency(*emom_eid, *emom_all);
  eeff_mom->SetMarkerStyle(20);
  eeff_mom->SetMarkerSize(1);
  TEfficiency *peff_mom = new TEfficiency(*pmom_eid, *pmom_all);
  peff_mom->SetMarkerStyle(20);
  peff_mom->SetMarkerSize(1);

  TCanvas*tc =new TCanvas("tc","tc");
  tc->Divide(2,2);
  tc->cd(1);
  eth_all->Draw();
  eth_eid->Draw("same");
  tc->cd(2);
  pth_all->Draw();
  pth_eid->Draw("same");
  tc->cd(3);
  eeff->Draw("");
  tc->cd(4);
  peff->Draw("");

  TCanvas*tc2 =new TCanvas("tc2","tc2");
  tc2->Divide(2,2);
  tc2->cd(1);
  emom_all->Draw();
  emom_eid->Draw("same");
  tc2->cd(2);
  pmom_all->Draw();
  pmom_eid->Draw("same");
  tc2->cd(3);
  eeff_mom->Draw("");
  tc2->cd(4);
  peff_mom->Draw("");

  out->cd();

  eth_all->Write();
  pth_all->Write();
  eth_eid->Write();
  pth_eid->Write();

  emom_all->Write();
  pmom_all->Write();
  emom_eid->Write();
  pmom_eid->Write();

  eeff->Write();
  peff->Write();
  eeff_mom->Write();
  peff_mom->Write();

  out->Save();

}

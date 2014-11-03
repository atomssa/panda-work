class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;


void tut_ana_mclist(int nevts=0)
{
	// *** some variables
	int i=0,j=0, k=0, l=0;

	TString OutFile="output.root";

	// *** the files coming from the simulation
	TString inPidFile = "pid_complete.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile = "simparams.root";

	gStyle->SetOptFit(1011);

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

	//
	// Now the analysis stuff comes...
	//


	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();

	// *** RhoCandLists for the analysis
	RhoCandList mctruth;

	// ***
	// the event loop
	// ***
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		cout<<"****** Evt " << i << endl;

		// *** the MC Truth objects
		theAnalysis->FillList(mctruth,"McTruth");

		//
		// Print MC Truth list with mother-daughter relations
		//
		for (j=0;j<mctruth.GetLength();++j)
		{
			RhoCandidate *mcmother = mctruth[j]->TheMother();

		  	int muid = -1;
		  	if (mcmother) muid = mcmother->GetTrackNumber();

		  	cout << "Track "<< mctruth[j]->GetTrackNumber()<<" (PDG:"<<mctruth[j]->PdgCode() <<") has mother "<<muid;

		  	if (mctruth[j]->NDaughters()>0) cout <<" and daughter(s) ";
		  	for (k=0;k<mctruth[j]->NDaughters();++k) cout <<mctruth[j]->Daughter(k)->GetTrackNumber()<<"  ";

		  	cout<<endl;
		}
		cout <<endl;
	}

}

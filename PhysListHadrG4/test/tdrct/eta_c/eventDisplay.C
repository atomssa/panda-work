

eventDisplay()
{
    //-----User Settings:-----------------------------------------------
  TString  SimEngine      ="TGeant3"; 
  TString  InputFile     ="evt_points_stt.root";
  TString  DigiFile      ="evt_digi_stt.root";
  TString  RecoFile		 = "evt_reco_stt.root";
  TString  ParFile       ="evt_params_stt.root";
  //------------------------------------------------------------------


// Load basic libraries
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  gSystem->Load("libEve");
  gSystem->Load("libEventDisplay");
  gSystem->Load("libPndEventDisplay");

                                     
  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(InputFile.Data());
  fRun->AddFriend(DigiFile.Data());
  fRun->AddFriend(RecoFile.Data());
  fRun->SetOutputFile("tst.root");

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(ParFile.Data());
  rtdb->setFirstInput(parInput1);
  FairEventManager *fMan= new FairEventManager();

  FairGeane *Geane = new FairGeane();
   fRun->AddTask(Geane);
 
 
 //----------------------Traks and points -------------------------------------
  FairMCTracks *Track =  new FairMCTracks ("Monte-Carlo Tracks");
  FairMCPointDraw *MvdPoints =   new FairMCPointDraw ("MVDPoint",kBlue,  kFullSquare);
  FairMCPointDraw *EMCPoints =   new FairMCPointDraw ("EmcHit",kOrange,  kFullSquare);
  FairMCPointDraw *TofPoint =    new FairMCPointDraw ("TofPoint",kYellow,  kFullSquare);
  FairMCPointDraw *TofSciFPoint= new FairMCPointDraw ("TofSciFPoint",kTeal, kFullSquare);
  FairMCPointDraw *MdtPoint =    new FairMCPointDraw ("MdtPoint",kAzure, kFullSquare);
  FairMCPointDraw *PndDrcPoint = new FairMCPointDraw ("PndDrcPoint",kViolet, kFullSquare);
  FairMCPointDraw *PndDchPoint = new FairMCPointDraw ("PndDchPoint",kPink, kFullSquare);
  FairMCPointDraw *PndTpcPoint = new FairMCPointDraw ("PndTpcPoint",kCyan,  kFullSquare);
  FairMCPointDraw *PndSTTPoint = new FairMCPointDraw ("STTPoint",kMagenta, kFullSquare);
  FairMCPointDraw *PndGEMPoint = new FairMCPointDraw ("GEMPoint",kRed, kFullSquare);
  FairMCPointDraw *PndDskPoint = new FairMCPointDraw ("DskCerenkov",kGreen, kFullSquare);

  FairHitDraw *MvdHitsPixel = new FairHitDraw("MVDHitsPixel");
  FairHitDraw *MvdHitsStrip = new FairHitDraw("MVDHitsStrip");

  FairHitDraw *GemHit = new FairHitDraw("GEMHit");
  FairHitDraw *PndTpcCluster = new FairHitDraw("PndTpcCluster");

  PndTrackCandDraw *LheTrackCand = new PndTrackCandDraw("LheTrackCand");

  PndTrackDraw *LheGenTrack = new PndTrackDraw("LheGenTrack", kTRUE);

  fMan->AddTask(Track);
  fMan->AddTask(MvdPoints);
  fMan->AddTask(EMCPoints);   
  fMan->AddTask(TofPoint);   
  fMan->AddTask(TofSciFPoint);
  fMan->AddTask(MdtPoint);
  fMan->AddTask(PndDrcPoint);
  fMan->AddTask(PndDchPoint);
  fMan->AddTask(PndTpcPoint);
  fMan->AddTask(PndSTTPoint);
  fMan->AddTask(PndGEMPoint);
  fMan->AddTask(PndDskPoint);
  
  fMan->AddTask(MvdHitsPixel);
  fMan->AddTask(MvdHitsStrip);
  fMan->AddTask(GemHit);
  fMan->AddTask(PndTpcCluster);

  fMan->AddTask(LheTrackCand);
  fMan->AddTask(LheGenTrack);

  fMan->Init();                     

}

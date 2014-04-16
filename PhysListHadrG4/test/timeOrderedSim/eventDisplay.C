

eventDisplay()
{
    //-----User Settings:-----------------------------------------------
  TString  SimEngine      ="TGeant3"; 
  TString  InputFile     ="Mvd_Sim.root";
  //------------------------------------------------------------------


// Load basic libraries
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  gSystem->Load("libEve");
  gSystem->Load("libEventDisplay");
  gSystem->Load("libPndEventDisplay");

  PndFileNameCreator creator(InputFile.Data());
  TString digiFile = creator.GetDigiFileName();
  TString recoFile = creator.GetRecoFileName();
  TString trackF = creator.GetTrackFindingFileName();
  TString ParFile = creator.GetParFileName();
                                     
  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(InputFile.Data());
  fRun->AddFriend(recoFile.Data());
  fRun->AddFriend(digiFile.Data());
  fRun->AddFriend(trackF.Data());
  fRun->SetOutputFile("tst.root");
//  fRun->RunWithTimeStamps();
   FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(ParFile.Data());
  rtdb->setFirstInput(parInput1);
  FairEventManager *fMan= new FairEventManager();

  FairGeane *Geane = new FairGeane();
  fRun->AddTask(Geane);
//
 //     fRun->Init();
 
 // SetPalette(1);

 //----------------------Traks and points -------------------------------------
  FairMCTracks *Track =  new FairMCTracks ("Monte-Carlo Tracks");
  FairMCPointDraw *MvdPoints =   new FairMCPointDraw ("MVDPoint",kBlue,  kFullSquare);
//  FairMCPointDraw *EMCPoints =   new FairMCPointDraw ("EmcHit",kOrange,  kFullSquare);
//  FairMCPointDraw *TofPoint =    new FairMCPointDraw ("TofPoint",kYellow,  kFullSquare);
//  FairMCPointDraw *TofSciFPoint= new FairMCPointDraw ("TofSciFPoint",kTeal, kFullSquare);
//  FairMCPointDraw *MdtPoint =    new FairMCPointDraw ("MdtPoint",kAzure, kFullSquare);
//  FairMCPointDraw *PndDrcPoint = new FairMCPointDraw ("PndDrcPoint",kViolet, kFullSquare);
//  FairMCPointDraw *PndDchPoint = new FairMCPointDraw ("PndDchPoint",kPink, kFullSquare);
//  FairMCPointDraw *PndTpcPoint = new FairMCPointDraw ("PndTpcPoint",kCyan,  kFullSquare);
  FairMCPointDraw *PndSTTPoint = new FairMCPointDraw ("STTPoint",kMagenta, kFullSquare);
//  FairMCPointDraw *PndGEMPoint = new FairMCPointDraw ("GEMPoint",kRed, kFullSquare);
//  FairMCPointDraw *PndDskPoint = new FairMCPointDraw ("DskCerenkov",kGreen, kFullSquare);
    FairMCPointDraw *FtsPoint = new FairMCPointDraw("FTSPoint", kRed, kFullSquare);
//  FairHitDraw *EMCRecoHit = new FairHitDraw("EmcRecoHit");
//
  PndTrackCandDraw* RiemannCand = new PndTrackCandDraw("MVDRiemannTrackCand");
  RiemannCand->SetTimeWindowPlus(10);
  RiemannCand->SetTimeWindowMinus(10);
  PndTrackDraw* PndTrackRiemann = new PndTrackDraw("MVDTrack");
//  PndRiemannTrackDraw* RiemannTrack = new PndRiemannTrackDraw("MVDRiemannTrack");
//  PndMvdDigiPixelDraw* MvdDigiPixel = new PndMvdDigiPixelDraw("MVDPixelDigis");
                                                            
  FairHitDraw *MvdRecoHit =   new FairHitDraw ("MVDHitsPixel");
  MvdRecoHit->SetTimeWindowPlus(10);
  MvdRecoHit->SetTimeWindowMinus(10);
  FairHitDraw *MvdRecoStrip = new FairHitDraw ("MVDHitsStrip");
  MvdRecoStrip->SetTimeWindowPlus(10);
  MvdRecoStrip->SetTimeWindowMinus(10);
  FairHitDraw *STTHits = new FairHitDraw ("STTHit");
  STTHits->SetTimeWindowPlus(200);
  STTHits->SetTimeWindowMinus(1);
  fMan->AddTask(Track);
  fMan->AddTask(MvdPoints);
  
  PndSttIsochroneDraw* STTIsochrone = new PndSttIsochroneDraw("STTHit");
  STTIsochrone->SetTimeWindowPlus(300);
  STTIsochrone->SetTimeWindowMinus(10);
  STTIsochrone->UseIsochroneTime();
  fMan->AddTask(STTIsochrone);
//  fMan->AddTask(MvdDigiPixel);
//  fMan->AddTask(EMCPoints);
//  fMan->AddTask(TofPoint);
//  fMan->AddTask( TofSciFPoint);
//  fMan->AddTask( MdtPoint);
//  fMan->AddTask( PndDrcPoint);
//  fMan->AddTask( PndDchPoint);
//  fMan->AddTask( PndTpcPoint);
  fMan->AddTask( PndSTTPoint);
//  fMan->AddTask( PndGEMPoint);
//  fMan->AddTask( PndDskPoint);
  fMan->AddTask( FtsPoint);
//
//  fMan->AddTask(EMCRecoHit);
  fMan->AddTask(MvdRecoHit);
  fMan->AddTask(MvdRecoStrip);
  fMan->AddTask(STTHits);
  
//  fMan->AddTask(RiemannCand);
  fMan->AddTask(PndTrackRiemann);

  fMan->Init();

}

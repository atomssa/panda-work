eventDisplay()
{
  //-----User Settings:-----------------------------------------------
  TString  SimEngine      ="TGeant3"; 
  TString  InputFile     ="output/data_QGSP_BERT_EMV_MOM_0.5_pim.root";
  TString  ParFile       ="output/simpar_QGSP_BERT_EMV_MOM_0.5_pim.root";
  //------------------------------------------------------------------


  // Load basic libraries
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  gSystem->Load("libEve");
  gSystem->Load("libEventDisplay");

                                     
  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(InputFile.Data());
  fRun->SetOutputFile("tst.root");

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(ParFile.Data());
  rtdb->setFirstInput(parInput1);
  FairEventManager *fMan= new FairEventManager();
 
 
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
  FairHitDraw *EMCRecoHit = new FairHitDraw("EmcRecoHit");
                                                            
  fMan->AddTask(Track);
  fMan->AddTask(MvdPoints);
  fMan->AddTask(EMCPoints);   
  fMan->AddTask(TofPoint);   
  fMan->AddTask( TofSciFPoint);
  fMan->AddTask( MdtPoint);
  fMan->AddTask( PndDrcPoint);
  fMan->AddTask( PndDchPoint);
  fMan->AddTask( PndTpcPoint);
  fMan->AddTask( PndSTTPoint);
  fMan->AddTask( PndGEMPoint);
  fMan->AddTask( PndDskPoint);

  fMan->AddTask(EMCRecoHit);
  
  fMan->Init();                     

}

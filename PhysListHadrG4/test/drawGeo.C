

drawGeo()
{
  
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();
  
  TFile* file = new TFile("simparams.root");
  file->Get("FairBaseParSet"); 
  
  gGeoManager->SetVisLevel(3);
  gGeoManager->GetMasterVolume()->Draw("ogl");


}


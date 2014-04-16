{
  // Macro to plot the magnetic field
 
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
  basiclibs();

  // Load this example libraries
  gSystem->Load("libGeoBase");
  gSystem->Load("libParBase");
  gSystem->Load("libBase");
  gSystem->Load("libField");

  PndMultiField *fField= new PndMultiField();

  PndTransMap *map_t= new PndTransMap("TransMap", "R");
  PndDipoleMap *map_d1= new PndDipoleMap("DipoleMap1", "R");
  PndDipoleMap *map_d2= new PndDipoleMap("DipoleMap2", "R");
  PndSolenoidMap *map_s1= new PndSolenoidMap("SolenoidMap1", "R");
  PndSolenoidMap *map_s2= new PndSolenoidMap("SolenoidMap2", "R");
  PndSolenoidMap *map_s3= new PndSolenoidMap("SolenoidMap3", "R");
  PndSolenoidMap *map_s4= new PndSolenoidMap("SolenoidMap4", "R");

  fField->AddField(map_t);
  fField->AddField(map_d1);
  fField->AddField(map_d2);
  fField->AddField(map_s1);
  fField->AddField(map_s2);
  fField->AddField(map_s3);
  fField->AddField(map_s4);

  fField->Init();
     
  Double_t x=0;
  Double_t y=0;
  Double_t z=0;
  Double_t po[3], BB[3];
  Double_t Btot=0;
     
  TH1F *Bx=new TH1F("Bx(z)","Bx(z)",850,-200,650);
  TH1F *By=new TH1F("By(z)","By(z)",850,-200,650);
  TH1F *Bz=new TH1F("Bz(z)","Bz(z)",850,-200,650);
  TH1F *Btotal=new TH1F("B(z)","B(z)",850,-200,650);

  for (Int_t iz=0; iz<=850; iz++) 
    {
      z = -200 + Double_t(iz);
      po[0]=x; po[1]=y; po[2]=z;
      BB[0]=0; BB[1]=0; BB[2]=0;
      fField->GetFieldValue(po,BB); //return valuse in KG (G3)
      //  cout << "Z =" << z <<" Bx "<< BB[0] << " By  " << BB[1] <<" Bz  " << BB[2] <<endl; 
      Bx->SetBinContent(iz+1,BB[0]/ 10.);
      By->SetBinContent(iz+1,BB[1]/ 10.);
      Bz->SetBinContent(iz+1,BB[2]/ 10.);
      Btot=TMath::Sqrt(BB[0]*BB[0]+BB[1]*BB[1]+BB[2]*BB[2]);
      Btotal->SetBinContent(iz+1,Btot/ 10.);
    }
    
     
  THStack *h1= new THStack("B(z)","Panda Field X=Y=0") ;
     
  TCanvas *c1 = new TCanvas("c1", "c1",4,31,900,600);
  c1->Range(-267,-0.285949,403,2.38354);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
    

  Bx->SetLineColor(2);
  By->SetLineColor(3);
  Bz->SetLineColor(4);

  h1->Add(Bx);
  h1->Add(By);
  h1->Add(Bz);
  h1->Add(Btotal);

  h1->Draw("nostack");
   
  TPaveText *pt = new TPaveText(0.01,0.945,0.311609,0.995,"blNDC");
  pt->SetName("title");
  pt->SetBorderSize(2);
  pt->SetFillColor(19);
  TText *text = pt->AddText("Panda Field X=Y=0");
  pt->Draw();
   
  pt = new TPaveText(148.439,1.2000,250.534,1.94786,"br");
  pt->SetFillColor(19);
  pt->Draw();
  TLine *line = new TLine(160.,1.84535,193.685,1.84535);
  line->Draw();
   
  tex = new TLatex(200.,1.81766,"Bz");
  tex->SetLineWidth(2);
  line->SetLineColor(4);
  tex->Draw();

  line = new TLine(160.,1.6868,196.,1.6868);
  line->SetLineColor(2);
  line->Draw();

  tex = new TLatex(200.,1.64905,"Bx");
  tex->SetLineWidth(2);
  tex->Draw();

  line = new TLine(160.,1.51818,196.,1.51818);
  line->SetLineColor(3);
  line->Draw();

  tex = new TLatex(200.,1.4900,"By");
  tex->SetLineWidth(2);
  tex->Draw();


  line = new TLine(160.,1.31818,196.,1.31818);
  tex->SetLineWidth(2);
  line->Draw();

  tex = new TLatex(200.,1.3100,"Bmod");
  tex->SetLineWidth(2);
  tex->Draw();
   


  tex = new TLatex(-224.51,0.970725,"Tesla");
  tex->SetTextSize(0.054717);
  tex->SetTextAngle(90.0);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(37.0985,-0.200324,"Z (cm)");
  tex->SetLineWidth(2);
  tex->Draw();
 
  TCanvas *c2 = new TCanvas("c2", "c2",4,31,1200,700);
  c2->Divide(1,3);
  c2->cd(1);
  Bx->Draw();
  c2->cd(2);
  By->Draw();
  c2->cd(3);
  Bz->Draw();
  y=0;

  TH2F *B=new TH2F("B mod","B mod y=0 plane ",425,-200,650,250,-250,250);

  for (Int_t iz=0; iz<425; iz++) 
    {
      z = -200 + 2*Double_t(iz) ;
      po[1]=y; po[2]=z;
      for (Int_t ix=0; ix<250; ix++) {
	x= -250 + 2*Double_t(ix) ;
	//  cout << "X = " << x << endl; 
	po[0]=x;
	fField->GetFieldValue(po,BB); //return valuse in KG (G3)
	Btot=TMath::Sqrt(BB[0]*BB[0]+BB[1]*BB[1]+BB[2]*BB[2]);
	B->SetBinContent(iz,ix,Btot/ 10.);
      }
    }
 
  TCanvas *c3 = new TCanvas("c3", "c3",4,31,1200,700);
  B->Draw("cont1");
	
};













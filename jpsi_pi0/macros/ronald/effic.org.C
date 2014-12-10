void effic(Int_t NTmax=100)
// Read back the new proba files
// Test with gosia's tuples
{
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  rootlogon();

  Float_t procut=0.99999;
 
   FILE *fp;
//   Int_t NTmax=1000;
//   Int_t NFmax=10;   // max number of files
//   Int_t NTmax=100;   // max number per file
   Int_t NEVcount[10]={0,0,0,0,0, 0,0,0,0,0}; ;   // max number per file
//   Int_t NTmax=10000;
//   Int_t NTmax=200;

  Int_t kLog=1;
  Double_t p0=0.2, p1=5.0;
  if(kLog>0) {p0=-0.7; p1=0.7;}

  TRandom *momran       = new TRandom();
  Float_t radeg=57.29578;

  TString Directory[10]={
     "electronsP_trunk","muonsP_trunk","pionsP_trunk","kaonsP_trunk","protonsP_trunk",
     "electronsN_trunk","muonsN_trunk","pionsN_trunk","kaonsN_trunk","protonsN_trunk"
                        };
   
  Int_t Type[10]={0,0,0,0,0, 1,1,1,1,1}; 

  TString particle[10]= {"elecP","muonP","pionP","kaonP","protP",
                         "elecN","muonN","pionN","kaonN","protN"};


  TString detector[7]= {"STT","DIS","DRC","MUO","MVD","EMC","ALL"};
  TString type[4]= {"e-cut","e-all","pi-cut","pi-all"};

  TH2F *hist[7][5];
  TH1F *hpp[7][5];
  TH1F *hth[7][5];
  TH1F *hprobe[7][2];
  TH1F *hfactor[7][2];
  TH1F *hkcut[7][2][2];
  TCanvas *cEFF[7];
  TCanvas *cCUT[7];

  int off=32;
  int start=250;
  int nbins=60;

  TString aline;

  for (Int_t k=0; k<7 ; k++){
    aline= "cEFF"+detector[k];
    cEFF[k] = new TCanvas(aline,aline,start+k*off,(k+1)*off, 1000,1000);
    cEFF[k]->Divide(3,3);
//    aline= "cCUT"+detector[k];
//    cCUT[k] = new TCanvas(aline,aline,start+k*off,(k+1)*off, 1000,1000);
//    cCUT[k]->Divide(2,2);

    aline= "hprobe0"+detector[k];
    hprobe[k][0] = new TH1F(aline,aline,100,0,1);
    aline= "hprobe2"+detector[k];
    hprobe[k][1] = new TH1F(aline,aline,100,0,1);
    aline= "hkcut e=>e "+detector[k];
    hkcut[k][0][0] = new TH1F(aline,aline,100,-10,10);
    aline= "hkcut e=>pi "+detector[k];
    hkcut[k][1][0] = new TH1F(aline,aline,100,-10,10);
    aline= "hkcut pi=>e "+detector[k];
    hkcut[k][0][1] = new TH1F(aline,aline,100,-10,10);
    aline= "hkcut pi=>pi "+detector[k];
    hkcut[k][1][1] = new TH1F(aline,aline,100,-10,10);
    aline= "hfactor0"+detector[k];
    hfactor[k][0] = new TH1F(aline,aline,100,-10,10);
    aline= "hfactor2"+detector[k];
    hfactor[k][1] = new TH1F(aline,aline,100,-10,10);
     for (Int_t l=0; l<5 ; l++){
       aline= "hist"+detector[k]+type[l];
       hist[k][l] = new TH2F(aline,aline,nbins,0,5,nbins,0,180);
       hist[k][l]->Sumw2(); 
       aline= "hpp"+detector[k]+type[l];
       hpp[k][l] = new TH1F(aline,aline,nbins,0,5);
       hpp[k][l]->Sumw2(); 
       aline= "hth"+detector[k]+type[l];
       hth[k][l] = new TH1F(aline,aline,nbins,0,180);
       hth[k][l]->Sumw2(); 
     }
  }
  TCanvas *cCHECK = new TCanvas("cCHECK","check",start,off, 1200,1200);
  cCHECK->Divide(2,2);

//  TCanvas *cFACTOR = new TCanvas("cFACTOR","factor",start,off, 1200,800);
//  cFACTOR->Divide(3,2);

// canvas and stuff
      gStyle->SetLabelSize(0.05,"X");
      gStyle->SetLabelSize(0.05,"Y");
      gStyle->SetLineWidth(2);
      gStyle->SetHistLineWidth(2);
      gStyle->SetLabelSize(0.05,"X");
      gStyle->SetLabelSize(0.05,"Y");
      gStyle->SetPalette(1);
//      gStyle->SetOptFit(1);
      gStyle->SetOptStat(0);


cout << " finished histos " << endl;

// Get the probabilities for the detectors (all but EMC)
  aline="Newproba.root";
  if(kLog>0) aline="NewprobaLog.root";
  TFile* hfile1 = new TFile (aline); 

  TH2D  *hprob[10][5];
  for (Int_t ip=0; ip<10 ; ip++){
    for (Int_t k=0; k<5 ; k++){
       aline="hprob"+ particle[ip]+detector[k];
       hprob[ip][k] = (TH2D)hfile1->Get(aline);  
//       cout << ip << " " << k << " " << hprob[ip][k]->GetEntries();         
//       cout << aline->Data() << endl;         

    }
  }

// Loop over files
  Double_t probRK[5][7];
//  for (Int_t Did = 0; Did < 4; Did++) {
  for (Int_t Did = 0; Did < 3; Did++) {
     if(Did == 1) continue;
     TString inFile = Directory[Did]+".root";
     cout << "filename:" << inFile;

     TFile *hfile1 = TFile::Open(inFile,"READ");
     TNtuple *NTev;
     NTev = (TNtuple*) hfile1->Get("ntuple");
// process the data
     NTevents=NTev->GetEntriesFast();
     cout << " NTevents: " << NTevents << endl;
     if(NTevents>NTmax && NTmax>0) NTevents=NTmax;
     for (Int_t j=0; j< NTevents; j++) {
       NTev->GetEntry(j); 
       NEVcount[Did]++;   
       Float_t    momM         = NTev->GetArgs()[ 0]; // MC momentum
       Float_t    thetaM       = NTev->GetArgs()[ 1]; // MC theta
       Float_t    phiM         = NTev->GetArgs()[ 2]; // MC phi
       Float_t    recQ         = NTev->GetArgs()[ 3]; // Reco charge
       Float_t    MOMBase      = NTev->GetArgs()[ 4]; // Reco momentum
       if(MOMBase<0) continue;
       Float_t    thetaR       = NTev->GetArgs()[ 5]; // Reco theta
       if(thetaR<0) continue;
       Float_t    phiR         = NTev->GetArgs()[ 6]; // Reco phi
       Float_t    qR           = NTev->GetArgs()[ 7]; // Reco q ??
       Float_t    Z20          = NTev->GetArgs()[ 8]; // Zernik moments: Z20
       Float_t    Z53          = NTev->GetArgs()[ 9]; // Zernik moments: Z53
       Float_t    Er           = NTev->GetArgs()[10]; // Raw energy from EMC
       Float_t    Ec           = NTev->GetArgs()[11]; // Calibrated energy from EMC
       Float_t    lat          = NTev->GetArgs()[12]; // Lateral momenta from EMC
       Float_t    emc_qa       = NTev->GetArgs()[13]; // matching QA from EMC
       Float_t    emc_index    = NTev->GetArgs()[14]; // index to EMC
       Float_t    emc_crystal  = NTev->GetArgs()[15]; // number of cristal
       Float_t    e1           = NTev->GetArgs()[16]; // E1 from EMC
       Float_t    e9           = NTev->GetArgs()[17]; // E9 from EMC
       Float_t    e25          = NTev->GetArgs()[18]; // E25 from EMC
       Float_t    e1e9         = NTev->GetArgs()[19]; // ratio E1/E9 from EMC
       Float_t    e9e25        = NTev->GetArgs()[20]; // ratio E9/E25 from EMC
       Float_t    stt_dedx     = NTev->GetArgs()[21]; // truncated dE/dx from STT
       Float_t    stt_hits     = NTev->GetArgs()[22]; // number of hits in STT
       Float_t    mvd_dedx     = NTev->GetArgs()[23]; // dE/dx from MVD
       Float_t    mvd_hits     = NTev->GetArgs()[24]; // number of hits in MVD
       Float_t    muo_index    = NTev->GetArgs()[25]; // index to MUO (muon detector)
       Float_t    muo_nbLayer  = NTev->GetArgs()[26]; // number of crossed layers in MUO
       Float_t    muo_module   = NTev->GetArgs()[27]; // number of modules in MUO
       Float_t    muo_iron     = NTev->GetArgs()[28]; // amount of crossed iron in MUO
       Float_t    muo_qa       = NTev->GetArgs()[29]; // matching parameters from MUO
       Float_t    muo_momIn    = NTev->GetArgs()[30]; // momentum at entrance of MUO
       Float_t    drc_thetaC   = NTev->GetArgs()[31]; // thetaC from DRC
       Float_t    drc_qa       = NTev->GetArgs()[32]; // matching QA from DRC
       Float_t    drc_nbPh     = NTev->GetArgs()[33]; // number of photons from DRC
       Float_t    drc_index    = NTev->GetArgs()[34]; // index to DRC
       Float_t    disc_thetaC  = NTev->GetArgs()[35]; // thetaC from DIS
       Float_t    disc_qs      = NTev->GetArgs()[36]; // matching QA from DIS
       Float_t    disc_nbPh    = NTev->GetArgs()[37]; // number of photons from DIS
       Float_t    disc_index   = NTev->GetArgs()[38]; // index to DIS

//       Float_t    THBase = radeg*thetaM; // Reco theta
       Float_t    THBase = thetaR; // Reco theta

// emc on the fly
       probRK[0][5] = NTev->GetArgs()[69]; // ele/emc
       probRK[1][5] = NTev->GetArgs()[70]; 
       probRK[2][5] = NTev->GetArgs()[71]; 
       probRK[3][5] = NTev->GetArgs()[72]; 
       probRK[4][5] = NTev->GetArgs()[73]; 
       if(j<5) cout << " Did: " << Did ; 
       if(j<5) cout << " MOMBase: " << MOMBase ; 
       if(j<5) cout << " THBase: " << THBase ; 
       if(j<5) cout << " EP: " << Er/MOMBase << endl;
       if(j<5) cout << " probEMC: " << probRK[0][5]; 
       if(j<5) cout << " " << probRK[1][5]; 
       if(j<5) cout << " " << probRK[2][5]; 
       if(j<5) cout << " " << probRK[3][5]; 
       if(j<5) cout << " " << probRK[4][5] << endl; 


/*
       probGS[0][0] = NTev->GetArgs()[64]; // ele/stt
       probGS[1][0] = NTev->GetArgs()[65]; 
       probGS[2][0] = NTev->GetArgs()[66]; 
       probGS[3][0] = NTev->GetArgs()[67]; 
       probGS[4][0] = NTev->GetArgs()[68]; 

       probGS[0][2] = NTev->GetArgs()[74]; // ele/drc
       probGS[1][2] = NTev->GetArgs()[75]; 
       probGS[2][2] = NTev->GetArgs()[76]; 
       probGS[3][2] = NTev->GetArgs()[77]; 
       probGS[4][2] = NTev->GetArgs()[78]; 

       probGS[0][1] = NTev->GetArgs()[79]; // ele/disc
       probGS[1][1] = NTev->GetArgs()[80]; 
       probGS[2][1] = NTev->GetArgs()[81]; 
       probGS[3][1] = NTev->GetArgs()[82]; 
       probGS[4][1] = NTev->GetArgs()[83]; 
*/
//  Fill histos
       if(MOMBase>5) MOMBase=4.999999;
       if(MOMBase<0.2) MOMBase=0.200001;
       Float_t momR=MOMBase;
       if(kLog>0) momR=TMath::Log10(MOMBase);

       EMCBase = Er/MOMBase;
       STTBase = stt_dedx;
       DISBase = radeg*disc_thetaC;
       DRCBase = radeg*drc_thetaC;
       MUOBase = muo_iron;
       MVDBase = 1000*mvd_dedx;

       if(j<5) cout << " Did: " << Did ; 
       if(j<5) cout << " momR: " << MOMBase; 
       if(j<5) cout << " EMC: " << EMCBase; 
       if(j<5) cout << " STT: " << STTBase ; 
       if(j<5) cout << " DIS: " << DISBase ; 
       if(j<5) cout << " DRC: " << DRCBase ; 
       if(j<5) cout << " MUO: " << MUOBase ; 
       if(j<5) cout << " MVD: " << MVDBase << endl; 

// calculate the probabilities for each detector
// first index is particle, second is detector
// STT data
       Double_t DETbase[5], DETtrafo[5];
       DETbase[0] = STTBase;    DETtrafo[0]= -2;
       DETbase[1] = DISBase;    DETtrafo[1]= -2;
       DETbase[2] = DRCBase;    DETtrafo[2]= -2;
       DETbase[3] = MUOBase;    DETtrafo[3]= -2;
       DETbase[4] = MVDBase;    DETtrafo[4]= -2;
       Int_t binx    = (hprob[0][0]->GetXaxis())->FindBin(momR);
       Int_t biny;
       for (Int_t id=0; id<5 ; id++){
         for (Int_t ip=0; ip<5 ; ip++){
           probRK[ip][id]=0.2;
         }
       }
       for (Int_t id=0; id<5 ; id++){
         if(j<5) cout << " check base: " << id << " " << DETbase[id] << endl; 
         if(DETbase[id]>0.01) {
           if(id==0) DETtrafo[id]= TMath::ATan(2.00 *(DETbase[id]-7));
           if(id==1) DETtrafo[id]= TMath::ATan(1.25 *(DETbase[id]-44));
           if(id==2) DETtrafo[id]= TMath::ATan(1.25 *(DETbase[id]-44));
           if(id==3) DETtrafo[id]= TMath::ATan(0.125*(DETbase[id]-40));
           if(id==4) DETtrafo[id]= TMath::ATan(2.00 *(DETbase[id]-4));
           biny = (hprob[0][0]->GetYaxis())->FindBin(DETtrafo[id]);   
           if(j<5) cout << " check biny: " << biny << endl; 
           for (Int_t ip=0; ip<5 ; ip++){
//             probRK[ip][id]=hprob[ip][id]->GetBinContent(binx,biny);
             if(abs(DETtrafo[id]) < 1.57) {
                 probRK[ip][id]=hprob[ip][id]-> Interpolate(momR, DETtrafo[id]);
             }
             if(j<5 && ip==0 && id==0) {
                cout << " momR: " << momR ; 
                cout << " base: " << DETbase[id] ; 
                cout << " traf: " << DETtrafo[id] ; 
                cout << " binx: " << binx ; 
                cout << " biny: " << biny ; 
                cout << " prob: " << probRK[ip][id] << endl; 
             }
           }
         }
       }  //  for (Int_t id=0; id<5 ; id++)

// Fill ALL using EMC (5), STT (0), DIS(1), DRC(2)
      Double_t sum=0;
      for (Int_t ip=0; ip<5 ; ip++){
         Double_t Kfactor= probRK[ip][5]/(1-probRK[ip][5]);    // emc
         Kfactor *= probRK[ip][0]/(1-probRK[ip][0]);    // stt
         Kfactor *= probRK[ip][1]/(1-probRK[ip][1]);    // dis
         Kfactor *= probRK[ip][2]/(1-probRK[ip][2]);    // drc
         probRK[ip][6]=Kfactor/(1+Kfactor);
         sum += probRK[ip][6];
      }  //  for (Int_t ip=0; ip<5 ; ip++)
// Normalisation
      for (Int_t ip=0; ip<5 ; ip++){
         probRK[ip][6] /= sum;
      }  //  for (Int_t ip=0; ip<5 ; ip++)

//  probRK[ip][id] contains now the proba for particle ip, detector id 
//  TString detector[7]= {"STT","DIS","DRC","MUO","MVD","EMC","ALL"};
//  TString type[4]= {"e-cut","e-all","pi-cut","pi-all"};

//  TH2F *hist[7][4];
//  TH1F *hpp[7][4];
//  TH1F *hth[7][4];

       Double_t probe;
       for (Int_t id=0; id<7 ; id++){
//  e->e efficiency
         if(Did==0) {
           hist[id][1]->Fill(MOMBase,THBase);   // all for e
           hpp[id][1]->Fill(MOMBase);  
           hth[id][1]->Fill(THBase);  
           probe = probRK[0][id];          // proba for e
           hprobe[id][0]->Fill(probe);
           factor=TMath::Log10(probe/(1-probe));
           hfactor[id][0]->Fill(factor);
           if( (id==6 && probe>procut) ||  
             (id<6 && probe > probRK[1][id] && probe > probRK[2][id] && 
               probe > probRK[3][id] && probe > probRK[4][id])) { 
             hist[id][0]->Fill(MOMBase,THBase);   // e best
             hpp[id][0]->Fill(MOMBase);  
             hth[id][0]->Fill(THBase);  
           }
         }
       
//  pi->e mis id
         if(Did==2) {
           hist[id][3]->Fill(MOMBase,THBase);   // all for e
           hpp[id][3]->Fill(MOMBase);  
           hth[id][3]->Fill(THBase);  
           probe = probRK[0][id];          // proba for e
           hprobe[id][1]->Fill(probe);
           factor=TMath::Log10(probe/(1-probe));
           hfactor[id][1]->Fill(factor);
           if( (id==6 && probe>procut) ||  
             (id<6 && probe > probRK[1][id] && probe > probRK[2][id] && 
               probe > probRK[3][id] && probe > probRK[4][id])) { 
             hist[id][2]->Fill(MOMBase,THBase);   // e best
             hpp[id][2]->Fill(MOMBase);  
             hth[id][2]->Fill(THBase);  
           }
         }

       }    // for (Int_t id=0; id<7 ; id++)

       Double_t probe = probRK[0][6];          // proba for e
       if(probe>procut) { 
         for (Int_t id=0; id<7 ; id++) {
/*
           cout << "Did: " << Did << " id: " << id;
           cout << " " << probRK[0][id] ;
           cout << " " << probRK[1][id] ;
           cout << " " << probRK[2][id] ;
           cout << " " << probRK[3][id] ;
           cout << " " << probRK[4][id] << endl;
*/
           if(Did==0) {
             hkcut[id][0][0]->Fill(TMath::Log10(probRK[0][id]));
             hkcut[id][1][0]->Fill(TMath::Log10(probRK[2][id]));
           }
           if(Did==2) {
             hkcut[id][0][1]->Fill(TMath::Log10(probRK[0][id]));
             hkcut[id][1][1]->Fill(TMath::Log10(probRK[2][id]));
           }
         }
       }       
       
     }   // for (Int_t j=0; j< NTevents; j++) 
   } //for (Int_t Did = 0; Did < 10; Did++) 

   cout << "finished files" << endl;

// plots

//    aline= "cEFF"+detector[k];
//    cEFF[k] = new TCanvas(aline,aline,start+k*off,(k+1)*off, 1200,900);
//    cEFF[k]->Divide(3,2);

// loop over detectors
       for (Int_t id=6; id<7 ; id++) {
         cEFF[id]->cd(1); hist[id][0]->Divide(hist[id][1]);
             aline="efficiency "+detector[id];hist[id][0]->SetTitle(aline);
             hist[id][0]->SetMaximum(1);hist[id][0]->SetMinimum(0.00001);
             hist[id][0]->Draw("COLZ");
         cEFF[id]->cd(4); hist[id][2]->Divide(hist[id][3]);
             aline="misid "+detector[id]; hist[id][2]->SetTitle(aline);
             hist[id][2]->SetMaximum(1);hist[id][2]->SetMinimum(0.00001);
             hist[id][2]->Draw("COLZ");
         cEFF[id]->cd(7);hist[id][4]->Divide(hist[id][0],hist[id][2],1,1,"");
             aline="S/B "+detector[id]; hist[id][4]->SetTitle(aline);
             hist[id][4]->Draw("COLZ");

         cEFF[id]->cd(2); hpp[id][0]->Divide(hpp[id][1]);
             aline="efficiency "+detector[id];hpp[id][0]->SetTitle(aline);
             hpp[id][0]->SetMaximum(1);hpp[id][0]->SetMinimum(0.00001);
             hpp[id][0]->Draw();
         cEFF[id]->cd(5); hpp[id][2]->Divide(hpp[id][3]);
             aline="misid "+detector[id]; hpp[id][2]->SetTitle(aline);
             hpp[id][2]->SetMaximum(1);hpp[id][2]->SetMinimum(0.00001);
             hpp[id][2]->Draw();
         cEFF[id]->cd(8);hpp[id][4]->Divide(hpp[id][0],hpp[id][2],1,1,"");
             aline="S/B "+detector[id]; hpp[id][4]->SetTitle(aline);
             hpp[id][4]->Draw();

         cEFF[id]->cd(3); hth[id][0]->Divide(hth[id][1]);
             aline="efficiency "+detector[id];hth[id][0]->SetTitle(aline);
             hth[id][0]->SetMaximum(1);hth[id][0]->SetMinimum(0.00001);
             hth[id][0]->Draw("COLZ");
         cEFF[id]->cd(6); hth[id][2]->Divide(hth[id][3]);
             aline="misid "+detector[id]; hth[id][2]->SetTitle(aline);
             hth[id][2]->SetMaximum(1);hth[id][2]->SetMinimum(0.00001);
             hth[id][2]->Draw("COLZ");
         cEFF[id]->cd(9);hth[id][4]->Divide(hth[id][0],hth[id][2],1,1,"");
             aline="S/B "+detector[id]; hth[id][4]->SetTitle(aline);
             hth[id][4]->Draw();


/*
         cCUT[id]->cd(1); hkcut[id][0][0]->Draw();
         cCUT[id]->cd(2); hkcut[id][1][0]->Draw();
         cCUT[id]->cd(3); hkcut[id][0][1]->Draw();
         cCUT[id]->cd(4); hkcut[id][1][1]->Draw();
*/
   cout << "finished cCUT: " << id << endl;

       }    // for (Int_t id=0; id<7 ; id++)


/*
   cCHECK->cd(1);hprobe[0][0]->Draw();hprobe[0][1]->SetLineColor(kRed);hprobe[0][1]->Draw("same");
   cCHECK->cd(2);hprobe[1][0]->Draw();hprobe[1][1]->SetLineColor(kRed);hprobe[1][1]->Draw("same");
   cCHECK->cd(3);hprobe[2][0]->Draw();hprobe[2][1]->SetLineColor(kRed);hprobe[2][1]->Draw("same");
   cCHECK->cd(4);hprobe[4][0]->Draw();hprobe[4][1]->SetLineColor(kRed);hprobe[4][1]->Draw("same");
   cCHECK->cd(5);hprobe[5][0]->Draw();hprobe[5][1]->SetLineColor(kRed);hprobe[5][1]->Draw("same");
   cCHECK->cd(6);hprobe[6][0]->Draw();hprobe[6][1]->SetLineColor(kRed);hprobe[6][1]->Draw("same");
   cCHECK->cd(0);

   cFACTOR->cd(1);hfactor[0][0]->Draw();hfactor[0][1]->SetLineColor(kRed);hfactor[0][1]->Draw("same");
   cFACTOR->cd(2);hfactor[1][0]->Draw();hfactor[1][1]->SetLineColor(kRed);hfactor[1][1]->Draw("same");
   cFACTOR->cd(3);hfactor[2][0]->Draw();hfactor[2][1]->SetLineColor(kRed);hfactor[2][1]->Draw("same");
   cFACTOR->cd(4);hfactor[4][0]->Draw();hfactor[4][1]->SetLineColor(kRed);hfactor[4][1]->Draw("same");
   cFACTOR->cd(5);hfactor[5][0]->Draw();hfactor[5][1]->SetLineColor(kRed);hfactor[5][1]->Draw("same");
   cFACTOR->cd(6);hfactor[6][0]->Draw();hfactor[6][1]->SetLineColor(kRed);hfactor[6][1]->Draw("same");
   cFACTOR->cd(0);
*/
// save some histos

   cout << "before heffelec" << endl;
   TH2F *heffelec = (TH2F*)hist[6][0]->Clone("heffelec");
   TH2F *heffpion = (TH2F*)hist[6][2]->Clone("heffpion");
   TH1F *hpppion = (TH1F*)hpp[6][2]->Clone("hpppion");
   TH1F *hthpion = (TH1F*)hth[6][2]->Clone("hthpion");
   cout << "after heffelec" << endl;
   cCHECK->cd(1); heffelec -> Draw("COLZ");
   cCHECK->cd(2); heffpion -> Draw("COLZ");
   cCHECK->cd(3); hpppion -> Draw("");
   cCHECK->cd(4); hthpion -> Draw("");


   TString outFile = "effic.root";
   cout << "filename:" << outFile << endl;
   TFile *out = TFile::Open(outFile,"RECREATE");

   out->cd();

   heffelec -> Write();
   heffpion -> Write();
   hpppion -> Write();
   hthpion -> Write();
   out->Save();
   out->Close("R");

   cout << "Yahoo!" << endl;


}



class TCandList;
class TCandidate;
class TFitParams;

void newana_check_eta(int nevts=0)
{
  TString OutFile="output.root";
  TString inPidFile  = "evt_pid_stt.root";
  TString inRecoFile  = "evt_reco_stt.root";
  TString inSimFile = "evt_points_stt.root";
  TString inParFile = "evt_params_stt.root";
  
  gStyle->SetOptFit(1011);
  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");
 
  FairLogger::GetLogger()->SetLogToFile(kFALSE);
  //FairLogger::GetLogger()->SetLogScreenLevel("DEBUG4");
  //FairLogger::GetLogger()->SetLogVerbosityLevel("MEDIUM");
  
  FairRunAna* fRun = new FairRunAna();
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  fRun->SetInputFile(inSimFile);
  fRun->AddFriend(inPidFile);
  fRun->AddFriend(inRecoFile);
  
  FairParRootFileIo* parIO = new FairParRootFileIo();
  parIO->open(inParFile);
  rtdb->setFirstInput(parIO);
  rtdb->setOutput(parIO);  
  
  fRun->SetOutputFile(OutFile);
  fRun->Init(); 

  TFile *out = TFile::Open("output_newmc.root","RECREATE");


  TH1F *h_etac_nocut=new TH1F("h_etac_nocut","m(eta_c), (no cuts);E, GeV",100,2.5,3.5);
  TH1F *h_mphi_nocuts=new TH1F("h_mphi_nocuts","#phi: m(K+ K-) (no cuts)",100,0.95,1.1);

  TH1F *h_etac_pid=new TH1F("h_etac_pid","m(eta_c), (MC PID);E, GeV",100,2.5,3.5);
  TH1F *h_mphi_pid=new TH1F("h_mphi_pid","#phi: m(K+ K-) (MC PID)",100,0.95,1.1);
  
  TH1F *h_etac_4c=new TH1F("h_etac_4c","m(eta_c), 4C-fit",100,2.5,3.5);
  TH1F *h_mphi_4c=new TH1F("h_mphi_4c","#phi: m(K+ K-) (4C-fit)",100,0.95,1.1);
  
  TH1F *h_etac_vtx=new TH1F("h_etac_vtx","m(eta_c), Vertex fit",100,2.5,3.5);
  TH1F *h_mphi_vtx=new TH1F("h_mphi_vtx","#phi: m(K+ K-) (Vertex fit)",100,0.95,1.1);
        
  TH1F *h_etac_phimass=new TH1F("h_etac_phimass","m(eta_c), (cut on #phi mass);E, GeV",100,2.5,3.5);
  TH1F *h_mphi_final=new TH1F("h_mphi_final","#phi: m(K+ K-)",100,0.95,1.1);

  TH1F *nc=new TH1F("nc","n charged",20,0,20);

  TH1F *h_chi2_4c=new TH1F("h_chi2_4c","#chi^{2} 4C-fit;#chi^{2}/N_{df}",100,0,100);
  TH1F *h_chi2b_4c=new TH1F("h_chi2b_4c","#chi^{2} 4C-fit;#chi^{2}/N_{df}",100,0,100);
  TH1F *h_chi2_vtx=new TH1F("h_chi2_vtx","#chi^{2} vertex;#chi^{2}/N_{df}",100,0,100);
  TH2F *hvpos = new TH2F("hvpos","(x,y) projection of fitted decay vertex",100,-5,5,100,-5,5);
  TH1F *hvzpos = new TH1F("hvzpos","z position of fitted decay vertex",100,-10,10);
  TH1F *hvtxresX = new TH1F("hvtxresX","X resolution of fitted decay vertex",100,-0.1,0.1);
  TH1F *hvtxresY = new TH1F("hvtxresY","Y resolution of fitted decay vertex",100,-0.1,0.1);
  TH1F *hvtxresZ = new TH1F("hvtxresZ","Z resolution of fitted decay vertex",100,-0.1,0.1);
  
  TPidMassSelector *phiMassSel=new TPidMassSelector("phi",1.02,0.02);

  // the candidates lists we need
  TCandList p1, p2, phi1, phi1_pid, etac, etac_pid,  etac_nocut, mctracks;
  
  int n_reco=0;
  // Number of events in file and number of reconstructed eta_c to store in root file
  TH1F *n_events=new TH1F("n_events","total number of events",1,0,1);
  TH1F *n_etac=new TH1F("n_etac","number of reconstructed eta_c",1,0,1);

  TLorentzVector ini(0,0,3.6772,4.7333);
  
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts= theAnalysis->GetEntries();
  n_events->SetBinContent(1,nevts);
  // cout << "nevts " << nevts << "\n";
  int i=0,j=0, k=0, l=0;

  // *************
  // this is the loop through the events ... as simple as this...
  // ****************
  while (theAnalysis->GetEvent() && i++<nevts)
    {
      // cout << "evt: " << i << endl;
      
      if ((i%100)==0)  cout<<"evt " << i << endl;
      //if (!((i+1)%100)) cout<<"evt " << i << "\n";
      
      theAnalysis->FillList(mctracks,"McTruth"); 
      theAnalysis->FillList(p1,"KaonVeryLoosePlus");
      theAnalysis->FillList(p2,"KaonVeryLooseMinus");
      
      int nchrg=p1.GetLength()+p2.GetLength();
      nc->Fill(nchrg);
      
      phi1.Combine(p1,p2);
      for (j=0;j<phi1.GetLength();++j) h_mphi_nocuts->Fill(phi1[j].M());
      
      phi1.Select(phiMassSel);
      etac_nocut.Combine(phi1,phi1);
      
      for (l=0;l<etac_nocut.GetLength();++l) {
       	h_etac_nocut->Fill(etac_nocut[l].M());
      }
    
      FairMCEventHeader*  evthead = theAnalysis->GetEventHeader();
      TVector3 mcVertex, mcD1Vertex, mcD2vertex;
      evthead->GetVertex(mcVertex);
      
      //MC PID
      // Leave only kaons in particle lists
      int n_removed=0;
      int ii=0;
      for (l=0;l<p1.GetLength();++l) {
       	ii=l-n_removed;
	Int_t mcIndex = p1[ii].GetMicroCandidate().GetMcIndex();
	if (mcIndex>-1){
	  if ((mctracks[mcIndex].PdgCode()!=321))// || (mctracks[mcIndex].GetMotherID()!=-1))
	    {
	      p1.Remove(p1[ii]);
	      n_removed++;
	    }
	  if (mcIndex==0) mcD1Vertex = mctracks[mcIndex].Pos();
	  if (mcIndex==2) mcD1Vertex = mctracks[mcIndex].Pos();
	}
	else
	  {
	    std::cout<<"stt h: " << p1[ii].GetMicroCandidate().GetSttHits() << std::endl;
	    std::cout<<"Kaon list 1, element "<<l<<" has no assosiated mcTRack"<<std::endl;
	  }
      }
      
      n_removed = 0;
      ii = 0;
      for (l=0;l<p2.GetLength();++l) {
       	ii=l-n_removed;
	Int_t mcIndex = p2[ii].GetMicroCandidate().GetMcIndex();
	if (mcIndex>-1){
	  if ((mctracks[mcIndex].PdgCode()!=-321))// || (mctracks[mcIndex].GetMotherID()!=-1))
	    {
	      p2.Remove(p2[ii]);
	      n_removed++;
	    }
	  if (mcIndex==1) mcD1Vertex = mctracks[mcIndex].Pos();
	  if (mcIndex==3) mcD1Vertex = mctracks[mcIndex].Pos();
	}
	else
	  {
	    std::cout<<"stt h: " << p2[ii].GetMicroCandidate().GetSttHits() << std::endl;
	    std::cout<<"Kaon list 2, element "<<l<<" has no assosiated mcTRack"<<std::endl;
	  }
      }
      
      phi1_pid.Combine(p1,p2);
      
      for (j=0;j<phi1_pid.GetLength();++j) h_mphi_pid->Fill(phi1_pid[j].M());
      phi1_pid.Select(phiMassSel);
      etac.Combine(phi1_pid,phi1_pid); 
      
      for (l=0;l<etac.GetLength();++l) {
       	h_etac_pid->Fill(etac[l].M());
      }
      //cout << " ";// << endl;
      //if (use4cfit)
      	{
      	  ////////////// 4C-fit ///////////////
      	  int best_i=0;
      	  double best_chi2=1000;
      	  TCandidate *ccfit;// = new TCandidate();
      	  double m_phi1, m_phi2;
      	  for (l=0;l<etac.GetLength();++l) {
      	    Pnd4CFitter fitter(etac[l],ini);
      	    fitter.FitConserveMasses();
      	    double chi2=fitter.GetChi2();
      	    if (chi2<best_chi2)
      	      {
      		best_chi2=chi2;
      		best_i = l;
      		ccfit = (TCandidate*)fitter.FittedCand(etac[l]);
      		TCandidate *phi1best = fitter.FittedCand(*(etac[l].Daughter(0)));
      		TCandidate *phi2best = fitter.FittedCand(*(etac[l].Daughter(1)));
      		TCandidate *k1best = fitter.FittedCand(*(etac[l].Daughter(0)->Daughter(0)));
      		TCandidate *k2best = fitter.FittedCand(*(etac[l].Daughter(0)->Daughter(1)));
      		TCandidate *k3best = fitter.FittedCand(*(etac[l].Daughter(1)->Daughter(0)));
      		TCandidate *k4best = fitter.FittedCand(*(etac[l].Daughter(1)->Daughter(1)));
      		TLorentzVector tlvk1 = k1best->P4();
      		TLorentzVector tlvk2 = k2best->P4();
      		TLorentzVector tlvk3 = k3best->P4();
      		TLorentzVector tlvk4 = k4best->P4();
      		m_phi1= (tlvk1+tlvk2).M();
      		m_phi2= (tlvk3+tlvk4).M();
      	      }
      	    h_chi2_4c->Fill(chi2/9); // Ndf=3N-3=9
      	  }
	  
      	  if (etac.GetLength()!=0) cout << "evt: " << i << "\tbest: " << best_chi2 << "\tlen " << etac.GetLength() << endl;  
	 
      	  if(/*(best_chi2<270)&&*/(etac.GetLength()!=0))
      	    {
      	      h_chi2b_4c->Fill(best_chi2/9); // Ndf=3N-3=9
      	      h_etac_4c->Fill(ccfit->M());
      	      h_mphi_4c->Fill(m_phi1);
      	      h_mphi_4c->Fill(m_phi2);
      	      h_mphi_final->Fill(m_phi1);
      	      h_mphi_final->Fill(m_phi2);
      	      if (((m_phi1>1.02-0.02)&&(m_phi1<1.02+0.02))&&((m_phi2>1.02-0.02)&&(m_phi2<1.02+0.02)))
      	      	{
      	      	  h_etac_phimass->Fill(etac[best_i].M());
      	      	  if ((etac[best_i].M()>2.9)&&(etac[best_i].M()<3.06))
      	      	    n_reco++;
      	      	}
      	    }
      	}
      	// else // use vertex fit
      	{ 
      	  TCandList etac_vtx;
      	  TCandidate *k1, *k2, *k3, *k4, *phi1tmp, *phi2tmp, *etac_tmp;
	  
      	  //Combine 4 kaons directly to candidates
      	  for (j=0;j<etac.GetLength();++j)
      	    {
      	      phi1tmp=etac[j].Daughter(0);
      	      phi2tmp=etac[j].Daughter(1);

      	      k1=phi1tmp->Daughter(0);
      	      k2=phi1tmp->Daughter(1);
      	      k3=phi2tmp->Daughter(0);
      	      k4=phi2tmp->Daughter(1);
			
      	      etac_tmp=k1->Combine(*k2,*k3,*k4);
      	      etac_vtx.Add(*etac_tmp);
      	    }

      	  for (l=0;l<etac.GetLength();++l) {
      	    h_etac_pid->Fill(etac_vtx[l].M());
      	  }

      	  ////////////// Vertex fit /////////////
      	  int best_i=0;
      	  double best_chi2=1000;
      	  TCandidate *etacfit_best=0;
      	  TCandidate *k1fit_best, *k2fit_best, *k3fit_best, *k4fit_best;
      	  TCandidate *phi1fit_best, *phi2fit_best;
      	  TVector3 bestPos;
      	  Float_t etacvtx_mass;
      	  for (j=0;j<etac_vtx.GetLength();++j)
      	    {
      	      PndKinVtxFitter vtxfitter(etac_vtx[j]);        // instantiate a vertex fitter
	      vtxfitter.AddMassConstraint(ini.M());
      	      vtxfitter.Fit();                          // do the vertex fit

      	      TCandidate *etacfit=vtxfitter.FittedCand(etac_vtx[j]);  // request the fitted EtaC candidate
      	      TVector3 etacVtx=etacfit->Pos();                    // and the decay vertex position
      	      double chi2_vtx=vtxfitter.GlobalChi2();
      	      h_chi2_vtx->Fill(chi2_vtx/5); // Number degree of freedom 2N-3=5
      	      // plot mass and vtx x,y projection after fit
      	      hvpos->Fill(etacVtx.X(),etacVtx.Y());
      	      hvzpos->Fill(etacVtx.Z());
      	      if(chi2_vtx<best_chi2)
      		{
      		  best_chi2=chi2;
      		  best_i=l;
      		  etacfit_best=etacfit;
      		  k1fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(0)));
      		  k2fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(1)));
      		  k3fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(2)));
      		  k4fit_best=vtxfitter.FittedCand(*(etacfit_best->Daughter(3)));
      		  etacvtx_mass = etacfit_best->M();
      		  bestPos = etacfit->Pos(); 
      		}

      	    }
      	  if(/*(best_chi2<150)&&*/(etac.GetLength()!=0))
      	    {
      	      //h_etac_vtx->Fill(etacfit_best->M());
      	      h_etac_vtx->Fill(etacvtx_mass);
      	      phi1fit_best=k1fit_best->Combine(*k2fit_best);
      	      phi2fit_best=k3fit_best->Combine(*k4fit_best);
      	      double m_phi1=phi1fit_best->M();
      	      double m_phi2=phi2fit_best->M();
      	      h_mphi_vtx->Fill(m_phi1);
      	      h_mphi_vtx->Fill(m_phi2);
      	      h_mphi_final->Fill(m_phi1);
      	      h_mphi_final->Fill(m_phi2);
      	      hvtxresX->Fill(mcVertex.X()-bestPos.X());
      	      hvtxresY->Fill(mcVertex.Y()-bestPos.Y());
      	      hvtxresZ->Fill(mcVertex.Z()-bestPos.Z());
 
	      if (((m_phi1>1.02-0.03)&&(m_phi1<1.02+0.03))&&((m_phi2>1.02-0.03)&&(m_phi2<1.02+0.03)))
      	       	{
		  h_etac_phimass->Fill(etacvtx_mass);
		  if ((etacvtx_mass>2.9)&&(etacvtx_mass<3.06))
      	       	    n_reco++;
      	       	}
      	    }
	  
      	} 
    }
  std::cout<<"Number of reconstructed eta_c = "<<n_reco<<std::endl;
  n_etac->SetBinContent(1,n_reco);


  out->cd();
  n_etac->Write();
  n_events->Write();
  h_etac_nocut->Write();
  h_etac_pid->Write();
  h_etac_phimass->Write();
  h_etac_vtx->Write();
  h_etac_4c->Write();

  h_mphi_nocuts->Write();
  h_mphi_pid->Write();
  h_mphi_vtx->Write();
  h_mphi_4c->Write();
  h_mphi_final->Write();

  nc->Write();

  h_chi2_4c->Write();
  h_chi2_vtx->Write();
  hvzpos->Write();
  hvpos->Write();
  hvtxresX->Write();
  hvtxresY->Write();
  hvtxresZ->Write();
  out->Save();

}

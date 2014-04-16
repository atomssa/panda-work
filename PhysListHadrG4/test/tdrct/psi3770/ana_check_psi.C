class TCandList;
class TCandidate;
class TFitParams;

void ana_check(int nevts=0)
{
  TString OutFile = "output.root";
  
  gStyle->SetOptFit(1011);

  TStopwatch timer;
  timer.Start();

  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");

  TString inPidFile  = "evt_pid.root";
  TString inSimFile = "evt_points.root";

  TFile *inFile = TFile::Open(inSimFile,"READ");
  TTree *tree=(TTree *) inFile->Get("cbmsim") ;
  tree->AddFriend("cbmsim",inPidFile);

  TClonesArray* mc_array=new TClonesArray("PndMCTrack");
  tree->SetBranchAddress("MCTrack",&mc_array);

  FairMCEventHeader* evthead;
  tree->SetBranchAddress("MCEventHeader.", &evthead);
  
  TFile *out = TFile::Open(OutFile,"RECREATE");

  // the PndEventReader takes care about file/event handling
  PndEventReader evr(inPidFile);

  TH1F *h_psi_nocut=new TH1F("h_psi_nocut","m(psi), (no cuts);E, GeV",100,3.2,4.4);
  TH1F *h_mD1_nocuts=new TH1F("h_mD1_nocuts","D+: m(#pi+ #pi+ K-) (no cuts)",100,1.5,2.2);
  TH1F *h_mD2_nocuts=new TH1F("h_mD2_nocuts","D-: m(#pi- #pi- K+) (no cuts)",100,1.5,2.2);
  TH1F *h_mD1_pid=new TH1F("h_mD1_pid","D+: m(#pi+ #pi+ K-) (pid)",100,1.5,2.2);
  TH1F *h_mD2_pid=new TH1F("h_mD2_pid","D-: m(#pi- #pi- K+) (pid)",100,1.5,2.2);

  TH1F *h_psi_pid=new TH1F("h_psi_pid","m(#psi), (MC PID);E, GeV",100,3.2,4.4);
  TH1F *h_psi_pid_vtx=new TH1F("h_psi_pid_vtx","m(#psi), (MC PID);E, GeV",100,3.2,4.4);
  
  TH1F *h_psi_4c=new TH1F("h_psi_4c","m(#psi), 4C-fit",100,3.2,4.4);
  TH1F *h_mD_4c=new TH1F("h_mD_4c","D: m(K #pi #pi) (4C-fit)",100,1.5,2.2);
  
  TH1F *h_psi_vtx=new TH1F("h_psi_vtx","m(#psi), Vertex fit",100,3.2,4.4);
  TH1F *h_D1_vtx_bef=new TH1F("h_D1_vtx_bef","#phi: m(K #pi #pi) (Vertex fit)",100,1.5,2.2);  
  TH1F *h_D2_vtx_bef=new TH1F("h_D2_vtx_bef","#phi: m(K #pi #pi) (Vertex fit)",100,1.5,2.2);  
  TH1F *h_D1_vtx=new TH1F("h_D1_vtx","#phi: m(K #pi #pi) (Vertex fit)",100,1.5,2.2);  
  TH1F *h_D2_vtx=new TH1F("h_D2_vtx","#phi: m(K #pi #pi) (Vertex fit)",100,1.5,2.2);
        
  TH1F *h_psi_Dmass=new TH1F("h_psi_Dmass","m(#psi), (cut on D mass);E, GeV",100,3.2,4.4);
  TH1F *h_mD_final=new TH1F("h_mD_final","D: m(K #pi #pi)",100,1.5,2.2);

  TH1F *nc=new TH1F("nc","n charged",20,0,20);

  TH1F *h_chi2_4c=new TH1F("h_chi2_4c","#chi^{2} 4C-fit;#chi^{2}/N_{df}",100,0,100);
  TH1F *h_chi2b_4c=new TH1F("h_chi2b_4c","#chi^{2} 4C-fit;#chi^{2}/N_{df}",100,0,100);
  TH1F *h_chi2_vtx_D1=new TH1F("h_chi2_vtx_D1","#chi^{2} vertex;#chi^{2}/N_{df}",100,0,100); 
  TH1F *h_chi2_vtx_D2=new TH1F("h_chi2_vtx_D2","#chi^{2} vertex;#chi^{2}/N_{df}",100,0,100);
  TH2F *hvpos_D1 = new TH2F("hvpos_D1","(x,y) projection of fitted decay vertex",100,-5,5,100,-5,5);
  TH1F *hvzpos_D1 = new TH1F("hvzpos_D1","z position of fitted decay vertex",100,-10,10);
  TH1F *hvtxresX_D1 = new TH1F("hvtxresX_D1","X resolution of fitted decay vertex",100,-0.1,0.1);
  TH1F *hvtxresY_D1 = new TH1F("hvtxresY_D1","Y resolution of fitted decay vertex",100,-0.1,0.1);
  TH1F *hvtxresZ_D1 = new TH1F("hvtxresZ_D1","Z resolution of fitted decay vertex",100,-0.1,0.1);
  TH2F *hvpos_D2 = new TH2F("hvpos_D2","(x,y) projection of fitted decay vertex",100,-5,5,100,-5,5);
  TH1F *hvzpos_D2 = new TH1F("hvzpos_D2","z position of fitted decay vertex",100,-10,10);
  TH1F *hvtxresX_D2 = new TH1F("hvtxresX_D2","X resolution of fitted decay vertex",100,-0.1,0.1);
  TH1F *hvtxresY_D2 = new TH1F("hvtxresY_D2","Y resolution of fitted decay vertex",100,-0.1,0.1);
  TH1F *hvtxresZ_D2 = new TH1F("hvtxresZ_D2","Z resolution of fitted decay vertex",100,-0.1,0.1);
  
  TPidMassSelector *psiMassSel=new TPidMassSelector("psi",1.87,0.07);

  TPidPlusSelector *pplusSel=new TPidPlusSelector("pplus");
  TPidMinusSelector *pminusSel=new TPidMinusSelector("pminus");

  // the candidates lists we need
  TCandList p1, p2, k1, k2, D1, D1_pid, D2, D2_pid, psi, psi_pid,  psi_nocut;
  
  int n_reco=0;
  // Number of events in file and number of reconstructed eta_c to store in root file
  TH1F *n_events=new TH1F("n_events","total number of events",1,0,1);
  TH1F *n_psi=new TH1F("n_psi","number of reconstructed eta_c",1,0,1);

  TLorentzVector ini(0,0,6.578800,7.583333);

  if (nevts==0) nevts=evr.GetEntries();

  n_events->SetBinContent(1,nevts);
  // cout << "nevts " << nevts << "\n";
  int i=0,j=0, k=0, l=0;

  // *************
  // this is the loop through the events ... as simple as this...
  // ****************
  while (evr.GetEvent() && i++<nevts)
    {
      //.cout << "evt: " << i << endl;
      if ((i%100)==0)  cout<<"evt " << i << "\n";
      //if (!((i+1)%100)) cout<<"evt " << i << "\n";
      evr.FillList(p1,"Charged");
      evr.FillList(p2,"Charged");
      evr.FillList(k1,"Charged");
      evr.FillList(k2,"Charged");

      p1.Select(pplusSel);
      p2.Select(pminusSel);
      k1.Select(pplusSel);
      k2.Select(pminusSel);
      
      int nchrg=p1.GetLength()+p2.GetLength();
      nc->Fill(nchrg);

      for (j=0;j<p1.GetLength();++j) {
	p1[j].SetMass(TRho::Instance()->GetPDG()->GetParticle(211)->Mass());
      }
      for (j=0;j<p2.GetLength();++j) {
	p2[j].SetMass(TRho::Instance()->GetPDG()->GetParticle(211)->Mass());
      }
      for (j=0;j<k1.GetLength();++j) {
	k1[j].SetMass(TRho::Instance()->GetPDG()->GetParticle(321)->Mass());
      }
      for (j=0;j<k2.GetLength();++j) {
	k2[j].SetMass(TRho::Instance()->GetPDG()->GetParticle(321)->Mass());
      }
      
      
      D1.Combine(p1,p1, k2);
      D2.Combine(p2,p2, k1);
      
      for (j=0;j<D1.GetLength();++j) h_mD1_nocuts->Fill(D1[j].M());
      for (j=0;j<D2.GetLength();++j) h_mD2_nocuts->Fill(D2[j].M());
      
      D1.Select(psiMassSel);
      D2.Select(psiMassSel);
      psi_nocut.Combine(D1,D2);
      
      for (l=0;l<psi_nocut.GetLength();++l) {
	h_psi_nocut->Fill(psi_nocut[l].M());
      }
      
      tree->GetEntry(i-1);
      TVector3 mcVertex, mcD1Vertex, mcD2vertex;
      
      if (((PndMCTrack*)mc_array->At(0))!=0) mcD1Vertex = ((PndMCTrack*)mc_array->At(0))->GetStartVertex();
      else
	cout << "Not found k1!" << endl;
      if (((PndMCTrack*)mc_array->At(3))!=0) mcD2Vertex = ((PndMCTrack*)mc_array->At(3))->GetStartVertex();
      else
	cout << "Not found k2!" << endl;
      
      evthead->GetVertex(mcVertex);
      
      // MC PID
      // Leave only kaons in particle lists
      int n_removed=0;
      int ii=0;
      for (l=0;l<p1.GetLength();++l) {
	ii=l-n_removed;
	if (p1[ii].GetMicroCandidate().GetMcIndex()>-1){
	  PndMCTrack *mcTrack = (PndMCTrack*)mc_array->At(p1[ii].GetMicroCandidate().GetMcIndex());
	  if (mcTrack!=0)
	    {
	      if ((mcTrack->GetPdgCode()!=211) || (mcTrack->GetMotherID()!=-1))
		{
		  p1.Remove(p1[ii]);
		  n_removed++;
		}
	    }
	  else
	    {
	      std::cout<<"stt h: " << p1[ii].GetMicroCandidate().GetSttHits() << std::endl;
	      std::cout<<"Kaon list 1, element "<<l<<" has no assosiated mcTRack"<<std::endl;
	    }
	}
      }

      n_removed=0;
      ii=0;
      for (l=0;l<p2.GetLength();++l) {
	ii=l-n_removed;
	if (p2[ii].GetMicroCandidate().GetMcIndex()>-1){
	  PndMCTrack *mcTrack = (PndMCTrack*)mc_array->At(p2[ii].GetMicroCandidate().GetMcIndex());
	  if (mcTrack!=0)
	    {
	      if ((mcTrack->GetPdgCode()!=-211) || (mcTrack->GetMotherID()!=-1))
		{
		  p2.Remove(p2[ii]);
		  n_removed++;
		} 
	    }
	  else
	    { 
	      std::cout<<"stt h: " << p2[ii].GetMicroCandidate().GetSttHits() << std::endl;
	      std::cout<<"Kaon list 2, element "<<l<<" has no assosiated mcTRack"<<std::endl;
	    }
	}
      }
      
      n_removed=0;
      ii=0;
      for (l=0;l<k1.GetLength();++l) {
	ii=l-n_removed;
	if (k1[ii].GetMicroCandidate().GetMcIndex()>-1){
	  PndMCTrack *mcTrack = (PndMCTrack*)mc_array->At(k1[ii].GetMicroCandidate().GetMcIndex());
	  if (mcTrack!=0)
	    {
	      if ((mcTrack->GetPdgCode()!=321) || (mcTrack->GetMotherID()!=-1))
		{
		  k1.Remove(k1[ii]);
		  n_removed++;
		} 
	    }
	  else
	    { 
	      std::cout<<"stt h: " << k1[ii].GetMicroCandidate().GetSttHits() << std::endl;
	      std::cout<<"Kaon list 2, element "<<l<<" has no assosiated mcTRack"<<std::endl;
	    }
	}
      }
      
      n_removed=0;
      ii=0;
      for (l=0;l<k2.GetLength();++l) {
	ii=l-n_removed;
	if (k2[ii].GetMicroCandidate().GetMcIndex()>-1){
	  PndMCTrack *mcTrack = (PndMCTrack*)mc_array->At(k2[ii].GetMicroCandidate().GetMcIndex());
	  if (mcTrack!=0)
	    {
	      if ((mcTrack->GetPdgCode()!=-321) || (mcTrack->GetMotherID()!=-1))
		{
		  k2.Remove(k2[ii]);
		  n_removed++;
		} 
	    }
	  else
	    { 
	      std::cout<<"stt h: " << k2[ii].GetMicroCandidate().GetSttHits() << std::endl;
	      std::cout<<"Kaon list 2, element "<<l<<" has no assosiated mcTRack"<<std::endl;
	    }
	}
      }
         
      D1_pid.Combine(p1,p1, k2);
      D2_pid.Combine(p2,p2, k1);
     
      for (j=0;j<D1_pid.GetLength();++j) h_mD1_pid->Fill(D1_pid[j].M());
      for (j=0;j<D2_pid.GetLength();++j) h_mD2_pid->Fill(D2_pid[j].M());
      D1_pid.Select(psiMassSel);
      D2_pid.Select(psiMassSel);
      psi.Combine(D1_pid,D2_pid); 
      
      for (l=0;l<psi.GetLength();++l) {
	h_psi_pid->Fill(psi[l].M());
      }
      
      cout << " ";// << endl;
      if (use4cfit)
	{
	  ////////////// 4C-fit ///////////////
	  int best_i=0;
	  double best_chi2=100000;
	  TCandidate *ccfit;// = new TCandidate();
	  double m_D1, m_D2;
	  for (l=0;l<psi.GetLength();++l) {
	    Pnd4CFitter fitter(psi[l],ini);
	    fitter.FitConserveMasses();
	    double chi2=fitter.GetChi2();
	    if (chi2<best_chi2)
	      {
		best_chi2=chi2;
		best_i = l;
		ccfit = (TCandidate*)fitter.FittedCand(psi[l]);
		TCandidate *D1best = fitter.FittedCand(*(psi[l].Daughter(0)));
		TCandidate *D2best = fitter.FittedCand(*(psi[l].Daughter(1)));
		TCandidate *p1best = fitter.FittedCand(*(psi[l].Daughter(0)->Daughter(0)));
		TCandidate *p2best = fitter.FittedCand(*(psi[l].Daughter(0)->Daughter(1)));
		TCandidate *k1best = fitter.FittedCand(*(psi[l].Daughter(0)->Daughter(2)));
		TCandidate *p3best = fitter.FittedCand(*(psi[l].Daughter(1)->Daughter(0)));
		TCandidate *p4best = fitter.FittedCand(*(psi[l].Daughter(1)->Daughter(1)));
		TCandidate *k2best = fitter.FittedCand(*(psi[l].Daughter(1)->Daughter(2)));
		TLorentzVector tlvk1 = p1best->P4();
		TLorentzVector tlvk2 = p2best->P4();
		TLorentzVector tlvk3 = k1best->P4();
		TLorentzVector tlvk4 = p3best->P4();
		TLorentzVector tlvk5 = p4best->P4();
		TLorentzVector tlvk6 = k2best->P4();
		m_D1= (tlvk1+tlvk2+tlvk3).M();
		m_D2= (tlvk4+tlvk5+tlvk6).M();
	      }
	    h_chi2_4c->Fill(chi2/9); // Ndf=3N-3=9
	  }
	  
	  if (psi.GetLength()!=0) cout << "evt: " << i << "\tbest: " << best_chi2 << "\tlen " << psi.GetLength() << endl;  
	 
	  if(/*(best_chi2<270)&&*/(psi.GetLength()!=0))
	    {
	      h_chi2b_4c->Fill(best_chi2/9); // Ndf=3N-3=9
	      h_psi_4c->Fill(ccfit->M());
	      h_mD_4c->Fill(m_D1);
	      h_mD_4c->Fill(m_D2);
	      h_mD_final->Fill(m_D1);
	      h_mD_final->Fill(m_D2);
	      if (((m_D1>1.87-0.07)&&(m_D1<1.87+0.07))&&((m_D2>1.87-0.07)&&(m_D2<1.87+0.07)))
	      	{
	      	  h_psi_Dmass->Fill(psi[best_i].M());
	      	  if ((psi[best_i].M()>3.7)&&(psi[best_i].M()<3.9))
	      	    n_reco++;
	      	}
	    }
	}
      else // use vertex fit
	{ 
	  TCandList psi_vtx, D1_vtx, D2_vtx;
	  TCandidate *ck1, *ck2, *cp1, *cp2, *cp3, *cp4, *D1tmp, *D2tmp, *D1_tmp, *D2_tmp, *psi_tmp;
	  
	  //Combine 4 kaons directly to candidates
	  for (j=0;j<psi.GetLength();++j)
	    {
	      D1tmp=psi[j].Daughter(0);
	      D2tmp=psi[j].Daughter(1);

	      cp1=D1tmp->Daughter(0);
	      cp2=D1tmp->Daughter(1); 
	      ck1=D1tmp->Daughter(2);
	      
	      cp3=D2tmp->Daughter(0);
	      cp4=D2tmp->Daughter(1);
	      ck2=D2tmp->Daughter(2);
	      D1_tmp=cp1->Combine(*cp2,*ck1); 
	      D2_tmp=cp3->Combine(*cp4,*ck2);	
	      //psi_tmp=cp1->Combine(*cp2,*ck1,*cp3, *cp4, *ck2);
	      D1_vtx.Add(*D1_tmp); 
	      D2_vtx.Add(*D2_tmp);
	      //psi_vtx.Add(*psi_tmp);
	    }

	  for (l=0;l<psi.GetLength();++l) {
	    // h_psi_pid_vtx->Fill(psi_vtx[l].M());
	  }

	  {
	    ////////////// Vertex fit D1 /////////////
	    int best_i=0;
	    double best_chi2=1000;
	    TCandidate *D1fit_best=0;
	    TCandidate *k1fit_best, *p1fit_best, *p2fit_best, *k1nofit, *p1nofit, *p2nofit;
	    TVector3 bestPos;
	    Float_t D1vtx_mass, pre_m_D1;
	    for (j=0;j<D1_vtx.GetLength();++j)
	      {
		PndKinVtxFitter vtxfitter(D1_vtx[j]);        // instantiate a vertex fitter
		vtxfitter.Fit();                          // do the vertex fit
		
		TCandidate *D1fit=vtxfitter.FittedCand(D1_vtx[j]);  // request the fitted Psi candidate
		TVector3 D1Vtx=D1fit->Pos();                    // and the decay vertex position
		double chi2_vtx_D1=vtxfitter.GlobalChi2();
		h_chi2_vtx_D1->Fill(chi2_vtx_D1/5); // Number degree of freedom 2N-3=5
		// plot mass and vtx x,y projection after fit
		hvpos_D1->Fill(D1Vtx.X(),D1Vtx.Y());
		hvzpos_D1->Fill(D1Vtx.Z());
		if(chi2_vtx_D1<best_chi2)
		{
		  best_chi2=chi2_vtx_D1;
		  best_i=l;
		  D1fit_best=D1fit;
		  p1fit_best=vtxfitter.FittedCand(*(D1fit_best->Daughter(0)));
		  p2fit_best=vtxfitter.FittedCand(*(D1fit_best->Daughter(1)));
		  k1fit_best=vtxfitter.FittedCand(*(D1fit_best->Daughter(2))); 
		  p1nofit=D1_vtx[j].Daughter(0);
		  p2nofit=D1_vtx[j].Daughter(1);
		  k1nofit=D1_vtx[j].Daughter(2);
		  TLorentzVector tlvp1 = p1nofit->P4();  
		  TLorentzVector tlvp2 = p2nofit->P4(); 
		  TLorentzVector tlvk1 = k1nofit->P4();
		  pre_m_D1= (tlvp1+tlvp2+tlvk1).M();
		  
		  D1vtx_mass = D1fit_best->M();
		  bestPos = D1fit->Pos(); 
		}
	      }
	    if((psi.GetLength()!=0))
	      {
		h_D1_vtx->Fill(D1vtx_mass);
		h_D1_vtx_bef->Fill(pre_m_D1);
		
		hvtxresX_D1->Fill(mcD1Vertex.X()-bestPos.X());
		hvtxresY_D1->Fill(mcD1Vertex.Y()-bestPos.Y());
		hvtxresZ_D1->Fill(mcD1Vertex.Z()-bestPos.Z());
 
	    }
	  }

	  
	  {
	    ////////////// Vertex fit D2 /////////////
	    int best_i=0;
	    double best_chi2=1000;
	    TCandidate *D2fit_best=0;
	    TCandidate *k1fit_best, *p1fit_best, *p2fit_best, *k1nofit, *p1nofit, *p2nofit;
	    TVector3 bestPos;
	    Float_t D2vtx_mass, pre_m_D2;
	    for (j=0;j<D2_vtx.GetLength();++j)
	      {
		PndKinVtxFitter vtxfitter(D2_vtx[j]);        // instantiate a vertex fitter
		vtxfitter.Fit();                          // do the vertex fit
		
		TCandidate *D2fit=vtxfitter.FittedCand(D2_vtx[j]);  // request the fitted Psi candidate
		TVector3 D2Vtx=D2fit->Pos();                    // and the decay vertex position
		double chi2_vtx_D2=vtxfitter.GlobalChi2();
		h_chi2_vtx_D2->Fill(chi2_vtx_D2/5); // Number degree of freedom 2N-3=5
		// plot mass and vtx x,y projection after fit
		hvpos_D2->Fill(D2Vtx.X(),D2Vtx.Y());
		hvzpos_D2->Fill(D2Vtx.Z());
		if(chi2_vtx_D2<best_chi2)
		{
		  best_chi2=chi2_vtx_D2;
		  best_i=l;
		  D2fit_best=D2fit;
		  p1fit_best=vtxfitter.FittedCand(*(D2fit_best->Daughter(0)));
		  p2fit_best=vtxfitter.FittedCand(*(D2fit_best->Daughter(1)));
		  k1fit_best=vtxfitter.FittedCand(*(D2fit_best->Daughter(2)));
		  p1nofit=D2_vtx[j].Daughter(0);
		  p2nofit=D2_vtx[j].Daughter(1);
		  k1nofit=D2_vtx[j].Daughter(2);
		  TLorentzVector tlvp1 = p1nofit->P4();  
		  TLorentzVector tlvp2 = p2nofit->P4(); 
		  TLorentzVector tlvk1 = k1nofit->P4();
		  pre_m_D2= (tlvp1+tlvp2+tlvk1).M();
		  
		  D2vtx_mass = D2fit_best->M();
		  bestPos = D2fit->Pos(); 
		}
	      }
	    if((psi.GetLength()!=0))
	      {
		h_D2_vtx->Fill(D2vtx_mass);
		h_D2_vtx_bef->Fill(pre_m_D2);
		hvtxresX_D2->Fill(mcD2Vertex.X()-bestPos.X());
		hvtxresY_D2->Fill(mcD2Vertex.Y()-bestPos.Y());
		hvtxresZ_D2->Fill(mcD2Vertex.Z()-bestPos.Z());
 
	    }
	  }
	}
      
    }
  std::cout<<"Number of reconstructed eta_c = "<<n_reco<<std::endl;
  n_psi->SetBinContent(1,n_reco);


  out->cd();
  
  h_mD1_nocuts->Write();
  h_mD1_pid->Write();
  h_mD2_nocuts->Write();
  h_mD2_pid->Write();
  h_psi_nocut->Write();
  h_psi_pid->Write();
  
  
  n_psi->Write();
  n_events->Write();
 
  h_psi_Dmass->Write();
  h_psi_vtx->Write(); 
  h_D1_vtx_bef->Write(); 
  h_D2_vtx_bef->Write();
  h_D1_vtx->Write(); 
  h_D2_vtx->Write();
  
  h_psi_4c->Write();

  //h_mD_vtx->Write();
  h_mD_4c->Write();
  h_mD_final->Write();

  nc->Write();

  h_chi2_4c->Write(); 
  h_chi2b_4c->Write();
  h_chi2_vtx_D1->Write();
  hvzpos_D1->Write();
  hvpos_D1->Write();
  hvtxresX_D1->Write();
  hvtxresY_D1->Write(); 
  hvtxresZ_D1->Write();
  hvzpos_D2->Write();
  hvpos_D2->Write();
  hvtxresX_D2->Write();
  hvtxresY_D2->Write(); 
  hvtxresZ_D2->Write();

  
  out->Save();

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

}

void removeCombinatoric(TCandList &psi_list)
{
  TCandidate *phi1, *phi2, *k1, *k2, *k3, *k4;
  bool isOverlap;

  int len1=psi_list.GetLength();
  int ii=0;

  for (int i=0; i<len1; ++i)
    {
      phi1=psi_list[ii].Daughter(0);
      phi2=psi_list[ii].Daughter(1);
      k1=phi1->Daughter(0);
      k2=phi1->Daughter(1);
      k3=phi2->Daughter(0);
      k4=phi2->Daughter(1);

      if ((k1->IsCloneOf(*k3))||(k2->IsCloneOf(*k4)))
	{
	  psi_list.Remove(psi_list[ii]);
	  std::cout<<"Overlap found, i="<<i<<std::endl;
	}
      else
	ii++;
    }

}



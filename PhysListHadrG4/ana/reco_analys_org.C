{
  // Macro loads a file after reconstruction and plots difference between initial direction of particle and angular position of cluster
	gROOT->SetStyle("Plain");
	gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
	gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
	rootlogon();
	basiclibs();
	
	TFile* f = new TFile("cluster_emc.root"); //file you want to analyse
	TTree *t=(TTree *) f->Get("cbmsim") ;
	TClonesArray* cluster_array=new TClonesArray("PndEmcCluster");
	t->SetBranchAddress("EmcCluster",&cluster_array);
	
	t->AddFriend("cbmsim", "digi_emc.root");
	TClonesArray* digi_array=new TClonesArray("PndEmcDigi");
	t->SetBranchAddress("EmcDigi",&digi_array);

	TFile* fsim = new TFile("sim_emc.root"); //file you want to analyse
	TTree *tsim=(TTree *) fsim->Get("cbmsim") ;
	
	PndEmcMapper::Init(1);

	TClonesArray* mctrack_array=new TClonesArray("PndMCTrack");
	tsim->SetBranchAddress("MCTrack",&mctrack_array);
	
	TVector3 photon_momentum;

	double cluster_energy;
	double cluster_theta, cluster_phi; //position of the cluster
	double theta, phi; // angular position of the initial particle
	double theta_diff, phi_diff;
	int ndigi, npoint;
	double max_energy=0;
	
	TH1F *h1= new TH1F("h1","Theta difference",200,-5.,5.);
	TH1F *h2= new TH1F("h2","Phi difference",200,-5.,5.);
	TH1F *h3= new TH1F("h3","Cluster energy",100,0.85,1.05);
	TH2F *h2theta= new TH2F("h2theta","Theta difference",200,0.,180.,200,-5.,5.);
	TH2F *h2phi= new TH2F("h2phi","Phi difference",200,0.,180.,200,-5.,5.);
	TH1F *hE1= new TH1F("hE1","E1",200,0.,1.05);
	TH1F *hE1E9= new TH1F("hE1E9","E1 / E9",200,0.,1.05);
	TH1F *hE9E25= new TH1F("hE9E25","E9 / E25",200,0.,1.05);

	// Cluster angular position
	// Entrance point is determined by minimal time
		
	// Cluster energy
	for (Int_t j=0; j< t->GetEntriesFast(); j++)
	{
		t->GetEntry(j);
		for (Int_t i=0; i<cluster_array->GetEntriesFast(); i++)
		{
			PndEmcCluster *cluster=(PndEmcCluster*)cluster_array->At(i);
			cluster_energy=cluster->energy();
			if ((cluster->NumberOfDigis()>1)&&(cluster_energy>0.02))
				h3->Fill(cluster_energy);
			PndEmcClusterEnergySums esum(*cluster, digi_array);
			hE1->Fill(esum.E1());
			hE1E9->Fill(esum.E1E9());
			hE9E25->Fill(esum.E9E25());
			
		}
	}

	for (Int_t j=0; j< t->GetEntriesFast(); j++)//t->GetEntriesFast()
	{
		t->GetEntry(j);
		tsim->GetEntry(j);
		
		PndMCTrack *mctrack=(PndMCTrack *) mctrack_array->At(0);
		photon_momentum=mctrack->GetMomentum();
		theta=photon_momentum.Theta();
		phi=photon_momentum.Phi();
		
	
		// Loop over clusters
		// If we have 1 initial particle and several cluster
		// we can separate cluster from the first interaction by maximum energy
		
		max_energy=0;
		
		for (Int_t i=0; i<cluster_array->GetEntriesFast(); i++)
		{
			PndEmcCluster *cluster=(PndEmcCluster*)cluster_array->At(i);
			cluster_energy=cluster->energy();
			if (cluster_energy>max_energy)
			{
				max_energy=cluster_energy;
				TVector3 cluster_pos=cluster->where();
				cluster_theta=cluster_pos.Theta();
				cluster_phi=cluster_pos.Phi();
			}
						
		}
		
		if (max_energy>0.6)
		{
			theta_diff=(cluster_theta-theta)*180./TMath::Pi();
			h1->Fill(theta_diff);
			h2theta->Fill(theta*TMath::RadToDeg(),theta_diff);
			
			phi_diff=(cluster_phi-phi)*180./TMath::Pi();
			h2->Fill(phi_diff);
			h2phi->Fill(phi*TMath::RadToDeg(),phi_diff);
		}
		
	}

	TCanvas* c1 = new TCanvas("c1", "Cluster Energy", 100, 100, 800, 800); 	
	h3->SetTitle("Cluster energy of 1 GeV photon");
	h3->GetXaxis()->SetTitle("Energy, GeV");
	h3->Draw();

	
	TCanvas* c2 = new TCanvas("c2", "#theta_{reco} - #theta_{truth}", 100, 100, 800, 800); 	
	h1->SetTitle("Difference between cluster and initial photon theta angle");
	h1->GetXaxis()->SetTitle("#theta_{reco} - #theta_{truth}, degree");
	h1->Draw();

// 	TF1 *f1 = new TF1("f1","gaus(0)",-2.,2.);
// 	Double_t par1[3]={60,0,0.5};
// 	f1->SetParameters(par1);
// 	f1->SetLineColor(2);
	// 
// 	h1->Fit("f1","RB");
// 	double mu1=f1->GetParameter(1);
// 	double sigma1=f1->GetParameter(2);

	TCanvas* c3 = new TCanvas("c3", "#phi_{reco} - #phi_{truth}", 100, 100, 800, 800); 	
	h2->SetTitle("Difference between cluster and initial photon phi angle");
	h2->GetXaxis()->SetTitle("#phi_{reco} - #phi_{truth}, degree");
	h2->Draw();
	
	TCanvas* c4 = new TCanvas("c4", "#theta_{reco} - #theta_{truth} vs #theta_{truth}", 100, 100, 800, 800); 	
	h2theta->SetTitle("Difference between cluster and initial photon theta angle");
	h2theta->GetXaxis()->SetTitle("#theta_{truth}, degree");
	h2theta->GetYaxis()->SetTitle("#theta_{reco} - #theta_{truth}, degree");
	h2theta->Draw();

	TCanvas* c5 = new TCanvas("c5", "#phi_{reco} - #phi_{truth} vs #phi_{truth}", 100, 100, 800, 800); 	
	h2phi->SetTitle("Difference between cluster and initial photon phi angle");
	h2phi->GetXaxis()->SetTitle("#phi_{truth}, degree");
	h2phi->GetYaxis()->SetTitle("#phi_{reco} - #phi_{truth}, degree");
	h2phi->Draw();
	
	TCanvas* c6 = new TCanvas("c6", "Cluster Properties", 100, 100, 800, 800); 
	c6->Divide(2,2);
	c6->cd(1); hE1->Draw();
	c6->cd(2); hE1E9->Draw();	
	c6->cd(3); hE9E25->Draw();


}


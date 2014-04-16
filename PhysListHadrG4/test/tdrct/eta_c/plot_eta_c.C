{
	bool vtxfit=0;
	
	if (vtxfit)
		TString inFile="etac_histo_vtx.root";
	else 
		TString inFile="etac_histo_4c.root";
	
// 	gStyle->SetLabelSize(0.08);
// 	gStyle->SetTitleSize(0.08);
	
	TFile *f = TFile::Open(inFile);
	TVectorD *n_etac = (TVectorD*)f->Get("n_etac");
	TVectorD *n_events = (TVectorD*)f->Get("n_events");

	TH1F *h_etac_nocut=(TH1F *)f->Get("h_etac_nocut");
	h_etac_nocut->UseCurrentStyle(); 
	TH1F *h_etac_pid=(TH1F *)f->Get("h_etac_pid");
	h_etac_pid->UseCurrentStyle(); 
	TH1F *h_etac_phimass=(TH1F *)f->Get("h_etac_phimass");
	h_etac_phimass->UseCurrentStyle(); 
	TH1F *h_etac_vtx=(TH1F *)f->Get("h_etac_vtx");
	h_etac_vtx->UseCurrentStyle(); 
	TH1F *h_etac_4c=(TH1F *)f->Get("h_etac_4c");
	h_etac_4c->UseCurrentStyle(); 

	TH1F *h_mphi_nocuts=(TH1F *)f->Get("h_mphi_nocuts");
	h_mphi_nocuts->UseCurrentStyle(); 
	TH1F *h_mphi_pid=(TH1F *)f->Get("h_mphi_pid");
	h_mphi_pid->UseCurrentStyle(); 
	TH1F *h_mphi_vtx=(TH1F *)f->Get("h_mphi_vtx");
	h_mphi_vtx->UseCurrentStyle(); 
	TH1F *h_mphi_4c=(TH1F *)f->Get("h_mphi_4c");
	h_mphi_4c->UseCurrentStyle(); 
	TH1F *h_mphi_final=(TH1F *)f->Get("h_mphi_final");
	h_mphi_final->UseCurrentStyle();
	
	TH1F *nc=(TH1F *)f->Get("nc");
	nc->UseCurrentStyle(); 
	 
	TH1F *h_chi2_4c=(TH1F *)f->Get("h_chi2_4c");
	h_chi2_4c->UseCurrentStyle(); 
	TH1F *h_chi2_vtx=(TH1F *)f->Get("h_chi2_vtx");
	h_chi2_vtx->UseCurrentStyle(); 
	TH1F *hvzpos=(TH1F *)f->Get("hvzpos");
	hvzpos->UseCurrentStyle(); 
	TH2F *hvpos=(TH2F *)f->Get("hvpos");
	hvpos->UseCurrentStyle(); 

	TCanvas *c1=new TCanvas("c1","c1",600,600);
	nc->Draw();

	TCanvas *c2=new TCanvas("c2","c2",600,600);
	c2->Divide(1,2);
	c2->cd(1);
	h_mphi_nocuts->Draw();
	c2->cd(2);
	h_etac_nocut->Draw();
	
	TCanvas *c3=new TCanvas("c3","c3",600,600);
	c3->Divide(1,2);
	c3->cd(1);
	h_mphi_pid->Draw();
	c3->cd(2);
	h_etac_pid->Draw();
	
	TCanvas *c4=new TCanvas("c4","c4",600,600);
	c4->Divide(1,2);
	c4->cd(1);
	h_mphi_4c->Draw();
	c4->cd(2);
	h_etac_4c->Draw();

	TCanvas *c5=new TCanvas("c5","c5",600,600);
	c5->Divide(1,2);
	c5->cd(1);
	h_mphi_final->Draw();
	
	double mean_phi=1.1;
	double rms_phi=0.1;
	
	TF1 *f1_phi = new TF1("f1_phi","gaus",0.9,1.1);
	h_mphi_final->Fit(f1_phi,"R","",mean_phi-1.6*rms_phi,mean_phi+1.6*rms_phi);
	
	double sigma1_phi=f1_phi->GetParameter(2);
	double mean1_phi=f1_phi->GetParameter(1);
	
	TF1 *f2_phi = new TF1("f2_phi","gaus",0.9,1.1);
	h_mphi_final->Fit(f2_phi,"R","",mean1_phi-1.6*sigma1_phi,mean1_phi+1.6*sigma1_phi);

	double sigma2_phi=f2_phi->GetParameter(2);
	std::cout<<"sigma phi="<<sigma2_phi<<std::endl;

	c5->cd(2);
	h_etac_phimass->Draw();
	// fit eta_c
	double mean=h_etac_phimass->GetMean();
	double rms=h_etac_phimass->GetRMS();
	
	TF1 *f1 = new TF1("f1","gaus",2.8,3.2);
	h_etac_phimass->Fit(f1,"R","",mean-1.6*rms,mean+1.6*rms);
	
	double sigma1=f1->GetParameter(2);
	double mean1=f1->GetParameter(1);
	
	TF1 *f2 = new TF1("f2","gaus",2.8,3.2);
	h_etac_phimass->Fit(f2,"R","",mean1-1.6*sigma1,mean1+1.6*sigma1);

	double sigma2=f2->GetParameter(2);
	std::cout<<"sigma="<<sigma2<<std::endl;

	TCanvas *c6=new TCanvas("c6","c6",600,600);
	h_chi2_4c->Draw();
	
	if (vtxfit)
	{
		TCanvas *c7=new TCanvas("c7","c7",600,600);
		hvpos->Draw();

		TCanvas *c8=new TCanvas("c8","c8",600,600);
		hvzpos->Draw();
		
		TCanvas *c9=new TCanvas("c9","c9",600,600);
		h_chi2_vtx->Draw();
		
		TCanvas *c10=new TCanvas("c10","c10",600,600);
		c10->Divide(1,2);
		c10->cd(1);
		h_mphi_vtx->Draw();
		c10->cd(2);
		h_etac_vtx->Draw();
	}
}
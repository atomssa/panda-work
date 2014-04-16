{	
		TString inFile="invariantmass_2pi_stt_a.root";
	

	TFile *f = TFile::Open(inFile);

	TCanvas *c2=new TCanvas("c2","c2",600,600);
	c2->cd();
	chivtx->Draw();

	TCanvas *c3=new TCanvas("c3","c3",600,600);
	c3->cd();
	hvpos->Draw();

	TCanvas *c4=new TCanvas("c4","c3",600,600);
	c4->cd();
	hvzpos->Draw();


	TCanvas *c5=new TCanvas("c5","c3",600,600);
        c5->cd();
        TF1 *f1 = new TF1("f1","gaus",-0.015,0.015);
        hvtxresX->Fit(f1,"R");

	TCanvas *c6=new TCanvas("c6","c6",600,600);
        c6->cd();
        TF1 *f1 = new TF1("f1","gaus",-0.016,0.016);
        hvtxresY->Fit(f1,"R");

	TCanvas *c7=new TCanvas("c7","c7",600,600);
        c7->cd();
        TF1 *f1 = new TF1("f1","gaus",0.50,-0.50);
        hvtxresZ->Fit(f1,"R");



	TCanvas *c8=new TCanvas("c8","c8",600,600);
	c8->cd();
	TF1 *f1 = new TF1("f1","gaus",2.97,3.18);
	invmasschicut_best->Fit(f1,"R");

	TCanvas *c10=new TCanvas("c10","c10",600,600);
        c10->cd();
        TF1 *f1 = new TF1("f1","gaus",2.97,3.18);
        invmasswithpid_sel->Fit(f1,"R");

	TCanvas *c11=new TCanvas("c11","c11",600,600);
        c11->cd();
        TF1 *f1 = new TF1("f1","gaus",2.97,3.18);
        invmassnocut->Fit(f1,"R");



	
}

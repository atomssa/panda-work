{
	bool saveHistos=false;
	TString inFile="etac_histo_stt.root";
	
	TFile *f = TFile::Open(inFile);
	
	TH1F *h_etac_nocut=(TH1F *)f->Get("h_etac_nocut");
	TH1F *h_etac_pid=(TH1F *)f->Get("h_etac_pid");
	TH1F *h_etac_vtx=(TH1F *)f->Get("h_etac_vtx");
	TH1F *h_etac_4c=(TH1F *)f->Get("h_etac_4c");
	TH1F *h_etac_phimass_4c=(TH1F *)f->Get("h_etac_phimass_4c");
	TH1F *h_etac_phimass_vtx=(TH1F *)f->Get("h_etac_phimass_vtx");
	
	TH1F *h_mphi_nocuts=(TH1F *)f->Get("h_mphi_nocuts");
	TH1F *h_mphi_pid=(TH1F *)f->Get("h_mphi_pid");
	TH1F *h_mphi_vtx=(TH1F *)f->Get("h_mphi_vtx");
	TH1F *h_mphi_4c=(TH1F *)f->Get("h_mphi_4c");
	TH1F *h_mphi_final_4c=(TH1F *)f->Get("h_mphi_final_4c");
	TH1F *h_mphi_final_vtx=(TH1F *)f->Get("h_mphi_final_vtx");
	
	TH1F *h_etac_phimass_vtx_2=(TH1F *)f->Get("h_etac_phimass_vtx_2");
	TH1F *h_mphi_vtx_2=(TH1F *)f->Get("h_mphi_vtx_2");
	TH1F *h_mphi_final_vtx_2=(TH1F *)f->Get("h_mphi_final_vtx_2");
	
	TH1F *n_etac_vtx_2=(TH1F *)f->Get("n_etac_vtx_2");
	TH1F *h_chi2_prefit=(TH1F *)f->Get("h_chi2_prefit");
	
	TH1F *nc=(TH1F *)f->Get("nc");
	 
	TH1F *h_chi2_4c=(TH1F *)f->Get("h_chi2_4c");
	TH1F *h_chi2b_4c=(TH1F *)f->Get("h_chi2b_4c");
	
	TH1F *h_chi2_vtx=(TH1F *)f->Get("h_chi2_vtx");
	TH1F *h_chi2b_vtx=(TH1F *)f->Get("h_chi2b_vtx");
	TH1F *hvzpos=(TH1F *)f->Get("hvzpos");
	TH2F *hvpos=(TH2F *)f->Get("hvpos");
	
	TH1F *hvtxresX=(TH1F *)f->Get("hvtxresX");
	TH1F *hvtxresY=(TH1F *)f->Get("hvtxresY");
	TH1F *hvtxresZ=(TH1F *)f->Get("hvtxresZ");
	
	TH1F *n_etac_4c=(TH1F *)f->Get("n_etac_4c");
	TH1F *n_etac_vtx=(TH1F *)f->Get("n_etac_vtx");
	TH1F *n_events=(TH1F *)f->Get("n_events");
	
	double eff1=n_etac_4c->GetBinContent(1)/n_events->GetBinContent(1);
	std::cout<<"Efficiency (4C-fit) = "<<eff1<<std::endl;
	double eff2=n_etac_vtx->GetBinContent(1)/n_events->GetBinContent(1);
	std::cout<<"Efficiency (vertex fit) = "<<eff2<<std::endl;

	double eff3=n_etac_vtx_2->GetBinContent(1)/n_events->GetBinContent(1);
	std::cout<<"Efficiency (vertex fit, prefit selection) = "<<eff3<<std::endl;
	
	// Max efficiency
	// Numebr of events with >=4 reconstructed tracks
	std::cout<<"Max efficiency="<<nc->Integral(5,20)/nc->GetEntries()<<std::endl;

	TCanvas *c1=new TCanvas("c1","N charged",600,600);
	nc->Draw();
	if (saveHistos) c1->SaveAs("n_charged_stt.png");

	TCanvas *c2=new TCanvas("c2","No cuts",600,600);
	c2->Divide(1,2);
	c2->cd(1);
	h_mphi_nocuts->SetTitleSize(18);
	h_mphi_nocuts->GetXaxis()->SetTitleSize(0.05);
	h_mphi_nocuts->GetXaxis()->SetLabelSize(0.06);
	h_mphi_nocuts->GetXaxis()->SetTitleSize(0.05);
	h_mphi_nocuts->GetYaxis()->SetLabelSize(0.06);
	h_mphi_nocuts->Draw();
	TLine *l1=new TLine(1.12,0,1.12,35000);
	l1->SetLineColor(2);
	l1->SetLineWidth(2);
	l1->Draw();
	
	c2->cd(2);
	h_etac_nocut->SetTitleSize(18);
	h_etac_nocut->GetXaxis()->SetTitleSize(0.05);
	h_etac_nocut->GetXaxis()->SetLabelSize(0.06);
	h_etac_nocut->GetXaxis()->SetTitleSize(0.05);
	h_etac_nocut->GetYaxis()->SetLabelSize(0.06);
	h_etac_nocut->Draw();
	if (saveHistos) c2->SaveAs("m_nocuts_stt.png");
	
	TCanvas *c3=new TCanvas("c3","MC PID",600,600);
	c3->Divide(1,2);
	c3->cd(1);
	h_mphi_pid->SetTitleSize(18);
	h_mphi_pid->GetXaxis()->SetTitleSize(0.05);
	h_mphi_pid->GetXaxis()->SetLabelSize(0.06);
	h_mphi_pid->GetXaxis()->SetTitleSize(0.05);
	h_mphi_pid->GetYaxis()->SetLabelSize(0.06);
	h_mphi_pid->Draw();

	TLine *l2=new TLine(1.12,0,1.12,35000);
	l2->SetLineColor(2);
	l2->SetLineWidth(2);
	l2->Draw();

	c3->cd(2);
	h_etac_pid->SetTitleSize(18);
	h_etac_pid->GetXaxis()->SetTitleSize(0.05);
	h_etac_pid->GetXaxis()->SetLabelSize(0.06);
	h_etac_pid->GetXaxis()->SetTitleSize(0.05);
	h_etac_pid->GetYaxis()->SetLabelSize(0.06);
	h_etac_pid->Draw();
	if (saveHistos) c3->SaveAs("m_pid_stt.png");

	//////////// 4C-fit fit ////////////////////////////////
// 	TCanvas *c4=new TCanvas("c4","chi2 (4C-fit)",600,600);
// 	c4->Divide(2,1);
// 	c4->cd(1);
// 	h_chi2_4c->Draw();
// 	c4->cd(2);
// 	h_chi2b_4c->Draw();
// 	if (saveHistos) c4->SaveAs("chi2_4c_stt.png");
// 
// 	TCanvas *c5=new TCanvas("c5","m (4C-fit)",600,600);
// 	c5->Divide(1,2);
// 	c5->cd(1);
// 	h_mphi_4c->Draw();
// 	c5->cd(2);
// 	h_etac_4c->Draw();
// 	if (saveHistos) c5->SaveAs("m_4c_stt.png");

	//////////// Vetrex fit ////////////////////////////////
	TCanvas *c6=new TCanvas("c6","Vertex position",600,600);
	c6->Divide(2,1);
	c6->cd(1);
	hvpos->Draw();
	c6->cd(2);
	hvzpos->Draw();
	if (saveHistos) c6->SaveAs("vertex_pos_stt.png");
	
// 	TCanvas *c7=new TCanvas("c7","Vertex fit chi2",600,600);
// 	c7->Divide(2,1);
// 	c7->cd(1);
// 	h_chi2_vtx->Draw();
// 	c7->cd(2);
// 	h_chi2b_vtx->Draw();
// 	if (saveHistos) c7->SaveAs("chi2_vtx_stt.png");
	
	TCanvas *c8=new TCanvas("c8","Vertex resolution",600,600);
	c8->Divide(2,2);
	c8->cd(1);
	hvtxresX->Draw();

	TF1 *f1_vtxx= new TF1("f1_vtxx","gaus",-1.,1.);
	hvtxresX->Fit(f1_vtxx,"R","",-1,1.);
	double mean_x=f1_vtxx->GetParameter(1);
	double sigma_x=f1_vtxx->GetParameter(2);
	TF1 *f2_vtxx= new TF1("f2_vtxx","gaus",-1.,1.);
	hvtxresX->Fit(f2_vtxx,"R","",mean_x-1.6*sigma_x,mean_x+1.6*sigma_x);
	double mean_x2=f2_vtxx->GetParameter(1);
	double sigma_x2=f2_vtxx->GetParameter(2);
	std::cout<<"!!!!!!!!!! vertex x resolution = "<<sigma_x2<<std::endl;
	
	c8->cd(2);
	hvtxresY->Draw();
	TF1 *f1_vtxy= new TF1("f1_vtxy","gaus",-1.,1.);
	hvtxresY->Fit(f1_vtxy,"R","",-1,1.);
	double mean_y=f1_vtxy->GetParameter(1);
	double sigma_y=f1_vtxy->GetParameter(2);
	TF1 *f2_vtxy= new TF1("f2_vtxy","gaus",-1.,1.);
	hvtxresY->Fit(f2_vtxy,"R","",mean_y-1.6*sigma_y,mean_y+1.6*sigma_y);
	double mean_y2=f2_vtxy->GetParameter(1);
	double sigma_y2=f2_vtxy->GetParameter(2);
	std::cout<<"!!!!!!!!!! vertex y resolution = "<<sigma_y2<<std::endl;

	c8->cd(3);
	hvtxresZ->Draw();
	TF1 *f1_vtxz= new TF1("f1_vtxz","gaus",-1.,1.);
	hvtxresZ->Fit(f1_vtxz,"R","",-1,1.);
	double mean_z=f1_vtxz->GetParameter(1);
	double sigma_z=f1_vtxz->GetParameter(2);
	TF1 *f2_vtxz= new TF1("f2_vtxz","gaus",-1.,1.);
	hvtxresZ->Fit(f2_vtxz,"R","",mean_z-1.6*sigma_z,mean_z+1.6*sigma_z);
	double mean_z2=f2_vtxz->GetParameter(1);
	double sigma_z2=f2_vtxz->GetParameter(2);
	std::cout<<"!!!!!!!!!! vertex z resolution = "<<sigma_z2<<std::endl;
	
	
	if (saveHistos) c8->SaveAs("vertex_res_stt.png");
	
	TCanvas *c9=new TCanvas("c9","m vertex",600,600);
	c9->Divide(1,2);
	c9->cd(1);
	h_mphi_vtx->SetTitleSize(18);
	h_mphi_vtx->GetXaxis()->SetTitleSize(0.05);
	h_mphi_vtx->GetXaxis()->SetLabelSize(0.06);
	h_mphi_vtx->GetXaxis()->SetTitleSize(0.05);
	h_mphi_vtx->GetYaxis()->SetLabelSize(0.06);
	h_mphi_vtx->Draw();
	c9->cd(2);
	h_etac_vtx->SetTitleSize(18);
	h_etac_vtx->GetXaxis()->SetTitleSize(0.05);
	h_etac_vtx->GetXaxis()->SetLabelSize(0.06);
	h_etac_vtx->GetXaxis()->SetTitleSize(0.05);
	h_etac_vtx->GetYaxis()->SetLabelSize(0.06);
	h_etac_vtx->Draw();
	if (saveHistos) c9->SaveAs("m_vtx_stt.png");


	double mean_phi, range_phi, sigma1_phi, mean1_phi, sigma2_phi, mean_etac, range_etac;
	double sigma1, mean1, sigma2;
	
	//////////////////// 4C fit ////////////////
// 	TCanvas *c10=new TCanvas("c10","Mass final (4C-fit)",600,600);
// 	c10->Divide(1,2);
// 	c10->cd(1);
// 	h_mphi_final_4c->Draw();
// 	
// 	mean_phi=1.02;
// 	range_phi=0.02;
// 	
// 	TF1 *f1_phi_4c = new TF1("f1_phi_4c","gaus",0.9,1.1);
// 	h_mphi_final_4c->Fit(f1_phi_4c,"R","",mean_phi-range_phi,mean_phi+range_phi);
// 	
// 	sigma1_phi=f1_phi_4c->GetParameter(2);
// 	mean1_phi=f1_phi_4c->GetParameter(1);
// 	
// 	TF1 *f2_phi_4c = new TF1("f2_phi_4c","gaus",0.9,1.1);
// 	h_mphi_final_4c->Fit(f2_phi_4c,"R","",mean1_phi-1.6*sigma1_phi,mean1_phi+1.6*sigma1_phi);
// 
// 	sigma2_phi=f2_phi_4c->GetParameter(2);
// 	std::cout<<"!!!!!!!!!!!!!!!! sigma phi (4c-fit)="<<sigma2_phi<<std::endl;
// 
// 	c10->cd(2);
// 	h_etac_phimass_4c->Draw();
// 	// fit eta_c
// 	mean_etac=2.98;
// 	range_etac=0.1;
// 	
// 	TF1 *f1_4c = new TF1("f1_4c","gaus",2.8,3.2);
// 	h_etac_phimass_4c->Fit(f1_4c,"R","",mean_etac-range_etac,mean_etac+range_etac);
// 	
// 	sigma1=f1_4c->GetParameter(2);
// 	mean1=f1_4c->GetParameter(1);
// 	
// 	TF1 *f2_4c = new TF1("f2_4c","gaus",2.8,3.2);
// 	h_etac_phimass_4c->Fit(f2_4c,"R","",mean1-1.6*sigma1,mean1+1.6*sigma1);
// 
// 	sigma2=f2_4c->GetParameter(2);
// 	std::cout<<"!!!!!!!!!!!!! sigma eta_c (4c-fit) ="<<sigma2<<std::endl;
// 	
// 	if (saveHistos) c10->SaveAs("m_final_4c_stt.png");
	
	//////////////////// Vertex fit ////////////////
	TCanvas *c11=new TCanvas("c11","Mass final (Vertex fit)",600,600);
	c11->Divide(1,2);
	c11->cd(1);
	h_mphi_final_vtx->SetTitleSize(18);
	h_mphi_final_vtx->GetXaxis()->SetTitleSize(0.05);
	h_mphi_final_vtx->GetXaxis()->SetLabelSize(0.06);
	h_mphi_final_vtx->GetXaxis()->SetTitleSize(0.05);
	h_mphi_final_vtx->GetYaxis()->SetLabelSize(0.06);

	h_mphi_final_vtx->Draw();
	
	mean_phi=1.02;
	range_phi=0.02;
	
	TF1 *f1_phi_vtx = new TF1("f1_phi_vtx","gaus",0.9,1.1);
	h_mphi_final_vtx->Fit(f1_phi_vtx,"R","",mean_phi-range_phi,mean_phi+range_phi);
	
	sigma1_phi=f1_phi_vtx->GetParameter(2);
	mean1_phi=f1_phi_vtx->GetParameter(1);
	
	TF1 *f2_phi_vtx = new TF1("f2_phi_vtx","gaus",0.9,1.1);
	h_mphi_final_vtx->Fit(f2_phi_vtx,"R","",mean1_phi-1.6*sigma1_phi,mean1_phi+1.6*sigma1_phi);

	sigma2_phi=f2_phi_vtx->GetParameter(2);
	std::cout<<"!!!!!!!!!!!! sigma phi (vertex fit)="<<sigma2_phi<<std::endl;
	
	TLine *l3=new TLine(1.0,0,1.0,7000);
	l3->SetLineColor(4);
	l3->SetLineWidth(2);
	l3->Draw();

	TLine *l4=new TLine(1.04,0,1.04,7000);
	l4->SetLineColor(4);
	l4->SetLineWidth(2);
	l4->Draw();
	

	c11->cd(2);
	h_etac_phimass_vtx->SetTitleSize(18);
	h_etac_phimass_vtx->GetXaxis()->SetTitleSize(0.05);
	h_etac_phimass_vtx->GetXaxis()->SetLabelSize(0.06);
	h_etac_phimass_vtx->GetXaxis()->SetTitleSize(0.05);
	h_etac_phimass_vtx->GetYaxis()->SetLabelSize(0.06);
	
	h_etac_phimass_vtx->Draw();
	// fit eta_c
	mean_etac=2.98;
	range_etac=0.1;
	
	TF1 *f1_vtx = new TF1("f1_vtx","gaus",2.8,3.2);
	h_etac_phimass_vtx->Fit(f1_vtx,"R","",mean_etac-range_etac,mean_etac+range_etac);
	
	sigma1=f1_vtx->GetParameter(2);
	mean1=f1_vtx->GetParameter(1);
	
	TF1 *f2_vtx = new TF1("f2_vtx","gaus",2.8,3.2);
	h_etac_phimass_vtx->Fit(f2_vtx,"R","",mean1-1.6*sigma1,mean1+1.6*sigma1);

	sigma2=f2_vtx->GetParameter(2);
	std::cout<<"!!!!!!!!!!!!!!! sigma eta_c (vertex fit) ="<<sigma2<<std::endl;
	
	TLine *l5=new TLine(2.9,0,2.9,1200);
	l5->SetLineColor(4);
	l5->SetLineWidth(2);
	l5->Draw();

	TLine *l6=new TLine(3.06,0,3.06,1200);
	l6->SetLineColor(4);
	l6->SetLineWidth(2);
	l6->Draw();
	
	if (saveHistos) c11->SaveAs("m_final_vtx_stt.png");
	
	//////////////////// Vertex fit ////////////////////////////////////
	/////////////////// Prefit best candidate selection ////////////////
// 	TCanvas *c54=new TCanvas("c54","chi2 (prefit selection)",600,600);
// 	h_chi2_prefit->Draw();
// 
// 	TCanvas *c51=new TCanvas("c51","Mass final (Vertex fit)",600,600);
// 	c51->Divide(1,2);
// 	c51->cd(1);
// 	h_mphi_final_vtx_2->SetTitleSize(18);
// 	h_mphi_final_vtx_2->GetXaxis()->SetTitleSize(0.05);
// 	h_mphi_final_vtx_2->GetXaxis()->SetLabelSize(0.06);
// 	h_mphi_final_vtx_2->GetXaxis()->SetTitleSize(0.05);
// 	h_mphi_final_vtx_2->GetYaxis()->SetLabelSize(0.06);
// 
// 	h_mphi_final_vtx_2->Draw();
// 	
// 	mean_phi=1.02;
// 	range_phi=0.02;
// 	
// 	TF1 *f1_phi_vtx_2 = new TF1("f1_phi_vtx_2","gaus",0.9,1.1);
// 	h_mphi_final_vtx_2->Fit(f1_phi_vtx_2,"R","",mean_phi-range_phi,mean_phi+range_phi);
// 	
// 	sigma1_phi=f1_phi_vtx_2->GetParameter(2);
// 	mean1_phi=f1_phi_vtx_2->GetParameter(1);
// 	
// 	TF1 *f2_phi_vtx_2 = new TF1("f2_phi_vtx_2","gaus",0.9,1.1);
// 	h_mphi_final_vtx_2->Fit(f2_phi_vtx_2,"R","",mean1_phi-1.6*sigma1_phi,mean1_phi+1.6*sigma1_phi);
// 
// 	sigma2_phi=f2_phi_vtx_2->GetParameter(2);
// 	std::cout<<"!!!!!!!!!!!! sigma phi (vertex fit)="<<sigma2_phi<<std::endl;
// 	
// 	TLine *l3=new TLine(1.0,0,1.0,7000);
// 	l3->SetLineColor(4);
// 	l3->SetLineWidth(2);
// 	l3->Draw();
// 
// 	TLine *l4=new TLine(1.04,0,1.04,7000);
// 	l4->SetLineColor(4);
// 	l4->SetLineWidth(2);
// 	l4->Draw();
// 	
// 
// 	c51->cd(2);
// 	h_etac_phimass_vtx_2->SetTitleSize(18);
// 	h_etac_phimass_vtx_2->GetXaxis()->SetTitleSize(0.05);
// 	h_etac_phimass_vtx_2->GetXaxis()->SetLabelSize(0.06);
// 	h_etac_phimass_vtx_2->GetXaxis()->SetTitleSize(0.05);
// 	h_etac_phimass_vtx_2->GetYaxis()->SetLabelSize(0.06);
// 	
// 	h_etac_phimass_vtx_2->Draw();
// 	// fit eta_c
// 	mean_etac=2.98;
// 	range_etac=0.1;
// 	
// 	TF1 *f1_vtx_2 = new TF1("f1_vtx_2","gaus",2.8,3.2);
// 	h_etac_phimass_vtx_2->Fit(f1_vtx_2,"R","",mean_etac-range_etac,mean_etac+range_etac);
// 	
// 	sigma1=f1_vtx_2->GetParameter(2);
// 	mean1=f1_vtx_2->GetParameter(1);
// 	
// 	TF1 *f2_vtx_2 = new TF1("f2_vtx_2","gaus",2.8,3.2);
// 	h_etac_phimass_vtx_2->Fit(f2_vtx_2,"R","",mean1-1.6*sigma1,mean1+1.6*sigma1);
// 
// 	sigma2=f2_vtx_2->GetParameter(2);
// 	std::cout<<"!!!!!!!!!!!!!!! sigma eta_c (vertex fit) ="<<sigma2<<std::endl;
// 	
// 	TLine *l5=new TLine(2.9,0,2.9,1200);
// 	l5->SetLineColor(4);
// 	l5->SetLineWidth(2);
// 	l5->Draw();
// 
// 	TLine *l6=new TLine(3.06,0,3.06,1200);
// 	l6->SetLineColor(4);
// 	l6->SetLineWidth(2);
// 	l6->Draw();
	

}

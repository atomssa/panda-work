{	
		TString inFile="out_test_2pi_stt_a.root";
	

	TFile *f = TFile::Open(inFile);
	
	
	gStyle->SetOptStat(kFALSE);
	gStyle->SetPalette(1);
	
	TCanvas *c1=new TCanvas("c1","",800,800);
	c1->Divide(1,3);
	c1->cd(1);
	hMcSttP->SetTitle("Geometrical coverage for the momentum");
	hMcSttP->GetXaxis()->SetTitle("Momentum (GeV)");
	hMcSttP->Draw();
	hMcSttP->SetLineColor(kRed);
	hMcSttP->Draw("same");
	hMcMvdP->SetLineColor(8);
	hMcMvdP->Draw("same");
	hMcGemP->SetLineColor(6);
	hMcGemP->Draw("same");
	
	leg=new TLegend(0.6,0.7,0.89,0.89);
	leg->SetFillColor(kWhite);
	leg->AddEntry(hMcP,"ALL","l");
	leg->AddEntry(hMcSttP,"STT","l");
	leg->AddEntry(hMcMvdP,"MVD","l");
	leg->AddEntry(hMcGemP,"GEM","l");
	leg->Draw();
	c1->cd(2);
	hMcSttTheta->Draw();
	hMcSttTheta->SetTitle("Geometrical coverage for the theta angle");
	hMcSttTheta->GetXaxis()->SetTitle("#theta (degree)");
	hMcSttTheta->SetLineColor(kRed);
	hMcSttTheta->Draw("same");
	hMcMvdTheta->SetLineColor(8);
	hMcMvdTheta->Draw("same");
	hMcGemTheta->SetLineColor(6);
	hMcGemTheta->Draw("same");
	
	leg=new TLegend(0.6,0.7,0.89,0.89);
	leg->SetFillColor(kWhite);
	leg->AddEntry(hMcTheta,"ALL","l");
	leg->AddEntry(hMcSttTheta,"STT","l");
	leg->AddEntry(hMcMvdTheta,"MVD","l");
	leg->AddEntry(hMcGemTheta,"GEM","l");
	leg->Draw();

	c1->cd(3);
	hMcSttPhi->Draw();
	hMcSttPhi->SetTitle("Geometrical coverage for the phi angle");
	hMcSttPhi->GetXaxis()->SetTitle("#phi (degree)");
	hMcSttPhi->SetLineColor(kRed);
	hMcSttPhi->Draw("same");
	hMcMvdPhi->SetLineColor(8);
	hMcMvdPhi->Draw("same");
	hMcGemPhi->SetLineColor(6);
	hMcGemPhi->Draw("same");
	
	leg=new TLegend(0.6,0.7,0.89,0.89);
	leg->SetFillColor(kWhite);
	leg->AddEntry(hMcPhi,"ALL","l");
	leg->AddEntry(hMcSttPhi,"STT","l");
	leg->AddEntry(hMcMvdPhi,"MVD","l");
	leg->AddEntry(hMcGemPhi,"GEM","l");
	leg->Draw();
		
	TCanvas *c4=new TCanvas("c4","",800,800);
	c4->Divide(1,3);
	c4->cd(1);
	nt->Draw("mc_p>>hMcP(100,0,5)");
	nt->Draw("mc_p>>hMcAccSttP(100,0,5)","mc_stt>0","E1");
	hMcAccSttP->Divide(hMcP);
	hMcAccSttP->SetTitle("Geometrical Acceptance");
	hMcAccSttP->GetXaxis()->SetTitle("Momentum (GeV)");
	hMcAccSttP->GetYaxis()->SetTitle("Acceptance");
	c4->cd(2);
	nt->Draw("mc_theta>>hMcTheta(50,0,180)");
	nt->Draw("mc_theta>>hMcAccSttTheta(50,0,180)","mc_stt>0","E1");
	hMcAccSttTheta->Divide(hMcTheta);
	hMcAccSttTheta->SetTitle("Geometrical Acceptance");
	hMcAccSttTheta->GetXaxis()->SetTitle("#theta (degree)");
	hMcAccSttTheta->GetYaxis()->SetTitle("Acceptance");
	c4->cd(3);
	nt->Draw("mc_phi>>hMcPhi(50,-200,200)");
	nt->Draw("mc_phi>>hMcAccSttPhi(50,-200,200)","mc_stt>0","E1");
	hMcAccSttPhi->Divide(hMcPhi);
	hMcAccSttPhi->SetTitle("Geometrical Acceptance");
	hMcAccSttPhi->GetXaxis()->SetTitle("#phi (degree)");
	hMcAccSttPhi->GetYaxis()->SetTitle("Acceptance");
	
	
	TCanvas *c5=new TCanvas("c5","",800,800);
	float numero=ntEvt->GetEntries(); 
	c5->cd();
	TH1F *f1=new TH1F("f1","",100,0,5);
	ntEvt.Draw("acc_stt>>f1");
	f1->SetTitle("Geometrical Global Acceptance");
	f1->GetXaxis()->SetTitle("Primary Tracks");
	f1->GetYaxis()->SetTitle("Acceptance");
	f1->GetYaxis()->SetTitleOffset(1.5);
	f1->Draw();
	f1->Scale(1/numero);


	TCanvas *c6=new TCanvas("c6","",800,800);
	c6->Divide(1,3);
	c6->cd(1);
	nt->Draw("mc_p>>hMcP(100,0,5)");
	nt->Draw("mc_p>>hRecoEffP(100,0,5)","mult>0","E1");
	hRecoEffP->Divide(hMcP);
	hRecoEffP->SetTitle("Pattern Recognition Efficiency");
	hRecoEffP->GetXaxis()->SetTitle("Momentum (GeV)");
	hRecoEffP->GetYaxis()->SetTitle("Efficiency");
	hRecoEffP->Draw();
	c6->cd(2);
	nt->Draw("mc_theta>>hMcTheta(50,0,180)");
	nt->Draw("mc_theta>>hRecoEffTheta(50,0,180)","mult>0","E1");
	hRecoEffTheta->Divide(hMcTheta);
	hRecoEffTheta->SetTitle("Pattern Recognition Efficiency");
	hRecoEffTheta->GetXaxis()->SetTitle("#theta (degree)");
	hRecoEffTheta->GetYaxis()->SetTitle("Efficiency");
	hRecoEffTheta->Draw();
	c6->cd(3);
	nt->Draw("mc_phi>>hMcPhi(50,-200,200)");
	nt->Draw("mc_phi>>hRecoEffPhi(50,-200,200)","mult>0","E1");
	hRecoEffPhi->Divide(hMcPhi);
	hRecoEffPhi->SetTitle("Pattern Recognition Efficiency");
	hRecoEffPhi->GetXaxis()->SetTitle("#phi (degree)");
	hRecoEffPhi->GetYaxis()->SetTitle("Efficiency");
	hRecoEffPhi->Draw();
	

	TCanvas *c14=new TCanvas("c14","",800,800);
        c14->Divide(1,3);
        c14->cd(1);
        nt->Draw("mc_p>>hMcP(100,0,5)");
        nt->Draw("mc_p>>hRecoEffCtP(100,0,5)","mult>0&&(stt>0||tpc>0)","E1");
        hRecoEffCtP->Divide(hMcP);
        hRecoEffCtP->SetTitle("Pattern Recognition Efficiency Only Central Tracker");
        hRecoEffCtP->GetXaxis()->SetTitle("Momentum (GeV)");
        hRecoEffCtP->GetYaxis()->SetTitle("Efficiency");
        hRecoEffCtP->Draw();
        c14->cd(2);
        nt->Draw("mc_theta>>hMcTheta(50,0,180)");
        nt->Draw("mc_theta>>hRecoEffCtTheta(50,0,180)","mult>0&&(stt>0||tpc>0)","E1");
        hRecoEffCtTheta->Divide(hMcTheta);
        hRecoEffCtTheta->SetTitle("Pattern Recognition Efficiency Only Central Tracker");
        hRecoEffCtTheta->GetXaxis()->SetTitle("#theta (degree)");
        hRecoEffCtTheta->GetYaxis()->SetTitle("Efficiency");
	hRecoEffCtTheta->Draw();
        c14->cd(3);
        nt->Draw("mc_phi>>hMcPhi(50,-200,200)");
        nt->Draw("mc_phi>>hRecoEffCtPhi(50,-200,200)","mult>0&&(stt>0||tpc>0)","E1");
        hRecoEffCtPhi->Divide(hMcPhi);
        hRecoEffCtPhi->SetTitle("Pattern Recognition Efficiency Only Central Tracker");
        hRecoEffCtPhi->GetXaxis()->SetTitle("#phi (degree)");
        hRecoEffCtPhi->GetYaxis()->SetTitle("Efficiency");
        hRecoEffCtPhi->Draw();



	TCanvas *c7=new TCanvas("c7","",800,800);

	c7->cd();
	TH1F *f1=new TH1F("f1","",100,0,5);
	ntEvt.Draw("cand>>f1");
	f1->SetTitle("Number of candidates");
	f1->GetXaxis()->SetTitle("Candidates");
	f1->GetYaxis()->SetTitle("1/evt");
	f1->GetYaxis()->SetTitleOffset(1.5);
	f1->Draw();
	f1->Scale(1/numero);

	TCanvas *c8=new TCanvas("c8","",800,800);

	c8->cd();
	TH1F *f1=new TH1F("f1","",100,0,5);
	ntEvt.Draw("eff>>f1");
      	f1->SetTitle("Number of reconstructed primary tracks STT+MVD+GEM");
	f1->GetXaxis()->SetTitle("Primary Tracks");
	f1->GetYaxis()->SetTitle("Efficiency");
	f1->GetYaxis()->SetTitleOffset(1.5);
	f1->Draw();
	f1->Scale(1/numero);

	TCanvas *c11=new TCanvas("c11","",800,800);

        c11->cd();
        TH1F *f1=new TH1F("f1","",100,0,5);
        ntEvt.Draw("effct>>f1");
        f1->SetTitle("Number of reconstructed primary tracks Only Central Tracker");
        f1->GetXaxis()->SetTitle("Primary Tracks");
        f1->GetYaxis()->SetTitle("Efficiency");
        f1->GetYaxis()->SetTitleOffset(1.5);
        f1->Draw();
        f1->Scale(1/numero);



	TCanvas *c9=new TCanvas("c9","",800,800);
	c9->Divide(3,2);
	c9->cd(1);
        nt->Draw("(mc_p-p):mc_p>>hresp_p(100,0,5,100,-0.5,0.5)","mult>0 && (stt>0||tpc>0)","colz");
	hresp_p->SetTitle("Resolution");
	hresp_p->GetXaxis()->SetTitle("MC Momentum (GeV)");
	hresp_p->GetYaxis()->SetTitle("(MC Mom - Reco Mom )");
	hresp_p->GetYaxis()->SetTitleOffset(1.5);

	c9->cd(2);
        nt->Draw("(mc_theta-theta):mc_theta>>hrestheta_theta(100,0,100,100,-2,2)","mult>0 && (stt>0||tpc>0)","colz");
	hrestheta_theta->SetTitle("Resolution");
	hrestheta_theta->GetXaxis()->SetTitle("MC #theta (degree)");
	hrestheta_theta->GetYaxis()->SetTitle("Mc #theta - Reco #theta ");
	hrestheta_theta->GetYaxis()->SetTitleOffset(1.5);

	c9->cd(3);
        nt->Draw("(mc_p-p):mc_theta>>hresp_theta(100,0,100,100,-0.5,0.5)","mult>0 &&  (stt>0||tpc>0)","colz");
	hresp_theta->SetTitle("Resolution");
	hresp_theta->GetXaxis()->SetTitle("MC #theta (degree)");
	hresp_theta->GetYaxis()->SetTitle("MC Mom - Reco Mom");
	hresp_theta->GetYaxis()->SetTitleOffset(1.5);

	c9->cd(4);
        nt->Draw("(mc_theta-theta):mc_p>>hrestheta_p(100,0,5,100,-2,2)","mult>0 &&  (stt>0||tpc>0)","colz");
	hrestheta_p->SetTitle("Resolution");
	hrestheta_p->GetXaxis()->SetTitle("MC Momentum (GeV)");
	hrestheta_p->GetYaxis()->SetTitle("Mc #theta - Reco #theta ");
	hrestheta_p->GetYaxis()->SetTitleOffset(1.5);

	c9->cd(5);
        nt->Draw("(mc_phi-phi):mc_phi>>hresphi_phi(100,-200,200,100,-2,2)","mult>0 &&  (stt>0||tpc>0)","colz");
	hresphi_phi->SetTitle("Resolution");
        hresphi_phi->GetXaxis()->SetTitle("MC Phi (degree)");
        hresphi_phi->GetYaxis()->SetTitle("Mc #phi - Reco #phi ");
        hresphi_phi->GetYaxis()->SetTitleOffset(1.5);

	c9->Update();

	TCanvas *c20=new TCanvas("c20","",800,800);
	c20->cd();
	TH1F *hist=new TH1F("hist","Total Momentum distribution",100,0,10);
        hist->GetXaxis()->SetTitle("MC Momentum (GeV)");
	nt->Draw("p>>hist","mult>0 && (stt>0 || tpc>0)");

	TCanvas *c21=new TCanvas("c21","",800,800);
        c21->cd();
        TH1F *hist2=new TH1F("hist2","#theta Angle Distribution",100,0,180);
        hist2->GetXaxis()->SetTitle("Theta (degree)");        
	nt->Draw("theta>>hist2","mult>0 && (stt>0 || tpc>0)");

	TCanvas *c22=new TCanvas("c22","",800,800);
        c22->cd();
        TH1F *hist3=new TH1F("hist3","(Reconstructed Momentum - MC Momentum )/MC Momentum",100,-1,1);
	hist3->SetStats(kTRUE);
        hist3->GetXaxis()->SetTitle("(Reco-MC)/MC");
        nt->Draw("(mc_p-p)/mc_p>>hist3","mult>0 && (stt>0 || tpc>0)");
	TF1 *f10 = new TF1("f10","gaus",-0.09,0.09);
        hist3->Fit(f10,"R");

	
}

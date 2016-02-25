void tvsthcm() {

  const char* bdir = "/Users/tujuba/panda/work/jpsi_pi0/";

  double umin[nplab] = {0.0}, umax[nplab] = {0.0};
  for (int iplab=0; iplab<nplab; ++iplab) {
    umin[iplab] = limits[iplab][0];
    umax[iplab] = limits[iplab][0]-(tmin[iplab]-tmax[iplab]);
  }

  int col[3] = {1,2,4};

  TLatex *tt[6][nplab];
  for (int iplab = 0; iplab < nplab; ++iplab) {
    for (int i=0; i<6; ++i) {
      tt[i][iplab] = new TLatex();
      tt[i][iplab]->SetNDC(true);
      tt[i][iplab]->SetTextColor(col[iplab]);
      if (tt>3) tt[i][iplab]->SetTextSize(0.7*tt[i][iplab]->GetTextSize());
      //tt[i][iplab]->SetTextColor(col[iplab]);
    }
  }

  TFile *f[3];
  TH2F* tvscost[3];
  TCanvas *tc[4];
  for (int icanv = 0; icanv < 4; ++icanv) tc[icanv] = new TCanvas(Form("tc%d",icanv),Form("tc%d",icanv),1400,1000);
  for (int iplab=nplab-1; iplab>=0; --iplab) {
    f[iplab]= TFile::Open(Form("%s/hists/sig%d_nofilt_noeff.root",bdir,iplab));
    tvscost[iplab] = (TH2F*)f[iplab]->Get("inv_p_t_pbar_pi0_cm_cost_pi0")->Clone(Form("tvscost%d",3-iplab));
    tvscost[iplab]->SetMarkerColor(col[iplab]);
    for (int icanv = 0; icanv < 4; ++icanv) {
      tc[icanv]->cd();
      tvscost[iplab]->SetTitle(";t[GeV^{2}];cos(#theta_{#pi^{0}}^{CM})");
      tvscost[iplab]->Draw(iplab==nplab-1?"":"same");
    }
  }

  tc[3]->cd();
  for (int iplab = 0; iplab < nplab; ++iplab) {
    tt[3][iplab]->DrawLatex(0.15,0.75-iplab*0.1,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));
  }
  tc[3]->Update();
  tc[3]->Print((Form("%s/figs/2015.09.15/validity_iplab3.png",bdir)));
  gSystem->Exec(Form("convert %s/figs/2015.09.15/validity_iplab3.png %s/figs/2015.09.15/validity_iplab3.pdf",bdir,bidr));

  TLine *tl_t[nplab][2][2]; // [plab][vert-horiz][lower-upper]
  TLine *tl_u[nplab][2][2]; // [plab][vert-horiz][lower-upper]
  TBox *tb_t[nplab][2];
  TBox *tb_u[nplab][2];
  for (int iplab=0; iplab<nplab; ++iplab) {
    cout << costh(tmin[iplab], plab[iplab]) << endl;
    tc[iplab]->cd();

    tl_t[iplab][0][0] = new TLine(tmin[iplab], -1.1, tmin[iplab], costh(tmin[iplab], plab[iplab]));
    tl_t[iplab][0][1] = new TLine(tmax[iplab], -1.1, tmax[iplab], costh(tmax[iplab], plab[iplab]));
    tl_t[iplab][1][0] = new TLine(-14.0, costh(tmin[iplab], plab[iplab]), tmin[iplab], costh(tmin[iplab], plab[iplab]));
    tl_t[iplab][1][1] = new TLine(-14.0, costh(tmax[iplab], plab[iplab]), tmax[iplab], costh(tmax[iplab], plab[iplab]));

    tl_u[iplab][0][0] = new TLine(umin[iplab], -1.1, umin[iplab], costh(umin[iplab], plab[iplab]));
    tl_u[iplab][0][1] = new TLine(umax[iplab], -1.1, umax[iplab], costh(umax[iplab], plab[iplab]));
    tl_u[iplab][1][0] = new TLine(-14.0, costh(umin[iplab], plab[iplab]), umin[iplab], costh(umin[iplab], plab[iplab]));
    tl_u[iplab][1][1] = new TLine(-14.0, costh(umax[iplab], plab[iplab]), umax[iplab], costh(umax[iplab], plab[iplab]));

    for (int hv=0; hv<2; ++hv) {
      for (int lu=0; lu<2; ++lu) {
	tl_t[iplab][hv][lu]->SetLineColor(col[iplab]);
	tl_t[iplab][hv][lu]->SetLineWidth(2);
	tl_t[iplab][hv][lu]->Draw();
	tl_u[iplab][hv][lu]->SetLineColor(col[iplab]);
	tl_u[iplab][hv][lu]->SetLineWidth(2);
	tl_u[iplab][hv][lu]->Draw();
      }
    }

    tb_t[iplab][0] = new TBox(-14.0, costh(tmin[iplab], plab[iplab]), tmax[iplab], costh(tmax[iplab], plab[iplab]));
    tb_t[iplab][1] = new TBox(tmin[iplab], -1.1, tmax[iplab], costh(tmax[iplab], plab[iplab]));
    for (int i=0; i<2; ++i) {
      tb_t[iplab][i]->SetFillStyle(3001);
      tb_t[iplab][i]->SetFillColor(col[iplab]);
    }
    tb_t[iplab][0]->Draw();
    tb_t[iplab][1]->Draw();

    tb_u[iplab][0] = new TBox(-14.0, costh(umin[iplab], plab[iplab]), umax[iplab], costh(umax[iplab], plab[iplab]));
    tb_u[iplab][1] = new TBox(umin[iplab], -1.1, umax[iplab], costh(umax[iplab], plab[iplab]));
    for (int i=0; i<2; ++i) {
      tb_u[iplab][i]->SetFillStyle(3003);
      tb_u[iplab][i]->SetFillColor(col[iplab]);
    }
    tb_u[iplab][0]->Draw();
    tb_u[iplab][1]->Draw();

    tt[0][iplab]->DrawLatex(0.15,0.45,iplab==0?Form("p^{LAB}_{#bar{p}} = %5.3f GeV/c",plab[iplab]):Form("p^{LAB}_{#bar{p}} = %3.1f GeV/c",plab[iplab]));

    tt[1][iplab]->DrawLatex(iplab==0?0.4:(iplab==1?0.35:0.25),iplab==0?0.3:(iplab==1?0.22:0.20),"|u|<<Q^{2}");
    tt[4][iplab]->DrawLatex(iplab==0?0.36:(iplab==1?0.31:0.21),iplab==0?0.25:(iplab==1?0.17:0.15),Form("%3.2f<|u|<%3.2f",tmin[iplab],tmax[iplab]));

    tt[2][iplab]->DrawLatex(iplab==0?0.35:(iplab==1?0.4:0.45),iplab==0?0.71:(iplab==1?0.82:0.78),"|t|<<Q^{2}");
    tt[5][iplab]->DrawLatex(iplab==0?0.31:(iplab==1?0.36:0.41),iplab==0?0.67:(iplab==1?0.78:0.74),Form("%3.2f<|t|<%3.2f",tmin[iplab],tmax[iplab]));

    tc[iplab]->Update();

    tc[iplab]->Print(Form("%s/figs/2015.09.15/validity_iplab%d.png",bdir, iplab));
    gSystem->Exec(Form("convert %s/figs/2015.09.15/validity_iplab%d.png %s/figs/2015.09.15/validity_iplab%d.pdf",bdir, iplab, bdir, iplab));
  }

}

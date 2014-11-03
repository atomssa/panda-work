double rms_diff(TH1* ho, TH1* hn) {
  if (ho->GetRMS()!=0) {
    return fabs(ho->GetRMS() - hn->GetRMS())/ho->GetRMS();
  } else if (hn->GetRMS()!=0) {
    return fabs(ho->GetRMS() - hn->GetRMS())/hn->GetRMS();
  } else {
    return 0;
  }
}

void draw(TH1* ho, TH1* hn) {

  cout << "bad: " << ho->GetName() << " rmsOld= " << ho->GetRMS() << " rmsNew= " << hn->GetRMS()
       << " rmsdiff= " << rms_diff(ho,hn) << " entOld= " << ho->GetEntries() << " entNew= " << hn->GetEntries() << endl;
  TCanvas *tc = new TCanvas(Form("tc_%s",ho->GetName()),ho->GetTitle());
  tc->Divide(2);
  tc->cd(1);
  ho->Draw();
  tc->cd(2);

  hn->Draw();
  tc->Update();
  //string dummy;
  //std::getline(std::cin,dummy);
  //if (dummy == "q") {
  //	exit(0);
  //}
}

void comp_th2(TFile *fo, TFile *fn, const char* name) {
  //cout << "comp_th2 " << name;
  TH2F *ho = (TH2F*) fo->Get(name);
  //cout << " old ok" << ho;
  TH2F *hn = (TH2F*) fn->Get(name);
  //cout << " new ok" << hn << endl;
  //if ( rms_diff(ho,hn) > 0.01 || ho->GetEntries()==0) draw(ho,hn);
  if ( rms_diff(ho,hn) > 0.01) draw(ho,hn);
}

void comp_th1(TFile *fo, TFile *fn, const char* name) {
  //cout << "comp_th1 " << name;
  TH1F *ho = (TH1F*) fo->Get(name);
  //cout << " old ok" << ho;
  TH1F *hn = (TH1F*) fn->Get(name);
  //cout << " new ok" << hn << endl;
  //if ( rms_diff(ho,hn) > 0.01 || ho->GetEntries()==0) draw(ho,hn);
  if ( rms_diff(ho,hn) > 0.01 ) draw(ho,hn);
}

void comp () {

  //TFile *fo= new TFile("hists/ana_bg_brem_0.root");
  //TFile *fn= new TFile("test/ana_bg_brem_0.root");

  TFile *fo= new TFile("hists/no_mcut/ana_jpsi_brem.root");
  TFile *fn= new TFile("test/ana_jpsi_brem.root");

  TIter nextkey(fo->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
    //cout << "name = " << key->GetName() << " title = " << key->GetTitle() << " classname = " << key->GetClassName() << endl;
    if (strcmp(key->GetClassName(),"TH1F")==0) {
      comp_th1(fo,fn,key->GetName());
    } else if (strcmp(key->GetClassName(),"TH2F")==0) {
      comp_th2(fo,fn,key->GetName());
    }
  }
}

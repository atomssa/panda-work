#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"

void SaveAndUpdateHisto(TH1D* currenthisto, TFile& storagefile)
{
	if (storagefile.Get(currenthisto->GetName()) != 0) {
		currenthisto->Add((TH1D*)storagefile.Get(currenthisto->GetName()));
	}
	cout << currenthisto->GetName() << ": " << currenthisto->GetEntries() << endl;
	currenthisto->Write();
}

void run_ana_mertens_evt7(TString fname, TString simfname, int nevts=0, TString outfilename="psiana.root")
{

TStopwatch timer;
timer.Start();

gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");rootlogon();

PndEventReader evr(fname);
if (nevts==0) nevts=evr.GetEntries();

// Event Statistics Histos
TH1D* hkpperevent = new TH1D("hkpperevent","K+ per Event", 20, 0, 20);
TH1D* hkmperevent = new TH1D("hkmperevent","K- per Event", 20, 0, 20);
TH1D* hpipperevent = new TH1D("hpipperevent","Pi+ per Event", 20, 0, 20);
TH1D* hpimperevent = new TH1D("hpimperevent","Pi- per Event", 20, 0, 20);
TH1D* hniceevents = new TH1D("hniceevents", "Nice Events", 10, 0, 10);

// D Histos
TH1D* hdprawmass = new TH1D("hdprawmass","D+ Raw Mass",2000,0,4);
TH1D* hdmrawmass = new TH1D("hdmrawmass","D- Raw Mass",2000,0,4);
TH1D* hdpmselmass = new TH1D("hdpmselmass","D+ Selected Mass",2000,0,4);
TH1D* hdmmselmass = new TH1D("hdmmselmass","D- Selected Mass",2000,0,4);
TH1D* hdpvtxacceptmass = new TH1D("hdpvtxacceptmass","D+ Vertex Fit Accepted Mass",2000,0,4);
TH1D* hdpvtxrejectmass = new TH1D("hdpvtxrejectmass","D+ Vertex Fit Rejected Mass",2000,0,4);
hdpvtxrejectmass->SetLineColor(kRed);
TH1D* hdmvtxacceptmass = new TH1D("hdmvtxacceptmass","D- Vertex Fit Accepted Mass",2000,0,4);
TH1D* hdmvtxrejectmass = new TH1D("hdmvtxrejectmass","D- Vertex Fit Rejected Mass",2000,0,4);
hdmvtxrejectmass->SetLineColor(kRed);
TH1D* hdpvtxchi2 = new TH1D("hdpvtxchi2","D+ Vertex Fit Chi2",400,0,200);
TH1D* hdmvtxchi2 = new TH1D("hdmvtxchi2","D- Vertex Fit Chi2",400,0,200);

// Raw Psi Histos
TH1D* hrawpsirawmass = new TH1D("hrawpsirawmass","Psi Raw Mass (Raw D)",4000,0,8);
TH1D* hrawpsirawpt = new TH1D("hrawpsirawpt","Psi Raw Pt (Raw D)",4000,0,8);
TH1D* hrawpsirawpz = new TH1D("hrawpsirawpz","Psi Raw Pz (Raw D)",6000,0,12);
TH1D* hrawpsirawdpt = new TH1D("hrawpsirawdpt","Psi Raw D Pt (Raw D)",4000,0,8);
TH1D* hrawpsirawdpz = new TH1D("hrawpsirawdpz","Psi Raw D Pz (Raw D)",6000,0,12);

TH1D* hrawpsimselmass = new TH1D("hrawpsimselmass","Psi Raw Mass (Mass Selected D)",4000,0,8);
TH1D* hrawpsimselpt = new TH1D("hrawpsimselpt","Psi Raw Pt (Mass Selected D)",4000,0,8);
TH1D* hrawpsimselpz = new TH1D("hrawpsimselpz","Psi Raw Pz (Mass Selected D)",6000,0,12);
TH1D* hrawpsimseldpt = new TH1D("hrawpsimseldpt","Psi Raw D Pt (Mass Selected D)",4000,0,8);
TH1D* hrawpsimseldpz = new TH1D("hrawpsimseldpz","Psi Raw D Pz (Mass Selected D)",6000,0,12);

TH1D* hrawpsivtxmass = new TH1D("hrawpsivtxmass","Psi Raw Mass (Vertex Fit D)",4000,0,8);
TH1D* hrawpsivtxpt = new TH1D("hrawpsivtxpt","Psi Raw Pt (Vertex Fit D)",4000,0,8);
TH1D* hrawpsivtxpz = new TH1D("hrawpsivtxpz","Psi Raw Pz (Vertex Fit D)",6000,0,12);
TH1D* hrawpsivtxdpt = new TH1D("hrawpsivtxdpt","Psi Raw D Pt (Vertex Fit D)",4000,0,8);
TH1D* hrawpsivtxdpz = new TH1D("hrawpsivtxdpz","Psi Raw D Pz (Vertex Fit D)",6000,0,12);

TH1D* hdpvertex = new TH1D("hdpvertex","D+ Vertex Resolution",800,0,8);
TH1D* hdpvertexfit = new TH1D("hdpvertexfit","D+ Vertex Resolution (After Vertex Fit)",800,0,8);
TH1D* hdmvertex = new TH1D("hdmvertex","D- Vertex Resolution",800,0,8);
TH1D* hdmvertexfit = new TH1D("hdmvertexfit","D- Vertex Resolution (After Vertex Fit)",800,0,8);
TH1D* hfinalpsimass = new TH1D("hfinalpsimass","Psi Mass after all cuts (Vertex Fit D)",800,0,8);

TH1D* hdpvertexx = new TH1D("hdpvertexx","D+ Vertex Resolution X",800,-4,4);
TH1D* hdpvertexxfit = new TH1D("hdpvertexxfit","D+ Vertex Resolution X (After Vertex Fit)",800,-4,4);
TH1D* hdmvertexx = new TH1D("hdmvertexx","D- Vertex Resolution X",800,-4,4);
TH1D* hdmvertexxfit = new TH1D("hdmvertexxfit","D- Vertex Resolution X (After Vertex Fit)",800,-4,4);

TH1D* hdpvertexy = new TH1D("hdpvertexy","D+ Vertex Resolution Y",800,-4,4);
TH1D* hdpvertexyfit = new TH1D("hdpvertexyfit","D+ Vertex Resolution Y (After Vertex Fit)",800,-4,4);
TH1D* hdmvertexy = new TH1D("hdmvertexy","D- Vertex Resolution Y",800,-4,4);
TH1D* hdmvertexyfit = new TH1D("hdmvertexyfit","D- Vertex Resolution Y (After Vertex Fit)",800,-4,4);

TH1D* hdpvertexz = new TH1D("hdpvertexz","D+ Vertex Resolution Z",800,-4,4);
TH1D* hdpvertexzfit = new TH1D("hdpvertexzfit","D+ Vertex Resolution Z (After Vertex Fit)",800,-4,4);
TH1D* hdmvertexz = new TH1D("hdmvertexz","D- Vertex Resolution Z",800,-4,4);
TH1D* hdmvertexzfit = new TH1D("hdmvertexzfit","D- Vertex Resolution Z (After Vertex Fit)",800,-4,4);

TH1D* hdpvertexxmc = new TH1D("hdpvertexxmc","D+ MC Vertex X",800,-4,4);
TH1D* hdpvertexymc = new TH1D("hdpvertexymc","D+ MC Vertex Y",800,-4,4);
TH1D* hdpvertexzmc = new TH1D("hdpvertexzmc","D+ MC Vertex Z",800,-4,4);
TH1D* hdmvertexxmc = new TH1D("hdmvertexxmc","D- MC Vertex X",800,-4,4);
TH1D* hdmvertexymc = new TH1D("hdmvertexymc","D- MC Vertex Y",800,-4,4);
TH1D* hdmvertexzmc = new TH1D("hdmvertexzmc","D- MC Vertex Z",800,-4,4);

TH1D* hdpvertexxpoca = new TH1D("hdpvertexxpoca","D+ Poca Vertex X",800,-4,4);
TH1D* hdpvertexypoca = new TH1D("hdpvertexypoca","D+ Poca Vertex Y",800,-4,4);
TH1D* hdpvertexzpoca = new TH1D("hdpvertexzpoca","D+ Poca Vertex Z",800,-4,4);
TH1D* hdmvertexxpoca = new TH1D("hdmvertexxpoca","D- Poca Vertex X",800,-4,4);
TH1D* hdmvertexypoca = new TH1D("hdmvertexypoca","D- Poca Vertex Y",800,-4,4);
TH1D* hdmvertexzpoca = new TH1D("hdmvertexzpoca","D- Poca Vertex Z",800,-4,4);

TH1D* hdpvertexxpocareso = new TH1D("hdpvertexxpocareso","D+ Poca Vertex Resolution X",800,-4,4);
TH1D* hdpvertexypocareso = new TH1D("hdpvertexypocareso","D+ Poca Vertex Resolution Y",800,-4,4);
TH1D* hdpvertexzpocareso = new TH1D("hdpvertexzpocareso","D+ Poca Vertex Resolution Z",800,-4,4);
TH1D* hdmvertexxpocareso = new TH1D("hdmvertexxpocareso","D- Poca Vertex Resolution X",800,-4,4);
TH1D* hdmvertexypocareso = new TH1D("hdmvertexypocareso","D- Poca Vertex Resolution Y",800,-4,4);
TH1D* hdmvertexzpocareso = new TH1D("hdmvertexzpocareso","D- Poca Vertex Resolution Z",800,-4,4);

TH1D* hdpvertexxreco = new TH1D("hdpvertexxreco","D+ Reco Vertex X",800,-4,4);
TH1D* hdpvertexyreco = new TH1D("hdpvertexyreco","D+ Reco Vertex Y",800,-4,4);
TH1D* hdpvertexzreco = new TH1D("hdpvertexzreco","D+ Reco Vertex Z",800,-4,4);
TH1D* hdmvertexxreco = new TH1D("hdmvertexxreco","D- Reco Vertex X",800,-4,4);
TH1D* hdmvertexyreco = new TH1D("hdmvertexyreco","D- Reco Vertex Y",800,-4,4);
TH1D* hdmvertexzreco = new TH1D("hdmvertexzreco","D- Reco Vertex Z",800,-4,4);

TH1D* hdpdecaylength = new TH1D("hdpdecaylength","D+ Decay Length",8000,0,8);
TH1D* hdmdecaylength = new TH1D("hdmdecaylength","D- Decay Length",8000,0,8);

// Mass Selected Psi Histos
TH1D* hmselpsirawmass = new TH1D("hmselpsirawmass","Psi Mass Selected Mass (Raw D)",4000,0,8);
TH1D* hmselpsirawpt = new TH1D("hmselpsirawpt","Psi Mass Selected Pt (Raw D)",4000,0,8);
TH1D* hmselpsirawpz = new TH1D("hmselpsirawpz","Psi Mass Selected Pz (Raw D)",6000,0,12);
TH1D* hmselpsirawdpt = new TH1D("hmselpsirawdpt","Psi Mass Selected D Pt (Raw D)",4000,0,8);
TH1D* hmselpsirawdpz = new TH1D("hmselpsirawdpz","Psi Mass Selected D Pz (Raw D)",6000,0,12);
hmselpsirawmass->SetLineColor(kBlue);
hmselpsirawpt->SetLineColor(kBlue);
hmselpsirawpz->SetLineColor(kBlue);
hmselpsirawdpt->SetLineColor(kBlue);
hmselpsirawdpz->SetLineColor(kBlue);

TH1D* hmselpsimselmass = new TH1D("hmselpsimselmass","Psi Mass Selected Mass (Mass Selected D)",4000,0,8);
TH1D* hmselpsimselpt = new TH1D("hmselpsimselpt","Psi Mass Selected Pt (Mass Selected D)",4000,0,8);
TH1D* hmselpsimselpz = new TH1D("hmselpsimselpz","Psi Mass Selected Pz (Mass Selected D)",6000,0,12);
TH1D* hmselpsimseldpt = new TH1D("hmselpsimseldpt","Psi Mass Selected D Pt (Mass Selected D)",4000,0,8);
TH1D* hmselpsimseldpz = new TH1D("hmselpsimseldpz","Psi Mass Selected D Pz (Mass Selected D)",6000,0,12);
hmselpsimselmass->SetLineColor(kBlue);
hmselpsimselpt->SetLineColor(kBlue);
hmselpsimselpz->SetLineColor(kBlue);
hmselpsimseldpt->SetLineColor(kBlue);
hmselpsimseldpz->SetLineColor(kBlue);

TH1D* hmselpsivtxmass = new TH1D("hmselpsivtxmass","Psi Mass Selected Mass (Vertex Fit D)",4000,0,8);
TH1D* hmselpsivtxpt = new TH1D("hmselpsivtxpt","Psi Mass Selected Pt (Vertex Fit D)",4000,0,8);
TH1D* hmselpsivtxpz = new TH1D("hmselpsivtxpz","Psi Mass Selected Pz (Vertex Fit D)",6000,0,12);
TH1D* hmselpsivtxdpt = new TH1D("hmselpsivtxdpt","Psi Mass Selected D Pt (Vertex Fit D)",4000,0,8);
TH1D* hmselpsivtxdpz = new TH1D("hmselpsivtxdpz","Psi Mass Selected D Pz (Vertex Fit D)",6000,0,12);
hmselpsivtxmass->SetLineColor(kBlue);
hmselpsivtxpt->SetLineColor(kBlue);
hmselpsivtxpz->SetLineColor(kBlue);
hmselpsivtxdpt->SetLineColor(kBlue);
hmselpsivtxdpz->SetLineColor(kBlue);

// D pt Selected Psi Histos
TH1D* hdptselpsirawmass = new TH1D("hdptselpsirawmass","Psi pt (D) Selected Mass (Raw D)",4000,0,8);
TH1D* hdptselpsirawpt = new TH1D("hdptselpsirawpt","Psi pt (D) Pt (Raw D)",4000,0,8);
TH1D* hdptselpsirawpz = new TH1D("hdptselpsirawpz","Psi pt (D) Selected Pz (Raw D)",6000,0,12);
TH1D* hdptselpsirawdpt = new TH1D("hdptselpsirawdpt","Psi pt (D) Selected D Pt (Raw D)",4000,0,8);
TH1D* hdptselpsirawdpz = new TH1D("hdptselpsirawdpz","Psi pt (D) Selected D Pz (Raw D)",6000,0,12);

TH1D* hdptselpsimselmass = new TH1D("hdptselpsimselmass","Psi pt (D) Selected Mass (Mass Selected D)",4000,0,8);
TH1D* hdptselpsimselpt = new TH1D("hdptselpsimselpt","Psi pt (D) Selected Pt (Mass Selected D)",4000,0,8);
TH1D* hdptselpsimselpz = new TH1D("hdptselpsimselpz","Psi pt (D) Selected Pz (Mass Selected D)",6000,0,12);
TH1D* hdptselpsimseldpt = new TH1D("hdptselpsimseldpt","Psi pt (D) Selected D Pt (Mass Selected D)",4000,0,8);
TH1D* hdptselpsimseldpz = new TH1D("hdptselpsimseldpz","Psi pt (D) Selected D Pz (Mass Selected D)",6000,0,12);

TH1D* hdptselpsivtxmass = new TH1D("hdptselpsivtxmass","Psi pt (D) Selected Mass (Vertex Fit D)",4000,0,8);
TH1D* hdptselpsivtxpt = new TH1D("hdptselpsivtxpt","Psi pt (D) Selected Pt (Vertex Fit D)",4000,0,8);
TH1D* hdptselpsivtxpz = new TH1D("hdptselpsivtxpz","Psi pt (D) Selected Pz (Vertex Fit D)",6000,0,12);
TH1D* hdptselpsivtxdpt = new TH1D("hdptselpsivtxdpt","Psi pt (D) Selected D Pt (Vertex Fit D)",4000,0,8);
TH1D* hdptselpsivtxdpz = new TH1D("hdptselpsivtxdpz","Psi pt (D) Selected D Pz (Vertex Fit D)",6000,0,12);

// Psi pt Selected Psi Histos
TH1D* hpptselpsirawmass = new TH1D("hpptselpsirawmass","Psi pt (Psi) Selected Mass (Raw D)",4000,0,8);
TH1D* hpptselpsirawpt = new TH1D("hpptselpsirawpt","Psi pt (Psi) Pt (Raw D)",4000,0,8);
TH1D* hpptselpsirawpz = new TH1D("hpptselpsirawpz","Psi pt (Psi) Selected Pz (Raw D)",6000,0,12);
TH1D* hpptselpsirawdpt = new TH1D("hpptselpsirawdpt","Psi pt (Psi) Selected D Pt (Raw D)",4000,0,8);
TH1D* hpptselpsirawdpz = new TH1D("hpptselpsirawdpz","Psi pt (Psi) Selected D Pz (Raw D)",6000,0,12);

TH1D* hpptselpsimselmass = new TH1D("hpptselpsimselmass","Psi pt (Psi) Selected Mass (Mass Selected D)",4000,0,8);
TH1D* hpptselpsimselpt = new TH1D("hpptselpsimselpt","Psi pt (Psi) Selected Pt (Mass Selected D)",4000,0,8);
TH1D* hpptselpsimselpz = new TH1D("hpptselpsimselpz","Psi pt (Psi) Selected Pz (Mass Selected D)",6000,0,12);
TH1D* hpptselpsimseldpt = new TH1D("hpptselpsimseldpt","Psi pt (Psi) Selected D Pt (Mass Selected D)",4000,0,8);
TH1D* hpptselpsimseldpz = new TH1D("hpptselpsimseldpz","Psi pt (Psi) Selected D Pz (Mass Selected D)",6000,0,12);

TH1D* hpptselpsivtxmass = new TH1D("hpptselpsivtxmass","Psi pt (Psi) Selected Mass (Vertex Fit D)",4000,0,8);
TH1D* hpptselpsivtxpt = new TH1D("hpptselpsivtxpt","Psi pt (Psi) Selected Pt (Vertex Fit D)",4000,0,8);
TH1D* hpptselpsivtxpz = new TH1D("hpptselpsivtxpz","Psi pt (Psi) Selected Pz (Vertex Fit D)",6000,0,12);
TH1D* hpptselpsivtxdpt = new TH1D("hpptselpsivtxdpt","Psi pt (Psi) Selected D Pt (Vertex Fit D)",4000,0,8);
TH1D* hpptselpsivtxdpz = new TH1D("hpptselpsivtxdpz","Psi pt (Psi) Selected D Pz (Vertex Fit D)",6000,0,12);

// D and Psi pt Selected Psi Histos
TH1D* hptselpsirawmass = new TH1D("hptselpsirawmass","Psi pt (Both) Selected Mass (Raw D)",4000,0,8);
TH1D* hptselpsirawpt = new TH1D("hptselpsirawpt","Psi pt (Both) Pt (Raw D)",4000,0,8);
TH1D* hptselpsirawpz = new TH1D("hptselpsirawpz","Psi pt (Both) Selected Pz (Raw D)",6000,0,12);
TH1D* hptselpsirawdpt = new TH1D("hptselpsirawdpt","Psi pt (Both) Selected D Pt (Raw D)",4000,0,8);
TH1D* hptselpsirawdpz = new TH1D("hptselpsirawdpz","Psi pt (Both) Selected D Pz (Raw D)",6000,0,12);

TH1D* hptselpsimselmass = new TH1D("hptselpsimselmass","Psi pt (Both) Selected Mass (Mass Selected D)",4000,0,8);
TH1D* hptselpsimselpt = new TH1D("hptselpsimselpt","Psi pt (Both) Selected Pt (Mass Selected D)",4000,0,8);
TH1D* hptselpsimselpz = new TH1D("hptselpsimselpz","Psi pt (Both) Selected Pz (Mass Selected D)",6000,0,12);
TH1D* hptselpsimseldpt = new TH1D("hptselpsimseldpt","Psi pt (Both) Selected D Pt (Mass Selected D)",4000,0,8);
TH1D* hptselpsimseldpz = new TH1D("hptselpsimseldpz","Psi pt (Both) Selected D Pz (Mass Selected D)",6000,0,12);

TH1D* hptselpsivtxmass = new TH1D("hptselpsivtxmass","Psi pt (Both) Selected Mass (Vertex Fit D)",4000,0,8);
TH1D* hptselpsivtxpt = new TH1D("hptselpsivtxpt","Psi pt (Both) Selected Pt (Vertex Fit D)",4000,0,8);
TH1D* hptselpsivtxpz = new TH1D("hptselpsivtxpz","Psi pt (Both) Selected Pz (Vertex Fit D)",6000,0,12);
TH1D* hptselpsivtxdpt = new TH1D("hptselpsivtxdpt","Psi pt (Both) Selected D Pt (Vertex Fit D)",4000,0,8);
TH1D* hptselpsivtxdpz = new TH1D("hptselpsivtxdpz","Psi pt (Both) Selected D Pz (Vertex Fit D)",6000,0,12);

// Vertex Fit Psi Histos (only taken for preselected Psi)
TH1D* hvtxpsichi2 = new TH1D("hvtxpsichi2","Psi Vertex Fit Chi2",400,0,200);
TH1D* hvtxpsiacceptmass = new TH1D("hvtxpsiacceptmass","Psi Vertex Fit Accepted Mass",4000,0,8);
TH1D* hvtxpsiacceptpt = new TH1D("hvtxpsiacceptpt","Psi Vertex Fit Accepted Pt",4000,0,8);
TH1D* hvtxpsiacceptpz = new TH1D("hvtxpsiacceptpz","Psi Vertex Fit Accepted Pz",6000,0,12);
TH1D* hvtxpsiacceptdpt = new TH1D("hvtxpsiacceptdpt","Psi Vertex Fit Accepted D Pt",4000,0,8);
TH1D* hvtxpsiacceptdpz = new TH1D("hvtxpsiacceptdpz","Psi Vertex Fit Accepted D Pz",6000,0,12);
TH1D* hvtxpsirejectmass = new TH1D("hvtxpsirejectmass","Psi Vertex Fit Rejected Mass",4000,0,8);
TH1D* hvtxpsirejectpt = new TH1D("hvtxpsirejectpt","Psi Vertex Fit Rejected Pt",4000,0,8);
TH1D* hvtxpsirejectpz = new TH1D("hvtxpsirejectpz","Psi Vertex Fit Rejected Pz",6000,0,12);
TH1D* hvtxpsirejectdpt = new TH1D("hvtxpsirejectdpt","Psi Vertex Fit Rejected D Pt",4000,0,8);
TH1D* hvtxpsirejectdpz = new TH1D("hvtxpsirejectdpz","Psi Vertex Fit Rejected D Pz",6000,0,12);
hvtxpsirejectmass->SetLineColor(kRed);
hvtxpsirejectpt->SetLineColor(kRed);
hvtxpsirejectpz->SetLineColor(kRed);
hvtxpsirejectdpt->SetLineColor(kRed);
hvtxpsirejectdpz->SetLineColor(kRed);

TH1D* hptselvtxpsimass = new TH1D("hptselvtxpsimass","Psi Vertex Fit pt (Both) Selected Mass",4000,0,8);
TH1D* hptselvtxpsipt = new TH1D("hptselvtxpsipt","Psi Vertex Fit pt (Both) Selected Pt",4000,0,8);
TH1D* hptselvtxpsipz = new TH1D("hptselvtxpsipz","Psi Vertex Fit pt (Both) Selected Pz",6000,0,12);
TH1D* hptselvtxpsidpt = new TH1D("hptselvtxpsidpt","Psi Vertex Fit pt (Both) Selected D Pt",4000,0,8);
TH1D* hptselvtxpsidpz = new TH1D("hptselvtxpsidpz","Psi Vertex Fit pt (Both) Selected D Pz",6000,0,12);

TH1D* dpmass = new TH1D("dpmass","D+ Raw",2000,0,4);
TH1D* dmmass = new TH1D("dmmass","D- Raw",2000,0,4);
TH1D* dp2mass = new TH1D("dp2mass","D+ after Selection",2000,0,4);
TH1D* dm2mass = new TH1D("dm2mass","D- after Selection",2000,0,4);
TH1D* dp3mass = new TH1D("dp3mass","D+ after Vertex Fit",2000,0,4);
TH1D* dm3mass = new TH1D("dm3mass","D- after Vertex Fit",2000,0,4);
TH1D* dp3massrej = new TH1D("dp3massrej","D+ after Vertex Fit",2000,0,4);
dp3massrej->SetLineColor(kRed);
TH1D* dm3massrej = new TH1D("dm3massrej","D- after Vertex Fit",2000,0,4);
dm3massrej->SetLineColor(kRed);

TH1D* dptfit = new TH1D("dptfit","Pt Difference",2000,-2,2);
TH1D* dptsel = new TH1D("dptsel","Pt Difference",2000,-2,2);
TH1D* psipt = new TH1D("psipt","Psi Pt",2000,-2,2);

TH1D* dpchi2 = new TH1D("dpchi2","D+ Chi2",400,0,200);
TH1D* dmchi2 = new TH1D("dmchi2","D- Chi2",400,0,200);
TH1D* psichi2 = new TH1D("psichi2","Psi Chi2",400,0,200);

TH1D* psimass = new TH1D("psimass","Psi Raw",4000,0,8);
TH1D* psimass2 = new TH1D("psimass2","Psi with selected D",4000,0,8);
TH1D* psimass3 = new TH1D("psimass3","Psi with Vertex Fit D",4000,0,8);
TH1D* psiptsel = new TH1D("psiptsel","Psi with small D Pt",4000,0,8);
TH1D* psiptselfine = new TH1D("psiptselfine","Psi with small Pt",4000,0,8);

TH1D* ppmass = new TH1D("ppmass","4C Fit pbarp input mass",8000,0,16);
TH1D* ppmass2 = new TH1D("ppmass2","4C Fit pbarp fit mass",8000,0,16);

// particle lists
TCandList pp;
TCandList pipbase, pimbase, kpbase, kmbase;
TCandList pip, pim, kp, km;
TCandList dpraw, dmraw;
TCandList dpmsel, dmmsel;
TCandList dpvtx, dmvtx;
TCandList rawpsiraw, rawpsimsel, rawpsivtx;
TCandList mselpsiraw, mselpsimsel, mselpsivtx;
TCandList dptselpsiraw, dptselpsimsel, dptselpsivtx;
TCandList pptselpsiraw, pptselpsimsel, pptselpsivtx;
TCandList ptselpsiraw, ptselpsimsel, ptselpsivtx;
TCandList vtxpsi, ptselvtxpsi, finalpsi;
TCandList mct;
TCandList dp,dm,pp,psi3770,psiraw,psi3770fit;

//TLorentzVector ini(0,0,6.23164,7.24015);
TPidMassSelector *jpsiMSel = new TPidMassSelector("jpsiMSel" , 3.096 , 0.3);
TPidMassSelector *dMSel = new TPidMassSelector("dMSel", TRho::Instance()->GetPDG()->GetParticle(411)->Mass(), 0.9);
TPidMassSelector *dMSelfine = new TPidMassSelector("dMSelfine", TRho::Instance()->GetPDG()->GetParticle(411)->Mass(), 0.5);
TPidMassSelector* psiMSel = new TPidMassSelector("psiMSel", TRho::Instance()->GetPDG()->GetParticle(40443)->Mass(), 0.9);

Double_t dvtxfitchi2limit = 18;
Double_t ptlimit = 0.5;
Double_t inipz = 6.5788;
Double_t pzlimit = 0.5;

int i=0,j=0;

TFile outputstorage(outfilename, "UPDATE");
TFile mcfile(simfname, "READ");
TTree* mctree = (TTree*)mcfile.Get("cbmsim");
TClonesArray* mcarray = new TClonesArray("PndMCTrack");
mctree->SetBranchAddress("MCTrack", &mcarray);

TVector3 dpstartvertex;
TVector3 dmstartvertex;

//PndAnalysis* theAnalysis = new PndAnalysis("BarrelGenTrack");

// **** loop over all _events_
//
while (evr.GetEvent() && ++i<nevts) {

	// sorry, main loop too long. will not indent code

if (!(i%100)) cout <<"evt "<<i<<endl;

mctree->GetEntry(i);

pp.Cleanup();
pip.Cleanup();
pim.Cleanup();
kp.Cleanup();
km.Cleanup();
dpraw.Cleanup();
dmraw.Cleanup();
dpmsel.Cleanup();
dmmsel.Cleanup();
dpvtx.Cleanup();
dmvtx.Cleanup();
rawpsiraw.Cleanup();
rawpsimsel.Cleanup();
rawpsivtx.Cleanup();
mselpsiraw.Cleanup();
mselpsimsel.Cleanup();
mselpsivtx.Cleanup();
dptselpsiraw.Cleanup();
dptselpsimsel.Cleanup();
dptselpsivtx.Cleanup();
pptselpsiraw.Cleanup();
pptselpsimsel.Cleanup();
pptselpsivtx.Cleanup();
ptselpsiraw.Cleanup();
ptselpsimsel.Cleanup();
ptselpsivtx.Cleanup();
ptselvtxpsi.Cleanup();
vtxpsi.Cleanup();
finalpsi.Cleanup();

evr.FillList(pipbase,"PionVeryLoosePlus");
evr.FillList(pimbase,"PionVeryLooseMinus");
evr.FillList(kpbase, "KaonVeryLoosePlus");
evr.FillList(kmbase, "KaonVeryLooseMinus");

// use Monte Carlo PID information to clean particle lists
for (j=0; j<kpbase.GetLength(); ++j) {
	int mcindex = kpbase[j].GetMcIdx();
	if (mcindex >= 0) PndMCTrack *mctrack = (PndMCTrack*)mcarray->At(mcindex);
	if ((mctrack != 0) && (mctrack->GetPdgCode() == 321)) kp.Add(kpbase[j]);
}
for (j=0; j<kmbase.GetLength(); ++j) {
	int mcindex = kmbase[j].GetMcIdx();
	if (mcindex >= 0) PndMCTrack *mctrack = (PndMCTrack*)mcarray->At(mcindex);
	if ((mctrack != 0) && (mctrack->GetPdgCode() == -321)) km.Add(kmbase[j]);
}

// store Monte Carlo Start Vertex information for later comparison
for (j=0; j<pipbase.GetLength(); ++j) {
	int mcindex = pipbase[j].GetMcIdx();
	if (mcindex >= 0) PndMCTrack *mctrack = (PndMCTrack*)mcarray->At(mcindex);
	if ((mctrack != 0) && (mctrack->GetPdgCode() == 211)) {
		Double_t x;
		Double_t y;
		Double_t z;
		std::stringstream strstrx;
		std::stringstream strstry;
		std::stringstream strstrz;
		strstrx << mctrack->GetStartVertex().x();
		strstrx >> x;
		strstry << mctrack->GetStartVertex().y();
		strstry >> y;
		strstrz << mctrack->GetStartVertex().z();
		strstrz >> z;
		dpstartvertex.SetXYZ(x, y, z);
		pip.Add(pipbase[j]);
		//dpstartvertex = mctrack->GetStartVertex();
	}
}
for (j=0; j<pimbase.GetLength(); ++j) {
	int mcindex = pimbase[j].GetMcIdx();
	if (mcindex >= 0) PndMCTrack *mctrack = (PndMCTrack*)mcarray->At(mcindex);
	if ((mctrack != 0) && (mctrack->GetPdgCode() == -211)) {
	Double_t x;
	Double_t y;
	Double_t z;
	std::stringstream strstrx;
	std::stringstream strstry;
	std::stringstream strstrz;
	strstrx << mctrack->GetStartVertex().x();
	strstrx >> x;
	strstry << mctrack->GetStartVertex().y();
	strstry >> y;
	strstrz << mctrack->GetStartVertex().z();
	strstrz >> z;
	dmstartvertex.SetXYZ(x, y, z);
	pim.Add(pimbase[j]);
	//dmstartvertex = mctrack->GetStartVertex();
	}
}

hkpperevent->Fill(kp.GetLength());
hkmperevent->Fill(km.GetLength());
hpipperevent->Fill(pip.GetLength());
hpimperevent->Fill(pim.GetLength());
if ((kp.GetLength() == 1) && (km.GetLength() == 1) && (pip.GetLength() == 2) && (pim.GetLength() == 2))		hniceevents->Fill(0);
if ((kp.GetLength() >= 1) && (km.GetLength() >= 1) && (pip.GetLength() >= 2) && (pim.GetLength() >= 2))		hniceevents->Fill(1);
if ((kp.GetLength() + pim.GetLength() == 3) && (km.GetLength() + pip.GetLength() == 3))				hniceevents->Fill(2);
if ((kp.GetLength() + pim.GetLength() == 3) || (km.GetLength() + pip.GetLength() == 3))				hniceevents->Fill(3);

if ((kp.GetLength()  > 0) && (km.GetLength()  > 0) && (pip.GetLength()  > 1) && (pim.GetLength()  > 1)) {
	dpraw.Combine(km,pip,pip);
	//cout << "D+ TCandList entries: " << dpraw.GetLength() << endl;
	dmraw.Combine(kp,pim,pim);
	//cout << "D- TCandList entries: " << dmraw.GetLength() << endl;
	dpmsel.Select(dpraw, dMSel);
	dmmsel.Select(dmraw, dMSel);
}

for (j=0;j<dpraw.GetLength();++j) {
	hdprawmass->Fill(dpraw[j].M());
}
for (j=0;j<dmraw.GetLength();++j) {
	hdmrawmass->Fill(dmraw[j].M());
}
for (j=0;j<dpmsel.GetLength();++j) {
	hdpmselmass->Fill(dpmsel[j].M());
}
for (j=0;j<dmmsel.GetLength();++j) {
	hdmmselmass->Fill(dmmsel[j].M());
}

//POCA test

for (j=0; j < dpmsel.GetLength(); ++j) {
	PndVtxPoca vtxfinder(dpmsel[j]);
	TVector3 pocavertex(0,0,0);
	double distval=-4444;
	distval = vtxfinder.GetPocaVtx(pocavertex);
	//reco vertex position
	hdpvertexxpoca->Fill(pocavertex.x());
	hdpvertexypoca->Fill(pocavertex.y());
	hdpvertexzpoca->Fill(pocavertex.z());
	//poca vertex resolution
	TVector3 vertexdisp = dpstartvertex - pocavertex;
	hdpvertexxpocareso->Fill(vertexdisp.x());
	hdpvertexypocareso->Fill(vertexdisp.y());
	hdpvertexzpocareso->Fill(vertexdisp.z());
}
for (j=0; j < dmmsel.GetLength(); ++j) {
	PndVtxPoca vtxfinder(dmmsel[j]);
	TVector3 pocavertex(0,0,0);
	double distval=-4444;
	distval = vtxfinder.GetPocaVtx(pocavertex);
	//reco vertex position
	hdmvertexxpoca->Fill(pocavertex.x());
	hdmvertexypoca->Fill(pocavertex.y());
	hdmvertexzpoca->Fill(pocavertex.z());
	//poca vertex resolution
	TVector3 vertexdisp = dmstartvertex - pocavertex;
	hdmvertexxpocareso->Fill(vertexdisp.x());
	hdmvertexypocareso->Fill(vertexdisp.y());
	hdmvertexzpocareso->Fill(vertexdisp.z());
}


//vertex fit for D+ and D-

dpvtx.Cleanup();
dmvtx.Cleanup();

int bestfitindex;
Double_t bestfitchi2 = dvtxfitchi2limit;
Double_t bestfitmass = 0;
for (j=0; j < dpmsel.GetLength(); ++j) {
	PndKinVtxFitter vtxfitter(dpmsel[j]);
	//vtxfitter.AddMassConstraint(TRho::Instance()->GetPDG()->GetParticle(411)->Mass());
	vtxfitter.Fit();
	TCandidate fitcand=*(vtxfitter.FittedCand(dpmsel[j]));
	TVector3 fitvertex=fitcand.Pos();
	hdpvtxchi2->Fill(vtxfitter.GlobalChi2());
	if (vtxfitter.GlobalChi2()<bestfitchi2) {
		if (bestfitchi2 < dvtxfitchi2limit) {
			hdpvtxrejectmass->Fill(bestfitmass);
		}
		bestfitchi2 = vtxfitter.GlobalChi2();
		bestfitindex = j;
		bestfitmass = fitcand.M();
	} else {
		hdpvtxrejectmass->Fill(fitcand.M());
	}
}
//process best candidate
if (bestfitchi2 < dvtxfitchi2limit) {
	PndKinVtxFitter vtxfitter(dpmsel[bestfitindex]);
	//vtxfitter.AddMassConstraint(TRho::Instance()->GetPDG()->GetParticle(411)->Mass());
	vtxfitter.Fit();
	TCandidate fitcand=*(vtxfitter.FittedCand(dpmsel[bestfitindex]));
	TVector3 fitvertex=fitcand.Pos();
	dpvtx.Add(fitcand);
	//distance between fitted vertex and MC vertex
	TVector3 vertexdisp = dpstartvertex - fitvertex;
	hdpvertexfit->Fill(vertexdisp.Mag());
	hdpvertexxfit->Fill(vertexdisp.x());
	hdpvertexyfit->Fill(vertexdisp.y());
	hdpvertexzfit->Fill(vertexdisp.z());
	//distance between initial vertex and MC vertex (initial vertex always 0,0,0)
	vertexdisp = dpstartvertex - dpmsel[bestfitindex].Pos();
	hdpvertex->Fill(vertexdisp.Mag());
	hdpvertexx->Fill(vertexdisp.x());
	hdpvertexy->Fill(vertexdisp.y());
	hdpvertexz->Fill(vertexdisp.z());
	//MC vertex position
	hdpvertexxmc->Fill(dpstartvertex.x());
	hdpvertexymc->Fill(dpstartvertex.y());
	hdpvertexzmc->Fill(dpstartvertex.z());
	//reco vertex position
	hdpvertexxreco->Fill(fitvertex.x());
	hdpvertexyreco->Fill(fitvertex.y());
	hdpvertexzreco->Fill(fitvertex.z());
	//D+ decay length
	hdpdecaylength->Fill(fitcand.M()*fitvertex.Mag()/(fitcand.P()*30));
}

bestfitchi2 = dvtxfitchi2limit;
for (j=0; j<dmmsel.GetLength(); ++j) {
	PndKinVtxFitter vtxfitter(dmmsel[j]);
	//vtxfitter.AddMassConstraint(TRho::Instance()->GetPDG()->GetParticle(411)->Mass());
	vtxfitter.Fit();
	TCandidate fitcand=*(vtxfitter.FittedCand(dmmsel[j]));
	TVector3 fitvertex=fitcand.Pos();
	hdmvtxchi2->Fill(vtxfitter.GlobalChi2());
	if (vtxfitter.GlobalChi2()<bestfitchi2) {
		if (bestfitchi2 < dvtxfitchi2limit) {
			hdmvtxrejectmass->Fill(bestfitmass);
		}
		bestfitchi2 = vtxfitter.GlobalChi2();
		bestfitindex = j;
		bestfitmass = fitcand.M();
	} else {
		hdmvtxrejectmass->Fill(fitcand.M());
	}
}
//process best candidate
if (bestfitchi2 < dvtxfitchi2limit) {
	PndKinVtxFitter vtxfitter(dmmsel[bestfitindex]);
	//vtxfitter.AddMassConstraint(TRho::Instance()->GetPDG()->GetParticle(411)->Mass());
	vtxfitter.Fit();
	TCandidate fitcand=*(vtxfitter.FittedCand(dmmsel[bestfitindex]));
	TVector3 fitvertex=fitcand.Pos();
	dmvtx.Add(fitcand);
	//distance between fitted vertex and MC vertex
	TVector3 vertexdisp = dmstartvertex - fitvertex;
	hdmvertexfit->Fill(vertexdisp.Mag());
	hdmvertexxfit->Fill(vertexdisp.x());
	hdmvertexyfit->Fill(vertexdisp.y());
	hdmvertexzfit->Fill(vertexdisp.z());
	vertexdisp = dmstartvertex - dmmsel[bestfitindex].Pos();
	hdmvertex->Fill(vertexdisp.Mag());
	hdmvertexx->Fill(vertexdisp.x());
	hdmvertexy->Fill(vertexdisp.y());
	hdmvertexz->Fill(vertexdisp.z());
	hdmvertexxmc->Fill(dmstartvertex.x());
	hdmvertexymc->Fill(dmstartvertex.y());
	hdmvertexzmc->Fill(dmstartvertex.z());
	hdmvertexxreco->Fill(fitvertex.x());
	hdmvertexyreco->Fill(fitvertex.y());
	hdmvertexzreco->Fill(fitvertex.z());
	hdmdecaylength->Fill(fitcand.M()*fitvertex.Mag()/(fitcand.P()*30));
}

dpvtx.Select(dMSelfine);
dmvtx.Select(dMSelfine);

for (j=0;j<dpvtx.GetLength();++j) {
	hdpvtxacceptmass->Fill(dpvtx[j].M());
}
for (j=0;j<dmvtx.GetLength();++j) {
	hdmvtxacceptmass->Fill(dmvtx[j].M());
}

// Make RAW Psi and prepare lists to hold pt selected Psi

TLorentzVector sum;

dptselpsiraw.Cleanup();
pptselpsiraw.Cleanup();
ptselpsiraw.Cleanup();
dptselpsimsel.Cleanup();
pptselpsimsel.Cleanup();
ptselpsimsel.Cleanup();
dptselpsivtx.Cleanup();
pptselpsivtx.Cleanup();
ptselpsivtx.Cleanup();

// Psi from all D mesons
rawpsiraw.Combine(dpraw,dmraw);
for (j=0;j<rawpsiraw.GetLength();++j) {
	hrawpsirawmass->Fill(rawpsiraw[j].M());
	hrawpsirawpt->Fill(rawpsiraw[j].Pt());
	hrawpsirawpz->Fill(rawpsiraw[j].Pz());
	sum = rawpsiraw[j].Daughter(0)->P4() + rawpsiraw[j].Daughter(1)->P4();
	hrawpsirawdpt->Fill(sum.Pt());
	hrawpsirawdpz->Fill(sum.Pz());
	if (sum.Pt() < ptlimit) {
		dptselpsiraw.Add(rawpsiraw[j]);
	}
	if (rawpsiraw[j].Pt() < ptlimit) {
		pptselpsiraw.Add(rawpsiraw[j]);
	}
	if ((sum.Pt() < ptlimit)&&(rawpsiraw[j].Pt() < ptlimit)) {
		ptselpsiraw.Add(rawpsiraw[j]);
	}
}

// Psi from mass selected D mesons
rawpsimsel.Combine(dpmsel,dmmsel);
for (j=0;j<rawpsimsel.GetLength();++j) {
	hrawpsimselmass->Fill(rawpsimsel[j].M());
	hrawpsimselpt->Fill(rawpsimsel[j].Pt());
	hrawpsimselpz->Fill(rawpsimsel[j].Pz());
	sum = rawpsimsel[j].Daughter(0)->P4() + rawpsimsel[j].Daughter(1)->P4();
	hrawpsimseldpt->Fill(sum.Pt());
	hrawpsimseldpz->Fill(sum.Pz());
	if (sum.Pt() < ptlimit) {
		dptselpsimsel.Add(rawpsimsel[j]);
	}
	if (rawpsimsel[j].Pt() < ptlimit) {
		pptselpsimsel.Add(rawpsimsel[j]);
	}
	if ((sum.Pt() < ptlimit)&&(rawpsimsel[j].Pt() < ptlimit)) {
		ptselpsimsel.Add(rawpsimsel[j]);
	}
}

// Psi from vertex fitted D mesons
rawpsivtx.Combine(dpvtx,dmvtx);
for (j=0;j<rawpsivtx.GetLength();++j) {
	hrawpsivtxmass->Fill(rawpsivtx[j].M());
	hrawpsivtxpt->Fill(rawpsivtx[j].Pt());
	hrawpsivtxpz->Fill(rawpsivtx[j].Pz());
	sum = rawpsivtx[j].Daughter(0)->P4() + rawpsivtx[j].Daughter(1)->P4();
	hrawpsivtxdpt->Fill(sum.Pt());
	hrawpsivtxdpz->Fill(sum.Pz());
	if (sum.Pt() < ptlimit) {
		dptselpsivtx.Add(rawpsivtx[j]);
	}
	if (rawpsivtx[j].Pt() < ptlimit) {
		pptselpsivtx.Add(rawpsivtx[j]);
	}
	if ((sum.Pt() < ptlimit)&&(rawpsivtx[j].Pt() < ptlimit)) {
		ptselpsivtx.Add(rawpsivtx[j]);
	}
	if ((sum.Pt() < ptlimit)&&(rawpsivtx[j].Pt() < ptlimit) && (rawpsivtx[j].Pz() > (inipz - pzlimit)) && (rawpsivtx[j].Pz() < (inipz + 	pzlimit))) {
		finalpsi.Add(rawpsivtx[j]);
	}
}

for (j=0;j<finalpsi.GetLength();++j) {
	hfinalpsimass->Fill(finalpsi[j].M());
}

// make mass selected Psi from Psi from all D mesons
mselpsiraw.Select(rawpsiraw, psiMSel);
for (j=0;j<mselpsiraw.GetLength();++j) {
	hmselpsirawmass->Fill(mselpsiraw[j].M());
	hmselpsirawpt->Fill(mselpsiraw[j].Pt());
	hmselpsirawpz->Fill(mselpsiraw[j].Pz());
	sum = mselpsiraw[j].Daughter(0)->P4() + mselpsiraw[j].Daughter(1)->P4();
	hmselpsirawdpt->Fill(sum.Pt());
	hmselpsirawdpz->Fill(sum.Pz());
}

// make mass selected Psi from Psi from mass selected D mesons
mselpsimsel.Select(rawpsimsel, psiMSel);
for (j=0;j<mselpsimsel.GetLength();++j) {
	hmselpsimselmass->Fill(mselpsimsel[j].M());
	hmselpsimselpt->Fill(mselpsimsel[j].Pt());
	hmselpsimselpz->Fill(mselpsimsel[j].Pz());
	sum = mselpsimsel[j].Daughter(0)->P4() + mselpsimsel[j].Daughter(1)->P4();
	hmselpsimseldpt->Fill(sum.Pt());
	hmselpsimseldpz->Fill(sum.Pz());
}

// make mass selected Psi from Psi from vertex fitted D mesons
mselpsivtx.Select(rawpsivtx, psiMSel);
for (j=0;j<mselpsivtx.GetLength();++j) {
	hmselpsivtxmass->Fill(mselpsivtx[j].M());
	hmselpsivtxpt->Fill(mselpsivtx[j].Pt());
	hmselpsivtxpz->Fill(mselpsivtx[j].Pz());
	sum = mselpsivtx[j].Daughter(0)->P4() + mselpsivtx[j].Daughter(1)->P4();
	hmselpsivtxdpt->Fill(sum.Pt());
	hmselpsivtxdpz->Fill(sum.Pz());
}


// plot pt (D) selected Psi
for (j=0;j<dptselpsiraw.GetLength();++j) {
	hdptselpsirawmass->Fill(dptselpsiraw[j].M());
	hdptselpsirawpt->Fill(dptselpsiraw[j].Pt());
	hdptselpsirawpz->Fill(dptselpsiraw[j].Pz());
	sum = dptselpsiraw[j].Daughter(0)->P4() + dptselpsiraw[j].Daughter(1)->P4();
	hdptselpsirawdpt->Fill(sum.Pt());
	hdptselpsirawdpz->Fill(sum.Pz());
}
for (j=0;j<dptselpsimsel.GetLength();++j) {
	hdptselpsimselmass->Fill(dptselpsimsel[j].M());
	hdptselpsimselpt->Fill(dptselpsimsel[j].Pt());
	hdptselpsimselpz->Fill(dptselpsimsel[j].Pz());
	sum = dptselpsimsel[j].Daughter(0)->P4() + dptselpsimsel[j].Daughter(1)->P4();
	hdptselpsimseldpt->Fill(sum.Pt());
	hdptselpsimseldpz->Fill(sum.Pz());
}
for (j=0;j<dptselpsivtx.GetLength();++j) {
	hdptselpsivtxmass->Fill(dptselpsivtx[j].M());
	hdptselpsivtxpt->Fill(dptselpsivtx[j].Pt());
	hdptselpsivtxpz->Fill(dptselpsivtx[j].Pz());
	sum = dptselpsivtx[j].Daughter(0)->P4() + dptselpsivtx[j].Daughter(1)->P4();
	hdptselpsivtxdpt->Fill(sum.Pt());
	hdptselpsivtxdpz->Fill(sum.Pz());
}

// plot pt (Psi, momentum should be the same as D sum) selected Psi
for (j=0;j<pptselpsiraw.GetLength();++j) {
	hpptselpsirawmass->Fill(pptselpsiraw[j].M());
	hpptselpsirawpt->Fill(pptselpsiraw[j].Pt());
	hpptselpsirawpz->Fill(pptselpsiraw[j].Pz());
	sum = pptselpsiraw[j].Daughter(0)->P4() + pptselpsiraw[j].Daughter(1)->P4();
	hpptselpsirawdpt->Fill(sum.Pt());
	hpptselpsirawdpz->Fill(sum.Pz());
}
for (j=0;j<pptselpsimsel.GetLength();++j) {
	hpptselpsimselmass->Fill(pptselpsimsel[j].M());
	hpptselpsimselpt->Fill(pptselpsimsel[j].Pt());
	hpptselpsimselpz->Fill(pptselpsimsel[j].Pz());
	sum = pptselpsimsel[j].Daughter(0)->P4() + pptselpsimsel[j].Daughter(1)->P4();
	hpptselpsimseldpt->Fill(sum.Pt());
	hpptselpsimseldpz->Fill(sum.Pz());
}
for (j=0;j<pptselpsivtx.GetLength();++j) {
	hpptselpsivtxmass->Fill(pptselpsivtx[j].M());
	hpptselpsivtxpt->Fill(pptselpsivtx[j].Pt());
	hpptselpsivtxpz->Fill(pptselpsivtx[j].Pz());
	sum = pptselpsivtx[j].Daughter(0)->P4() + pptselpsivtx[j].Daughter(1)->P4();
	hpptselpsivtxdpt->Fill(sum.Pt());
	hpptselpsivtxdpz->Fill(sum.Pz());
}

// plot pt (both D and Psi, should be redundant) selected Psi
for (j=0;j<ptselpsiraw.GetLength();++j) {
	hptselpsirawmass->Fill(ptselpsiraw[j].M());
	hptselpsirawpt->Fill(ptselpsiraw[j].Pt());
	hptselpsirawpz->Fill(ptselpsiraw[j].Pz());
	sum = ptselpsiraw[j].Daughter(0)->P4() + ptselpsiraw[j].Daughter(1)->P4();
	hptselpsirawdpt->Fill(sum.Pt());
	hptselpsirawdpz->Fill(sum.Pz());
}
for (j=0;j<ptselpsimsel.GetLength();++j) {
	hptselpsimselmass->Fill(ptselpsimsel[j].M());
	hptselpsimselpt->Fill(ptselpsimsel[j].Pt());
	hptselpsimselpz->Fill(ptselpsimsel[j].Pz());
	sum = ptselpsimsel[j].Daughter(0)->P4() + ptselpsimsel[j].Daughter(1)->P4();
	hptselpsimseldpt->Fill(sum.Pt());
	hptselpsimseldpz->Fill(sum.Pz());
}
for (j=0;j<ptselpsivtx.GetLength();++j) {
	hptselpsivtxmass->Fill(ptselpsivtx[j].M());
	hptselpsivtxpt->Fill(ptselpsivtx[j].Pt());
	hptselpsivtxpz->Fill(ptselpsivtx[j].Pz());
	sum = ptselpsivtx[j].Daughter(0)->P4() + ptselpsivtx[j].Daughter(1)->P4();
	hptselpsivtxdpt->Fill(sum.Pt());
	hptselpsivtxdpz->Fill(sum.Pz());
}

// Make Psi Vertex Fit

for (j=0;j<ptselpsivtx.GetLength();++j) {
	PndKinVtxFitter fitter(ptselpsivtx[j]);
	fitter.Fit();
	TCandidate fitcand=*(fitter.FittedCand(ptselpsivtx[j]));
	TCandidate dpfit=*(fitter.FittedCand(*(ptselpsivtx[j].Daughter(0))));
	TCandidate dmfit=*(fitter.FittedCand(*(ptselpsivtx[j].Daughter(1))));
	//TCandidate* dpfit = fitcand->Daughter(0);
	//TCandidate* dmfit = fitcand->Daughter(1);
	TLorentzVector sum=dpfit.P4()+dmfit.P4();
	hvtxpsichi2->Fill(fitter.GlobalChi2());
	if (fitter.GlobalChi2()<50) {
		hvtxpsiacceptmass->Fill(fitcand.M());
		hvtxpsiacceptpt->Fill(fitcand.Pt());
		hvtxpsiacceptpz->Fill(fitcand.Pz());
		hvtxpsiacceptdpt->Fill(sum.Pt());
		hvtxpsiacceptdpz->Fill(sum.Pz());
		vtxpsi.Add(fitcand);
		//pt selection here for now
		if (fitcand.Pt() < 0.5) {
			hptselvtxpsimass->Fill(fitcand.M());
			hptselvtxpsipt->Fill(fitcand.Pt());
			hptselvtxpsipz->Fill(fitcand.Pz());
			hptselvtxpsidpt->Fill(sum.Pt());
			hptselvtxpsidpz->Fill(sum.Pz());
		}
	} else {
		hvtxpsirejectmass->Fill(fitcand.M());
		hvtxpsirejectpt->Fill(fitcand.Pt());
		hvtxpsirejectpz->Fill(fitcand.Pz());
		hvtxpsirejectdpt->Fill(sum.Pt());
		hvtxpsirejectdpz->Fill(sum.Pz());
	}
}

} // end main loop

//TLorentzVector ini(0,0,6.5788,7.58364);

//plot histograms    

TCanvas *eventstatisticscanvas = new TCanvas("eventstatisticscanvas","eventstatisticscanvas",600,600);
eventstatisticscanvas->Divide(2,3);
	eventstatisticscanvas->cd(1); hniceevents->Draw();
	eventstatisticscanvas->cd(3); hkpperevent->Draw();
	eventstatisticscanvas->cd(4); hkmperevent->Draw();
	eventstatisticscanvas->cd(5); hpipperevent->Draw();
	eventstatisticscanvas->cd(6); hpimperevent->Draw();

TCanvas *dcanvas=new TCanvas("dcanvas","dcanvas",600,600);
dcanvas->Divide(2,4);
	dcanvas->cd(1); hdprawmass->Draw();
	dcanvas->cd(2); hdmrawmass->Draw();
	dcanvas->cd(3); hdpmselmass->Draw();
	dcanvas->cd(4); hdmmselmass->Draw();
	dcanvas->cd(5); hdpvtxacceptmass->Draw(); hdpvtxrejectmass->Draw("same");
	dcanvas->cd(6); hdmvtxacceptmass->Draw(); hdmvtxrejectmass->Draw("same");
	dcanvas->cd(7); hdpvtxchi2->Draw(); //psimass3->Draw(); //ppmass->Draw();
	dcanvas->cd(8); hdmvtxchi2->Draw(); //ppmass->Draw();
	dcanvas->cd();

TCanvas *rawpsirawcanvas = new TCanvas("rawpsirawcanvas","rawpsirawcanvas",600,600);
rawpsirawcanvas->Divide(3,2);
	rawpsirawcanvas->cd(1); hrawpsirawmass->Draw(); hmselpsirawmass->Draw("same");
	rawpsirawcanvas->cd(2); hrawpsirawpt->Draw(); hmselpsirawpt->Draw("same");
	rawpsirawcanvas->cd(3); hrawpsirawpz->Draw(); hmselpsirawpz->Draw("same");
	rawpsirawcanvas->cd(5); hrawpsirawdpt->Draw(); hmselpsirawdpt->Draw("same");
	rawpsirawcanvas->cd(6); hrawpsirawdpz->Draw(); hmselpsirawdpz->Draw("same");

TCanvas *rawpsimselcanvas = new TCanvas("rawpsimselcanvas","rawpsimselcanvas",600,600);
rawpsimselcanvas->Divide(3,2);
	rawpsimselcanvas->cd(1); hrawpsimselmass->Draw(); hmselpsimselmass->Draw("same");
	rawpsimselcanvas->cd(2); hrawpsimselpt->Draw(); hmselpsimselpt->Draw("same");
	rawpsimselcanvas->cd(3); hrawpsimselpz->Draw(); hmselpsimselpz->Draw("same");
	rawpsimselcanvas->cd(5); hrawpsimseldpt->Draw(); hmselpsimseldpt->Draw("same");
	rawpsimselcanvas->cd(6); hrawpsimseldpz->Draw(); hmselpsimseldpz->Draw("same");

TCanvas *rawpsivtxcanvas = new TCanvas("rawpsivtxcanvas","rawpsivtxcanvas",600,600);
rawpsivtxcanvas->Divide(3,2);
	rawpsivtxcanvas->cd(1); hrawpsivtxmass->Draw(); hmselpsivtxmass->Draw("same");
	rawpsivtxcanvas->cd(2); hrawpsivtxpt->Draw(); hmselpsivtxpt->Draw("same");
	rawpsivtxcanvas->cd(3); hrawpsivtxpz->Draw(); hmselpsivtxpz->Draw("same");
	rawpsivtxcanvas->cd(4); hfinalpsimass->Draw();
	rawpsivtxcanvas->cd(5); hrawpsivtxdpt->Draw(); hmselpsivtxdpt->Draw("same");
	rawpsivtxcanvas->cd(6); hrawpsivtxdpz->Draw(); hmselpsivtxdpz->Draw("same");

TCanvas *dptselpsirawcanvas = new TCanvas("dptselpsirawcanvas","dptselpsirawcanvas",600,600);
dptselpsirawcanvas->Divide(3,2);
	dptselpsirawcanvas->cd(1); hdptselpsirawmass->Draw();
	dptselpsirawcanvas->cd(2); hdptselpsirawpt->Draw();
	dptselpsirawcanvas->cd(3); hdptselpsirawpz->Draw();
	dptselpsirawcanvas->cd(5); hdptselpsirawdpt->Draw();
	dptselpsirawcanvas->cd(6); hdptselpsirawdpz->Draw();

TCanvas *dptselpsimselcanvas = new TCanvas("dptselpsimselcanvas","dptselpsimselcanvas",600,600);
dptselpsimselcanvas->Divide(3,2);
	dptselpsimselcanvas->cd(1); hdptselpsimselmass->Draw();
	dptselpsimselcanvas->cd(2); hdptselpsimselpt->Draw();
	dptselpsimselcanvas->cd(3); hdptselpsimselpz->Draw();
	dptselpsimselcanvas->cd(5); hdptselpsimseldpt->Draw();
	dptselpsimselcanvas->cd(6); hdptselpsimseldpz->Draw();

TCanvas *dptselpsivtxcanvas = new TCanvas("dptselpsivtxcanvas","dptselpsivtxcanvas",600,600);
dptselpsivtxcanvas->Divide(3,2);
	dptselpsivtxcanvas->cd(1); hdptselpsivtxmass->Draw();
	dptselpsivtxcanvas->cd(2); hdptselpsivtxpt->Draw();
	dptselpsivtxcanvas->cd(3); hdptselpsivtxpz->Draw();
	dptselpsivtxcanvas->cd(5); hdptselpsivtxdpt->Draw();
	dptselpsivtxcanvas->cd(6); hdptselpsivtxdpz->Draw();

TCanvas *pptselpsirawcanvas = new TCanvas("pptselpsirawcanvas","pptselpsirawcanvas",600,600);
pptselpsirawcanvas->Divide(3,2);
	pptselpsirawcanvas->cd(1); hpptselpsirawmass->Draw();
	pptselpsirawcanvas->cd(2); hpptselpsirawpt->Draw();
	pptselpsirawcanvas->cd(3); hpptselpsirawpz->Draw();
	pptselpsirawcanvas->cd(5); hpptselpsirawdpt->Draw();
	pptselpsirawcanvas->cd(6); hpptselpsirawdpz->Draw();

TCanvas *pptselpsimselcanvas = new TCanvas("pptselpsimselcanvas","pptselpsimselcanvas",600,600);
pptselpsimselcanvas->Divide(3,2);
	pptselpsimselcanvas->cd(1); hpptselpsimselmass->Draw();
	pptselpsimselcanvas->cd(2); hpptselpsimselpt->Draw();
	pptselpsimselcanvas->cd(3); hpptselpsimselpz->Draw();
	pptselpsimselcanvas->cd(5); hpptselpsimseldpt->Draw();
	pptselpsimselcanvas->cd(6); hpptselpsimseldpz->Draw();

TCanvas *pptselpsivtxcanvas = new TCanvas("pptselpsivtxcanvas","pptselpsivtxcanvas",600,600);
pptselpsivtxcanvas->Divide(3,2);
	pptselpsivtxcanvas->cd(1); hpptselpsivtxmass->Draw();
	pptselpsivtxcanvas->cd(2); hpptselpsivtxpt->Draw();
	pptselpsivtxcanvas->cd(3); hpptselpsivtxpz->Draw();
	pptselpsivtxcanvas->cd(5); hpptselpsivtxdpt->Draw();
	pptselpsivtxcanvas->cd(6); hpptselpsivtxdpz->Draw();

TCanvas *ptselpsirawcanvas = new TCanvas("ptselpsirawcanvas","ptselpsirawcanvas",600,600);
ptselpsirawcanvas->Divide(3,2);
	ptselpsirawcanvas->cd(1); hptselpsirawmass->Draw();
	ptselpsirawcanvas->cd(2); hptselpsirawpt->Draw();
	ptselpsirawcanvas->cd(3); hptselpsirawpz->Draw();
	ptselpsirawcanvas->cd(5); hptselpsirawdpt->Draw();
	ptselpsirawcanvas->cd(6); hptselpsirawdpz->Draw();

TCanvas *ptselpsimselcanvas = new TCanvas("ptselpsimselcanvas","ptselpsimselcanvas",600,600);
ptselpsimselcanvas->Divide(3,2);
	ptselpsimselcanvas->cd(1); hptselpsimselmass->Draw();
	ptselpsimselcanvas->cd(2); hptselpsimselpt->Draw();
	ptselpsimselcanvas->cd(3); hptselpsimselpz->Draw();
	ptselpsimselcanvas->cd(5); hptselpsimseldpt->Draw();
	ptselpsimselcanvas->cd(6); hptselpsimseldpz->Draw();

TCanvas *ptselpsivtxcanvas = new TCanvas("ptselpsivtxcanvas","ptselpsivtxcanvas",600,600);
ptselpsivtxcanvas->Divide(3,2);
	ptselpsivtxcanvas->cd(1); hptselpsivtxmass->Draw();
	ptselpsivtxcanvas->cd(2); hptselpsivtxpt->Draw();
	ptselpsivtxcanvas->cd(3); hptselpsivtxpz->Draw();
	ptselpsivtxcanvas->cd(5); hptselpsivtxdpt->Draw();
	ptselpsivtxcanvas->cd(6); hptselpsivtxdpz->Draw();

TCanvas* vtxpsicanvas = new TCanvas("vtxpsicanvas","vtxpsicanvas",600,600);
vtxpsicanvas->Divide(3,2);
	vtxpsicanvas->cd(1); hvtxpsichi2->Draw();
	vtxpsicanvas->cd(2); hvtxpsiacceptmass->Draw(); hvtxpsirejectmass->Draw("same");
	vtxpsicanvas->cd(3); hvtxpsiacceptpt->Draw(); hvtxpsirejectpt->Draw("same");
	vtxpsicanvas->cd(4); hvtxpsiacceptpz->Draw(); hvtxpsirejectpz->Draw("same");
	vtxpsicanvas->cd(5); hvtxpsiacceptdpt->Draw(); hvtxpsirejectdpt->Draw("same");
	vtxpsicanvas->cd(6); hvtxpsiacceptdpz->Draw(); hvtxpsirejectdpz->Draw("same");

TCanvas* ptselvtxpsicanvas = new TCanvas("ptselvtxpsicanvas","ptselvtxpsicanvas",600,600);
ptselvtxpsicanvas->Divide(2,3);
	ptselvtxpsicanvas->cd(1); hptselvtxpsimass->Draw();
	ptselvtxpsicanvas->cd(3); hptselvtxpsipt->Draw();
	ptselvtxpsicanvas->cd(4); hptselvtxpsipz->Draw();
	ptselvtxpsicanvas->cd(5); hptselvtxpsidpt->Draw();
	ptselvtxpsicanvas->cd(6); hptselvtxpsidpz->Draw();;

TCanvas* vertexresolutioncanvas = new TCanvas("vertexresolutioncanvas","vertexresolutioncanvas",600,600);
vertexresolutioncanvas->Divide(2,2);
	vertexresolutioncanvas->cd(1); hdpvertex->Draw();
	vertexresolutioncanvas->cd(2); hdpvertexfit->Draw();
	vertexresolutioncanvas->cd(3); hdmvertex->Draw();
	vertexresolutioncanvas->cd(4); hdmvertexfit->Draw();

TCanvas* vertexresolutioncanvas2 = new TCanvas("vertexresolutioncanvas2","vertexresolutioncanvas2",600,600);
vertexresolutioncanvas2->Divide(3,3);
	vertexresolutioncanvas2->cd(1); hdmvertexx->Draw();
	vertexresolutioncanvas2->cd(2); hdmvertexy->Draw();
	vertexresolutioncanvas2->cd(3); hdmvertexz->Draw();
	vertexresolutioncanvas2->cd(4); hdmvertexxfit->Draw();
	vertexresolutioncanvas2->cd(5); hdmvertexyfit->Draw();
	vertexresolutioncanvas2->cd(6); hdmvertexzfit->Draw();
	vertexresolutioncanvas2->cd(7); hdmvertexxpocareso->Draw();
	vertexresolutioncanvas2->cd(8); hdmvertexypocareso->Draw();
	vertexresolutioncanvas2->cd(9); hdmvertexzpocareso->Draw();

TCanvas* vertexresolutioncanvas3 = new TCanvas("vertexresolutioncanvas3","vertexresolutioncanvas3",600,600);
vertexresolutioncanvas3->Divide(3,3);
	vertexresolutioncanvas3->cd(1); hdpvertexx->Draw();
	vertexresolutioncanvas3->cd(2); hdpvertexy->Draw();
	vertexresolutioncanvas3->cd(3); hdpvertexz->Draw();
	vertexresolutioncanvas3->cd(4); hdpvertexxfit->Draw();
	vertexresolutioncanvas3->cd(5); hdpvertexyfit->Draw();
	vertexresolutioncanvas3->cd(6); hdpvertexzfit->Draw();
	vertexresolutioncanvas3->cd(7); hdpvertexxpocareso->Draw();
	vertexresolutioncanvas3->cd(8); hdpvertexypocareso->Draw();
	vertexresolutioncanvas3->cd(9); hdpvertexzpocareso->Draw();

TCanvas* vertexpositionpcanvas = new TCanvas("vertexpositionpcanvas","vertexpositionpcanvas",600,600);
vertexpositionpcanvas->Divide(3,3);
	vertexpositionpcanvas->cd(1); hdpvertexxmc->Draw();
	vertexpositionpcanvas->cd(2); hdpvertexymc->Draw();
	vertexpositionpcanvas->cd(3); hdpvertexzmc->Draw();
	vertexpositionpcanvas->cd(4); hdpvertexxreco->Draw();
	vertexpositionpcanvas->cd(5); hdpvertexyreco->Draw();
	vertexpositionpcanvas->cd(6); hdpvertexzreco->Draw();
	vertexpositionpcanvas->cd(7); hdpvertexxpoca->Draw();
	vertexpositionpcanvas->cd(8); hdpvertexypoca->Draw();
	vertexpositionpcanvas->cd(9); hdpvertexzpoca->Draw();

TCanvas* vertexpositionmcanvas = new TCanvas("vertexpositionmcanvas","vertexpositionmcanvas",600,600);
vertexpositionmcanvas->Divide(3,3);
	vertexpositionmcanvas->cd(1); hdmvertexxmc->Draw();
	vertexpositionmcanvas->cd(2); hdmvertexymc->Draw();
	vertexpositionmcanvas->cd(3); hdmvertexzmc->Draw();
	vertexpositionmcanvas->cd(4); hdmvertexxreco->Draw();
	vertexpositionmcanvas->cd(5); hdmvertexyreco->Draw();
	vertexpositionmcanvas->cd(6); hdmvertexzreco->Draw();
	vertexpositionmcanvas->cd(7); hdmvertexxpoca->Draw();
	vertexpositionmcanvas->cd(8); hdmvertexypoca->Draw();
	vertexpositionmcanvas->cd(9); hdmvertexzpoca->Draw();


TCanvas* decaylengthcanvas = new TCanvas("decaylengthcanvas","decaylengthcanvas",600,600);
decaylengthcanvas->Divide(2,1);
	decaylengthcanvas->cd(1); hdpdecaylength->Draw();
	decaylengthcanvas->cd(2); hdmdecaylength->Draw();

//save and update the interesting histos
outputstorage.cd();
SaveAndUpdateHisto(hniceevents, outputstorage);
SaveAndUpdateHisto(hkpperevent, outputstorage);
SaveAndUpdateHisto(hkmperevent, outputstorage);
SaveAndUpdateHisto(hpipperevent, outputstorage);
SaveAndUpdateHisto(hpimperevent, outputstorage);

SaveAndUpdateHisto(hdprawmass, outputstorage);
SaveAndUpdateHisto(hdmrawmass, outputstorage);
SaveAndUpdateHisto(hdpmselmass, outputstorage);
SaveAndUpdateHisto(hdmmselmass, outputstorage);
SaveAndUpdateHisto(hdpvtxacceptmass, outputstorage);
SaveAndUpdateHisto(hdpvtxrejectmass, outputstorage);
SaveAndUpdateHisto(hdmvtxacceptmass, outputstorage);
SaveAndUpdateHisto(hdmvtxrejectmass, outputstorage);
SaveAndUpdateHisto(hdpvtxchi2, outputstorage);
SaveAndUpdateHisto(hdmvtxchi2, outputstorage);

SaveAndUpdateHisto(hrawpsirawmass, outputstorage);
SaveAndUpdateHisto(hmselpsirawmass, outputstorage);
SaveAndUpdateHisto(hrawpsirawpt, outputstorage);
SaveAndUpdateHisto(hmselpsirawpt, outputstorage);
SaveAndUpdateHisto(hrawpsirawpz, outputstorage);
SaveAndUpdateHisto(hmselpsirawpz, outputstorage);
SaveAndUpdateHisto(hrawpsirawdpt, outputstorage);
SaveAndUpdateHisto(hmselpsirawdpt, outputstorage);
SaveAndUpdateHisto(hrawpsirawdpz, outputstorage);
SaveAndUpdateHisto(hmselpsirawdpz, outputstorage);
SaveAndUpdateHisto(hrawpsimselmass, outputstorage);
SaveAndUpdateHisto(hmselpsimselmass, outputstorage);
SaveAndUpdateHisto(hrawpsimselpt, outputstorage);
SaveAndUpdateHisto(hmselpsimselpt, outputstorage);
SaveAndUpdateHisto(hrawpsimselpz, outputstorage);
SaveAndUpdateHisto(hmselpsimselpz, outputstorage);
SaveAndUpdateHisto(hrawpsimseldpt, outputstorage);
SaveAndUpdateHisto(hmselpsimseldpt, outputstorage);
SaveAndUpdateHisto(hrawpsimseldpz, outputstorage);
SaveAndUpdateHisto(hmselpsimseldpz, outputstorage);
SaveAndUpdateHisto(hrawpsivtxmass, outputstorage);
SaveAndUpdateHisto(hmselpsivtxmass, outputstorage);
SaveAndUpdateHisto(hrawpsivtxpt, outputstorage);
SaveAndUpdateHisto(hmselpsivtxpt, outputstorage);
SaveAndUpdateHisto(hrawpsivtxpz, outputstorage);
SaveAndUpdateHisto(hmselpsivtxpz, outputstorage);
SaveAndUpdateHisto(hfinalpsimass, outputstorage);
SaveAndUpdateHisto(hrawpsivtxdpt, outputstorage);
SaveAndUpdateHisto(hmselpsivtxdpt, outputstorage);
SaveAndUpdateHisto(hrawpsivtxdpz, outputstorage);
SaveAndUpdateHisto(hmselpsivtxdpz, outputstorage);

SaveAndUpdateHisto(hdptselpsirawmass, outputstorage);
SaveAndUpdateHisto(hdptselpsirawpt, outputstorage);
SaveAndUpdateHisto(hdptselpsirawpz, outputstorage);
SaveAndUpdateHisto(hdptselpsirawdpt, outputstorage);
SaveAndUpdateHisto(hdptselpsirawdpz, outputstorage);
SaveAndUpdateHisto(hdptselpsimselmass, outputstorage);
SaveAndUpdateHisto(hdptselpsimselpt, outputstorage);
SaveAndUpdateHisto(hdptselpsimselpz, outputstorage);
SaveAndUpdateHisto(hdptselpsimseldpt, outputstorage);
SaveAndUpdateHisto(hdptselpsimseldpz, outputstorage);
SaveAndUpdateHisto(hdptselpsivtxmass, outputstorage);
SaveAndUpdateHisto(hdptselpsivtxpt, outputstorage);
SaveAndUpdateHisto(hdptselpsivtxpz, outputstorage);
SaveAndUpdateHisto(hdptselpsivtxdpt, outputstorage);
SaveAndUpdateHisto(hdptselpsivtxdpz, outputstorage);

SaveAndUpdateHisto(hpptselpsirawmass, outputstorage);
SaveAndUpdateHisto(hpptselpsirawpt, outputstorage);
SaveAndUpdateHisto(hpptselpsirawpz, outputstorage);
SaveAndUpdateHisto(hpptselpsirawdpt, outputstorage);
SaveAndUpdateHisto(hpptselpsirawdpz, outputstorage);
SaveAndUpdateHisto(hpptselpsimselmass, outputstorage);
SaveAndUpdateHisto(hpptselpsimselpt, outputstorage);
SaveAndUpdateHisto(hpptselpsimselpz, outputstorage);
SaveAndUpdateHisto(hpptselpsimseldpt, outputstorage);
SaveAndUpdateHisto(hpptselpsimseldpz, outputstorage);
SaveAndUpdateHisto(hpptselpsivtxmass, outputstorage);
SaveAndUpdateHisto(hpptselpsivtxpt, outputstorage);
SaveAndUpdateHisto(hpptselpsivtxpz, outputstorage);
SaveAndUpdateHisto(hpptselpsivtxdpt, outputstorage);
SaveAndUpdateHisto(hpptselpsivtxdpz, outputstorage);

SaveAndUpdateHisto(hptselpsirawmass, outputstorage);
SaveAndUpdateHisto(hptselpsirawpt, outputstorage);
SaveAndUpdateHisto(hptselpsirawpz, outputstorage);
SaveAndUpdateHisto(hptselpsirawdpt, outputstorage);
SaveAndUpdateHisto(hptselpsirawdpz, outputstorage);
SaveAndUpdateHisto(hptselpsimselmass, outputstorage);
SaveAndUpdateHisto(hptselpsimselpt, outputstorage);
SaveAndUpdateHisto(hptselpsimselpz, outputstorage);
SaveAndUpdateHisto(hptselpsimseldpt, outputstorage);
SaveAndUpdateHisto(hptselpsimseldpz, outputstorage);
SaveAndUpdateHisto(hptselpsivtxmass, outputstorage);
SaveAndUpdateHisto(hptselpsivtxpt, outputstorage);
SaveAndUpdateHisto(hptselpsivtxpz, outputstorage);
SaveAndUpdateHisto(hptselpsivtxdpt, outputstorage);
SaveAndUpdateHisto(hptselpsivtxdpz, outputstorage);

SaveAndUpdateHisto(hvtxpsichi2, outputstorage);
SaveAndUpdateHisto(hvtxpsiacceptmass, outputstorage);
SaveAndUpdateHisto(hvtxpsirejectmass, outputstorage);
SaveAndUpdateHisto(hvtxpsiacceptpt, outputstorage);
SaveAndUpdateHisto(hvtxpsirejectpt, outputstorage);
SaveAndUpdateHisto(hvtxpsiacceptpz, outputstorage);
SaveAndUpdateHisto(hvtxpsirejectpz, outputstorage);
SaveAndUpdateHisto(hvtxpsiacceptdpt, outputstorage);
SaveAndUpdateHisto(hvtxpsirejectdpt, outputstorage);
SaveAndUpdateHisto(hvtxpsiacceptdpz, outputstorage);
SaveAndUpdateHisto(hvtxpsirejectdpz, outputstorage);

SaveAndUpdateHisto(hptselvtxpsimass, outputstorage);
SaveAndUpdateHisto(hptselvtxpsipt, outputstorage);
SaveAndUpdateHisto(hptselvtxpsipz, outputstorage);
SaveAndUpdateHisto(hptselvtxpsidpt, outputstorage);
SaveAndUpdateHisto(hptselvtxpsidpz, outputstorage);
SaveAndUpdateHisto(hdpvertex, outputstorage);
SaveAndUpdateHisto(hdpvertexfit, outputstorage);
SaveAndUpdateHisto(hdmvertex, outputstorage);
SaveAndUpdateHisto(hdmvertexfit, outputstorage);

SaveAndUpdateHisto(hdmvertexx, outputstorage);
SaveAndUpdateHisto(hdmvertexy, outputstorage);
SaveAndUpdateHisto(hdmvertexz, outputstorage);
SaveAndUpdateHisto(hdmvertexxfit, outputstorage);
SaveAndUpdateHisto(hdmvertexyfit, outputstorage);
SaveAndUpdateHisto(hdmvertexzfit, outputstorage);
SaveAndUpdateHisto(hdpvertexx, outputstorage);
SaveAndUpdateHisto(hdpvertexy, outputstorage);
SaveAndUpdateHisto(hdpvertexz, outputstorage);
SaveAndUpdateHisto(hdpvertexxfit, outputstorage);
SaveAndUpdateHisto(hdpvertexyfit, outputstorage);
SaveAndUpdateHisto(hdpvertexzfit, outputstorage);

SaveAndUpdateHisto(hdpvertexxmc, outputstorage);
SaveAndUpdateHisto(hdpvertexymc, outputstorage);
SaveAndUpdateHisto(hdpvertexzmc, outputstorage);
SaveAndUpdateHisto(hdpvertexxreco, outputstorage);
SaveAndUpdateHisto(hdpvertexyreco, outputstorage);
SaveAndUpdateHisto(hdpvertexzreco, outputstorage);
SaveAndUpdateHisto(hdmvertexxmc, outputstorage);
SaveAndUpdateHisto(hdmvertexymc, outputstorage);
SaveAndUpdateHisto(hdmvertexzmc, outputstorage);
SaveAndUpdateHisto(hdmvertexxreco, outputstorage);
SaveAndUpdateHisto(hdmvertexyreco, outputstorage);
SaveAndUpdateHisto(hdmvertexzreco, outputstorage);
SaveAndUpdateHisto(hdpdecaylength, outputstorage);
SaveAndUpdateHisto(hdmdecaylength, outputstorage);

SaveAndUpdateHisto(hdpvertexxpoca, outputstorage);
SaveAndUpdateHisto(hdpvertexypoca, outputstorage);
SaveAndUpdateHisto(hdpvertexzpoca, outputstorage);
SaveAndUpdateHisto(hdmvertexxpoca, outputstorage);
SaveAndUpdateHisto(hdmvertexypoca, outputstorage);
SaveAndUpdateHisto(hdmvertexzpoca, outputstorage);
SaveAndUpdateHisto(hdpvertexxpocareso, outputstorage);
SaveAndUpdateHisto(hdpvertexypocareso, outputstorage);
SaveAndUpdateHisto(hdpvertexzpocareso, outputstorage);
SaveAndUpdateHisto(hdmvertexxpocareso, outputstorage);
SaveAndUpdateHisto(hdmvertexypocareso, outputstorage);
SaveAndUpdateHisto(hdmvertexzpocareso, outputstorage);

timer.Stop();
Double_t rtime = timer.RealTime();
Double_t ctime = timer.CpuTime();
    
printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
    
} // end macro

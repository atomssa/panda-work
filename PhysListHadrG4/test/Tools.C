// This is a collection of useful tool functions for root macros
// mostly it is to make it more eye-candy and easier to work
// R. Kliemt, 2008
// Thanks go to Thomas Goepfert & Markus Warsinsky from IKTP Dresden for BetterStatBox & BuildLegend_THStack
// Found set_nicer_2d_plot_style on the web
// Content:
//   DrawText(posX,posY,text,size,color)
//     R. Kliemt, 2008
//     Drawing some TLatex text in a root canvas, overlaying the content
//   DrawNice2DHisto(TH2*,opt="",range)
//     R. Kliemt, 2008
//     Set manually the maximum range & better labeling
//   set_nicer_2d_plot_style()
//     From the web somewhere
//     setting a smoother rainbow color profile to gStyle
//   BetterStatBox(TPad*)
//     Thanks go to Thomas Goepfert & Markus Warsinsky from IKTP Dresden
//     Disentangle & colorize statboxes from different histograms in one plot
//     1) Draw histograms with "sames" option, 2) call Updade() the canvas
//     3) use TPad* p=(TPad*)gPad; BetterStatBox(p);
//   TLegend* BuildLegend_THStack(THStack*,x1,y1,x2,y2 )
//     Thanks go to Thomas Goepfert & Markus Warsinsky from IKTP Dresden
//     Making a legend for staccked histograms. It still has to be drawn.
//   LoadPandaStyle()
//     Sets the formerly official Panda styling to histograms, fonts etc.
//     Note that this is a modified version to remove clipping bugs and to set
//     the nicer rainbow plots.
//   TH1D TransformHisto(TH2*,min,max)
//     T.Stockmanns, 2009
//     transforming from a 2D plot in a 1D plot, using the root sorting of the bins
//   plothistosfromfile(TString filename = "histos.root", TString ext=".ps", Int_t divx=2, Int_t divy=2, Int_t pix = 300)
//     R.Kliemt, 2012
//     Plots histograms from a root file to, e.g., a pdf. Several pages are created
//   TString InitDefaultRun(TString filetag)
//     R.Kliemt, 2011
//     Initialize a default named set of files and a FairRunAna object
//   DrawHistSecondScale(TH1* histo, int color)
//     R.Kliemt, 2013
//     Draw a histogram overlaying as with the "same" option, but scaling it to
//     the right size and adding the proper axis on the right
//


#include <TLatex.h>
#include <TColor.h>
#include <TString.h>
#include <TPaveStats.h>
#include <TH1.h>
#include <TLegend.h>
#include <THStack.h>
#include <TList.h>
#include <TStyle.h>
#include <TROOT.h>

void Tools()
{
  cout<<"pandaroot/macro/run/Tools.C loaded. Enjoy it."<<endl;
  return;
}

void DrawText(Double_t posX = 0., Double_t posY = 0., const char* text = "",
              Double_t size=0.08, Int_t col=1
              //                      Int_t align=0, Double_t angle=0.,
              //                      Int_t font, Bool_t bNDC=kTRUE
              )
{
  // Drawing some TLatex text in a root canvas, overlaying the content

  TLatex* pText=new TLatex(posX,posY,text);
  pText->SetNDC(kTRUE);
  //   pText->SetNDC(bNDC);
  pText->SetTextColor(col); pText->SetTextSize(size);
  //   pText->SetTextAngle(angle);
  //   pText->SetTextFont(font); pText->SetTextAlign(align);
  pText->Draw();
}

void DrawNice2DHisto(TH2* h,const char* opt="",double range = 10.)
{
  // Draw a 2D histo with the rainbow colors and the palette besides
  TString options = "colz"; options += opt;
  if(gPad->GetLogz() == 1 && h->GetMaximum()<range) h->SetAxisRange(0.,range,"Z");
  else {
    // we squeeze here to see better the mean value
    //range = 10*h->Integral()/(h->GetNbinsX()/h->GetNbinsY());
    cout << "range "<<range<<"   max "<<h->GetMaximum()<<"   min "<<h->GetMinimum()<<endl;
    if( range < h->GetMaximum()) range = h->GetMaximum();
    if(-range > h->GetMinimum()) range = -1*h->GetMinimum();
    range *= 1.1; // zoom a bit out to see the max better
    if(h->GetMinimum() < 0.) h->SetAxisRange(-range,range,"Z");
    else h->SetAxisRange(0,range,"Z");
  }
  h->SetStats(kFALSE);
  h->SetTitleOffset(0.8,"T");
  gPad->SetRightMargin(0.15);
  h->DrawCopy(options.Data());
}

void set_nicer_2d_plot_style()
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 99;//255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //     gStyle->CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void BetterStatBox( TPad* pad ){
  /// Vorgehensweise:
  /// * Histogramme in ein TPad zeichnen (1. Histogramm: Draw(), ntes Histogramm:
  ///   Draw("sames"))
  /// * die Update() Funktion vom Canvas aufrufen (die Statistikboxen werden erst
  ///   erzeugt, wenn das Histogramm wirklich gezeichnet wird)
  /// * der Funktion das aktuelle TPad uebergeben
  ///   (z.B.: TPad* p=(TPad*)gPad; BetterStatBox(p);)

  int entries  = pad->GetListOfPrimitives()->GetEntries();
  TString name = "";
  TPaveStats *st;
  int obj = 1;

  for(int i=0; i<entries; i++){
    name = pad->GetListOfPrimitives()->At(i)->IsA()->GetName();
    if( name.Contains("TH1") ){
      TH1* hist = (TH1*)pad->GetListOfPrimitives()->At(i);
      st = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
      if(st){
        st->SetOptStat(1);
        st->SetX1NDC(0.80);
        st->SetY1NDC(1-0.16*obj) ;
        st->SetY2NDC( st->GetY1NDC()+0.15 ) ;
        st->SetTextColor( hist->GetLineColor() );
        st->SetLineColor( hist->GetLineColor() );
        st->Draw("same");
        obj++;
      }
    }
  }
}

// TLegend discontinued in ROOT ???
//TLegend* BuildLegend_THStack( THStack* stack, float x1, float y1, float x2, float y2 ){
//
//  TLegend* legend = new TLegend(x1,y1,x2,y2);
//  TList*   list = stack->GetHists();
//  TIter    next( list );
//  TH1*     hist;
//
//  while ( hist = (TH1*)next() ) legend->AddEntry(hist,"","F");
//
//  return legend;
//}

void ImproveDefaultStyle()
{
  gStyle->SetOptFit ( 1011 );
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.14);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetTextFont(22);//changed to bold (R.K.)
  gStyle->SetTextSize(0.08);
  gStyle->SetLabelFont(22,"x");//changed to bold (R.K.)
  gStyle->SetLabelFont(22,"y");//changed to bold (R.K.)
  gStyle->SetLabelFont(22,"z");//changed to bold (R.K.)
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");
  // use bold lines and markers
  gStyle->SetMarkerStyle(8);
  gStyle->SetHistLineWidth(1.85);//1.85);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //R.K. avoid clumsy axis lables
  gStyle->SetNdivisions(509); // default root value is 510

}

void LoadPandaStyle(void)
{

  //--------------------------------------------------------------------------
  // File and Version Information:
  //      $Id: PBase.C,v 1.8 2006/09/22 12:04:13 kliemt Exp $
  //
  // Description:
  //      Initialization code executed at the start of a ROOT session.
  //      Set up the Panda style for approved plots.
  //
  // Environment:
  //      Software developed for the PANDA Detector at GSI, Darmstadt
  //
  // Author List:
  //      Sergey Ganzhur                Original Author
  //      Ralf Kliemt (2008)            Small adjustments for PandaRoot use
  //
  // Copyright Information:
  //     Copyright (C) 2001-2002        Ruhr Universitaet Bochum
  //
  //------------------------------------------------------------------------

  if(gStyle->GetName() == "PANDA") return;
  // use the 'plain' style for plots (white backgrounds, etc)
  //cout << "...using style 'Plain'" << endl;
  //gROOT->SetStyle("Plain");

  // Create the 'PANDA' style for approved plots. Note that this style may need
  // some fine tuning in your macro depending on what you are plotting, e.g.
  //
  //  gStyle->SetMarkerSize(0.75);  // use smaller markers in a histogram with many bins
  //  gStyle->SetTitleOffset(0.65,"y");  // bring y axis label closer to narrow values


  // Ralf Kliemt:
  // I changed a bit for myself here

  TStyle *pandaStyle= new TStyle("PANDA","PANDA approved plots style, modified by ralfk");

  // use plain black on white colors
  pandaStyle->SetFrameBorderMode(0);
  pandaStyle->SetCanvasBorderMode(0);
  pandaStyle->SetPadBorderMode(0);
  pandaStyle->SetPadColor(0);
  pandaStyle->SetCanvasColor(0);
  pandaStyle->SetTitleColor(1);
  pandaStyle->SetStatColor(0);
  //   pandaStyle->SetFillColor(0);// conflict with 2D plots
  //   pandaStyle->SetPalette(1);
  //R.K: Now remove the title box
  //pandaStyle->SetTitleAlign(33);
  pandaStyle->SetTitleBorderSize(0);
  pandaStyle->SetTitleColor(1);
  pandaStyle->SetTitleFillColor(0);
  pandaStyle->SetTitleFontSize(0.07);


  // set the paper & margin sizes
  pandaStyle->SetPaperSize(20,26);
  pandaStyle->SetPadTopMargin(0.14);
  pandaStyle->SetPadRightMargin(0.05);
  pandaStyle->SetPadBottomMargin(0.14);
  pandaStyle->SetPadLeftMargin(0.14);
//  pandaStyle->SetPaperSize(20,26);
//  pandaStyle->SetPadTopMargin(0.05);
//  pandaStyle->SetPadRightMargin(0.05);
//  pandaStyle->SetPadBottomMargin(0.16);
//  pandaStyle->SetPadLeftMargin(0.12);

  // use large Times-Roman fonts
  //   pandaStyle->SetTextFont(132);
  //   pandaStyle->SetLabelFont(132,"x");
  //   pandaStyle->SetLabelFont(132,"y");
  //   pandaStyle->SetLabelFont(132,"z");
  pandaStyle->SetTextFont(22);//changed to bold (R.K.)
  pandaStyle->SetTextSize(0.08);
  pandaStyle->SetLabelFont(22,"x");//changed to bold (R.K.)
  pandaStyle->SetLabelFont(22,"y");//changed to bold (R.K.)
  pandaStyle->SetLabelFont(22,"z");//changed to bold (R.K.)
  pandaStyle->SetLabelSize(0.05,"x");
  pandaStyle->SetTitleSize(0.06,"x");
  pandaStyle->SetLabelSize(0.05,"y");
  pandaStyle->SetTitleSize(0.06,"y");
  pandaStyle->SetLabelSize(0.05,"z");
  pandaStyle->SetTitleSize(0.06,"z");

  // use bold lines and markers
  pandaStyle->SetMarkerStyle(8);
  pandaStyle->SetHistLineWidth(1.85);//1.85);
  pandaStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // do not display any of the standard histogram decorations
  pandaStyle->SetOptTitle(1);
  pandaStyle->SetOptStat(1);
  pandaStyle->SetOptFit(1);

  // put tick marks on top and RHS of plots
  pandaStyle->SetPadTickX(1);
  pandaStyle->SetPadTickY(1);

  //R.K. avoid clumsy axis lables
  pandaStyle->SetNdivisions(509); // default root value is 510

  //cout <<"    For approved plots use: gROOT->SetStyle(\"PANDA\");"<< endl;

  pandaStyle->cd();
  gROOT->ForceStyle();
  set_nicer_2d_plot_style();//[R.K.] use nicer 2D plots
  return;
}//void PBase::LoadPandaStyle(void)




TH1D TransformHisto(TH2* h2, double min, double max)
{
  /*
   *  Created on: Feb 25, 2009
   *      Author: stockman
   */

	TH1D result("h1","h1", 1000, min, max);
	int nbins = h2->GetNbinsX() * h2->GetNbinsY();
	for (int i = 0; i < nbins; i++){
		//std::cout << h2->GetBinContent(i) << std::endl;
		result->Fill(h2->GetBinContent(i));
		if (i == 10)
			cout << h2->GetBinContent(i) << endl;
	}
	return result;
}


plothistosfromfile(TString filename = "histos.root", TString ext=".ps", Int_t divx=2, Int_t divy=2, Int_t pix = 300)
{ // Plot all histograms into a ps file
  // works with TH1, TH2, & TProfile
  //LoadPandaStyle();
  set_nicer_2d_plot_style();
  TFile* file = new TFile(filename.Data());
  if (!file) {cout<<"File \""<<filename.Data()<<"\" is not there..."<<endl;return;}
  TCanvas* can = new TCanvas();
  Int_t pixx = ceil(1.4*pix*divx);
  Int_t pixy = pix*divy;
  can->SetCanvasSize(pixx,pixy);
  can->Divide(divx, divy);
  TString picname = filename;
  ext="."+ext;
  ext.ReplaceAll("..",".");
  picname.ReplaceAll(".root",ext); // ps, png, pdf ...
  TString pic = picname + "["; // open empty ps
  cout << "opening: " << pic.Data()<<endl;
  can->Print(pic);
  pic=picname;

  TList* list = file->GetListOfKeys();
  if (!list) {cout<<"List not there..."<<endl;return;}
  int padcount = 1;
  TString keyclass="";
  for(int i=0;i<list->GetEntries();i++)
  {
    if(padcount > divx*divy)
    {
      can->Print(pic.Data());
      can->Clear();
      can->SetCanvasSize(pixx,pixy);
      can->Divide(divx, divy);
      padcount=1;
    }
    can->cd(padcount);
    TKey* key = (TKey*)list->At(i);
    keyclass = key->GetClassName();
    //cout<<keyclass.Data()<<endl;
    if(keyclass.Contains("TH1"))
    {
      //cout<<"try plotting a TH1"<<endl;
      TH1* his = (TH1*)key->ReadObj();
      his->GetXaxis()->SetNoExponent(); // put exponents to numbers directly
      his->GetYaxis()->SetNoExponent(); // put exponents to numbers directly
      TString options(his->GetOption());
      if(options.Contains("log")){
        gPad->SetLogy();
        options.ReplaceAll("log","");
        his->SetOption(options.Data());
      }
      his->Draw();
    }else if(keyclass.Contains("TH2"))
    {
      //cout<<"try plotting a TH2"<<endl;
      TH2* his2 = (TH2*)key->ReadObj();
      his2->GetXaxis()->SetNoExponent(); // put exponents to numbers directly
      his2->GetYaxis()->SetNoExponent(); // put exponents to numbers directly
      TString options(his->GetOption());
      if(options.Contains("log")){
        gPad->SetLogz();
        options.ReplaceAll("log","");
        his->SetOption(options.Data());
      }
      if(options.Contains("nice")){
        options.ReplaceAll("nice","");
        his->SetOption(options.Data());
        DrawNice2DHisto(his2);
      } else {
        his2->Draw("colz");
      }
    }else if(keyclass.Contains("TProfile"))
    {
      //cout<<"try plotting a TH2"<<endl;
      TProfile* hpro = (TProfile*)key->ReadObj();
      hpro->GetXaxis()->SetNoExponent(); // put exponents to numbers directly
      hpro->GetYaxis()->SetNoExponent(); // put exponents to numbers directly
      TString options(hpro->GetOption());
      if(options.Contains("log")){
        gPad->SetLogy();
        options.ReplaceAll("log","");
        hpro->SetOption(options.Data());
      }
      hpro->Draw();
    } else continue;

    padcount++;
  }

  can->Print(pic.Data());
  pic = picname + "]"; // close ps
  can->Print(pic.Data());
  cout << "closed: " << pic.Data()<<endl;
  TString convertcmd = "test -r ps2pdf && ps2pdf ";
  convertcmd += pic.Data();
  gSystem->Exec(convertcmd.Data());
  delete can;
  return;
}

void LoadManySimFiles(TString treename="cbmsim")
{ // to use that method you should have opened some files
  // containing the same tree structure, like splitted files of
  // mass production simulations. (like "root -f data/sim01*.root"
  // tree examination is available via TTree::Draw(...)

  TIter next(gROOT->GetListOfFiles());
  TFile *fi=0;
  TChain *R=new TChain(treename.Data());
  while (fi=(TFile*)next()) R->Add(fi->GetName());
  cout<<(Int_t)R->GetEntries()<<endl;
}

TString InitDefaultRun(TString filetag)
{
  cout << "-I- Using InitDefaultRun() from macro/run/Tools.C with the sim file " << filetag.Data() << endl;
  FairRunAna* fRun = new FairRunAna();
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();

  filetag.ReplaceAll("_sim.root",".root");
  PndFileNameCreator namecreator(filetag.Data());
  TString simFile = namecreator.GetSimFileName();
  TString parFile = namecreator.GetParFileName();
  TString recoFile = namecreator.GetRecoFileName();
  TString tracksFile = namecreator.GetCustomFileName("tracks");
  TString pidFile = namecreator.GetCustomFileName("pid");
  TString histoFile = namecreator.GetCustomFileName("histos");
  TString evrdummy = namecreator.GetCustomFileName("evrdummy");
  std::cout<<"simFile="<<simFile.Data()<<std::endl;
  std::cout<<"parFile="<<parFile.Data()<<std::endl;
  std::cout<<"recoFile="<<recoFile.Data()<<std::endl;
  std::cout<<"trkFile="<<tracksFile.Data()<<std::endl;
  std::cout<<"pidFile="<<pidFile.Data()<<std::endl;
  fRun->SetInputFile(simFile);
  fRun->AddFriend(recoFile);//,"rec");
  fRun->AddFriend(tracksFile);//,"trk");
  fRun->AddFriend(pidFile);//,"pid");
  FairParRootFileIo* parIO = new FairParRootFileIo();
  parIO->open(parFile.Data());
  rtdb->setFirstInput(parIO);
  rtdb->setOutput(parIO);

  TString outFile = evrdummy;
  TString sysFile = gSystem->Getenv("VMCWORKDIR");

  fRun->SetOutputFile(outFile.Data());
  FairGeane* geane = new FairGeane();
  fRun->AddTask(geane);
  fRun->Init();
  return histoFile;
}


DrawHistSecondScale(TH1* hist, int color=kRed)
{
    // scale hint1 to the pad coordinates
   Float_t rightmax = 1.1*hist->GetMaximum();
   Float_t scale = gPad->GetUymax()/rightmax;
   hist->SetLineColor(color);
   hist->Scale(scale);
   hist->Draw("same");

   // draw an axis on the right side
   TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
   gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetLineColor(color);
   axis->SetTextColor(color);
   axis->Draw();
   return;
}



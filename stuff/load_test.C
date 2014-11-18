load_test_idx(int i) {

  if (i==0) {  gSystem->Load("libAnalysisTools"); }

  if (i==1) {  gSystem->Load("libBase"); } // OK

  if (i==2) {  gSystem->Load("libDpmEvtGen"); } // OK

  if (i==3) {  gSystem->Load("libDrc"); }

  if (i==4) {  gSystem->Load("libDsk"); }

  if (i==5) {  gSystem->Load("libEMFFGEN"); }

  if (i==6) {  gSystem->Load("libEmc"); }

  if (i==7) {  gSystem->Load("libEventDisplay"); } // OK

  if (i==8) {  gSystem->Load("libEvtGen"); } // OK

  if (i==9) {  gSystem->Load("libEvtGenDirect"); }

  if (i==10) {  gSystem->Load("libEvtGenExternal"); } //OK

  if (i==11) {  gSystem->Load("libFairDB"); } // OK

  if (i==12) {  gSystem->Load("libFairTools"); } // OK

  if (i==13) {  gSystem->Load("libField"); } // OK

  if (i==14) {  gSystem->Load("libFlukaResults"); } //OK

  if (i==15) {  gSystem->Load("libFtof"); }

  if (i==16) {  gSystem->Load("libFts"); }

  if (i==17) {  gSystem->Load("libGeane"); } // OK

  if (i==18) {  gSystem->Load("libGem"); }

  if (i==19) {  gSystem->Load("libGen"); } //OK

  if (i==10) {  gSystem->Load("libGeoBase"); } // OK

  if (i==21) {  gSystem->Load("libGlobal"); }

  if (i==22) {  gSystem->Load("libHyp"); }

  if (i==23) {  gSystem->Load("libLmd"); }

  if (i==24) {  gSystem->Load("libLmdFit"); }

  if (i==25) {  gSystem->Load("libLmdReco"); }

  if (i==26) {  gSystem->Load("libLmdTool"); }

  if (i==27) {  gSystem->Load("libLmdTrk"); }

  if (i==28) {  gSystem->Load("libMCMatch"); }

  if (i==29) {  gSystem->Load("libMCMatchExamples"); }

  if (i==30) {  gSystem->Load("libMbsAPI"); } // OK

  if (i==31) {  gSystem->Load("libMdt"); }

  if (i==32) {  gSystem->Load("libMva"); }

  if (i==33) {  gSystem->Load("libMvd"); }

  if (i==34) {  gSystem->Load("libMvdReco"); }

  if (i==35) {  gSystem->Load("libMvdTrk"); }

  if (i==36) {  gSystem->Load("libPGen"); } // OK

  if (i==37) {  gSystem->Load("libParBase"); } // OK

  if (i==38) {  gSystem->Load("libPassive"); } // OK

  if (i==39) {  gSystem->Load("libPhotosCxxInterface"); } // OK

  if (i==40) {  gSystem->Load("libPhotosFortran"); } // OK

  if (i==41) {  gSystem->Load("libPid"); }

  if (i==42) {  gSystem->Load("libPndData"); }

  if (i==43) {  gSystem->Load("libPndEventDisplay"); }

  if (i==44) {  gSystem->Load("libRecoHits"); }

  if (i==45) {  gSystem->Load("libRecoTasks"); }

  if (i==46) {  gSystem->Load("libRho"); }

  if (i==47) {  gSystem->Load("libRich"); }

  if (i==48) {  gSystem->Load("libSciT"); }

  if (i==49) {  gSystem->Load("libSds"); }

  if (i==50) {  gSystem->Load("libSdsReco"); }

  if (i==51) {  gSystem->Load("libSofTrig"); }

  if (i==52) {  gSystem->Load("libStt"); }

  if (i==53) {  gSystem->Load("libSttCellTrackFinder"); }

  if (i==54) {  gSystem->Load("libSttMvdTracking"); }

  if (i==55) {  gSystem->Load("libTracking"); }

  if (i==56) {  gSystem->Load("libTrackingQA"); }

  if (i==57) {  gSystem->Load("libTrkBase"); } // OK

  if (i==58) {  gSystem->Load("libbuffers"); }

  if (i==59) {  gSystem->Load("libfsim"); }

  if (i==60) {  gSystem->Load("libgeneralTools"); }

  if (i==61) {  gSystem->Load("libgenfit"); } // OK

  if (i==62) {  gSystem->Load("libgenfitAdapters"); }

  if (i==63) {  gSystem->Load("libhelixpropagator"); }

  if (i==64) {  gSystem->Load("libriemann"); }

  if (i==65) {  gSystem->Load("libtrackrep"); }

}

void load_test() {
  for (int i=0; i<66; ++i ) {
    cout << "Testing loadin of modlue " << i << endl;
    load_test_idx(i);
  }
}

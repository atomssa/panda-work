#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class  PndEmc+;
#pragma link C++ class  PndEmcPoint+;
#pragma link C++ class  PndGeoEmc+;
#pragma link C++ class  PndEmcGeoPar+;
#pragma link C++ class  PndEmcApd+;
#pragma link C++ class  PndEmcApdPoint+;
#pragma link C++ class  PndGeoEmcApd+;
#pragma link C++ class  PndEmcContFact+;
#pragma link C++ class  PndEmcHit+;
#pragma link C++ class  PndEmcApdHit+;
#pragma link C++ class  PndEmcHitProducer+;
#pragma link C++ class  PndEmcApdHitProducer+;
#pragma link C++ class  PndEmcTwoCoordIndex+;
#pragma link C++ class  PndEmcXtal+;
#pragma link C++ class  PndEmcMapper+;
#pragma link C++ class  PndEmcStructure+;
#pragma link C++ class  PndEmcWaveformData+;
#pragma link C++ class  PndEmcWaveform+;
#pragma link C++ class  PndEmcMultiWaveform+;
#pragma link C++ class  PndEmcDigi+;
#pragma link C++ class  PndEmcSharedDigi;
#pragma link C++ class  PndEmcHitsToWaveform+;
#pragma link C++ class  PndEmcFWEndcapTimebasedWaveforms+;
#pragma link C++ class  PndEmcAbsWaveformSimulator+;
#pragma link C++ class  PndEmcMultiWaveformSimulator+;
#pragma link C++ class  PndEmcAbsWaveformModifier+;
#pragma link C++ class  PndEmcWaveformDigitizer+;
#pragma link C++ class  PndEmcFullStackedWaveformSimulator+;
#pragma link C++ class  PndEmcShapingNoiseAdder+;
#pragma link C++ class  PndEmcWaveformToDigi+;
#pragma link C++ class  PndEmcWaveformToCalibratedDigi+;
#pragma link C++ class  PndEmcMultiWaveformToCalibratedDigi+;
#pragma link C++ class  PndEmcFadcFilter+;
#pragma link C++ class  PndEmcAsicPulseshape+;
#pragma link C++ class  PndEmcAbsPulseshape+;
#pragma link C++ class  PndEmcFittedPulseshape+;
#pragma link C++ class  PndEmcFullDigiTask+;
#pragma link C++ class  PndEmcDigiPar+;
#pragma link C++ class  PndEmcDigiNonuniformityPar+;
#pragma link C++ class  PndEmcDigiNonuniParObject+;
#pragma link C++ class  PndEmcCluster+;
#pragma link C++ class  PndEmcMakeDigi+;
#pragma link C++ class  PndEmcMakeCluster+;
#pragma link C++ class  PndEmcHeader+;
#pragma link C++ class  PndEmcHdrFiller+;
#pragma link C++ class  PndEmcRecoPar+;
#pragma link C++ class  PndEmcBump+;
#pragma link C++ class  PndEmcMakeBump+;
#pragma link C++ class  PndEmc2DLocMaxFinder+;
#pragma link C++ class  PndEmcExpClusterSplitter+;
#pragma link C++ class  PndEmcPhiBumpSplitter+;
#pragma link C++ class  PndEmcRecoHit+;
#pragma link C++ class  PndEmcMakeRecoHit+;
#pragma link C++ class  PndEmcClusterProperties+;
#pragma link C++ class  PndEmcClusterDistances+;
#pragma link C++ class  PndEmcClusterEnergySums+;
#pragma link C++ class  PndEmcClusterMoments+;
#pragma link C++ class  PndEmcXClMoments+;
#pragma link C++ class  PndEmcMakeCorr+;
#pragma link C++ class  PndEmcCorrection+;
#pragma link C++ class  PndEmcErrorMatrix+;
#pragma link C++ class  PndEmcErrorMatrixPar+;
#pragma link C++ class  PndEmcErrorMatrixParObject+;
#pragma link C++ class  PndEmcClusterCalibrator+;
#pragma link C++ class  PndEmcAbsClusterCalibrator+;
#pragma link C++ class  PndEmcClusterHistCalibrator+;
#pragma link C++ class  PndEmcClusterSimpleCalibrator+;
#pragma link C++ class  PndEmcClusterCalibrationPar+;
#pragma link C++ class  PndEmcClusterCalibrationParObject+;

#pragma link C++ class  PndEmcWaveformBuffer+;
#pragma link C++ class  PndEmcDigiWriteoutBuffer+;
#pragma link C++ class  PndEmcDigiRingSorter+;
#pragma link C++ class  PndEmcDigiSorterTask+;

//#pragma link C++ class  ReadMainzProto60+;
//#pragma link C++ class  ReadMainzProto60v4+;
//#pragma link C++ class  ReadMainzProto60v6+;
//#pragma link C++ class  PndEmcReadProtoData+;
//
//#pragma link C++ class  PndEmcReadProto192Data+;
//#pragma link C++ class  TProtoUnpackEvent+;
//#pragma link C++ class  TGo4EventElement+;

#pragma link C++ class PndEmcAbsPSA+;
#pragma link C++ class PndEmcPSAParabolicBaseline+;
#pragma link C++ class PndEmcPSAFPGADigitalFilterAnalyser+;
#pragma link C++ class PndEmcPSAFPGAFilterCF+;
#pragma link C++ class PndEmcPSAFPGAFilterDelay+;
#pragma link C++ class PndEmcPSAFPGAFilterLine+;
#pragma link C++ class PndEmcPSAFPGAFilterMA+;
#pragma link C++ class PndEmcPSAFPGAFilterMWD+;
#pragma link C++ class PndEmcPSAFPGAIntegratingAnalyser+;
#pragma link C++ class PndEmcPSAFPGALinFitter+;
#pragma link C++ class PndEmcPSAFPGAMLinFitter+;
#pragma link C++ class PndEmcPSAFPGASampleAnalyser+;


#pragma link C++ class PndEmcAbsCrystalCalibrator+;
#pragma link C++ class PndEmcDummyCrystalCalibrator+;
#pragma link C++ class PndEmcFileCrystalCalibrator+;
#pragma link C++ class PndEmcSimCrystalCalibrator+;

#pragma link C++ class  vector<PndEmcHit*>;
#pragma link C++ class  vector<PndEmcPoint*>;
#pragma link C++ class  vector<PndEmcDigi*>;

#pragma link C++ class PndEmcWaveformWriteoutBuffer+;
#pragma link C++ class PndEmcWaveformSorterTask+;
#pragma link C++ class PndEmcWaveformRingSorter+;
#pragma link C++ class PndEmcFpgaPar+;
#pragma link C++ class PndEmcDigiCalibrator+;
#pragma link C++ class PndEmcCorrBump+;
#pragma link C++ class PndEmcAnalysis+;

#endif

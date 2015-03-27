# Create a library called "" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR}
${BASE_INCLUDE_DIRECTORIES}
${SIMPATH}/include
${CMAKE_SOURCE_DIR}/passive
${CMAKE_SOURCE_DIR}/field
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/emc/EmcData
${CMAKE_SOURCE_DIR}/emc/EmcDigi
${CMAKE_SOURCE_DIR}/emc/EmcMC
${CMAKE_SOURCE_DIR}/emc/EmcReco
${CMAKE_SOURCE_DIR}/emc/EmcTools
#${CMAKE_SOURCE_DIR}/emc/EmcProtoTools
${CMAKE_SOURCE_DIR}/emc/EmcCorr
${CMAKE_SOURCE_DIR}/emc/EmcTimeBased
${CMAKE_SOURCE_DIR}/emc
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set(SRCS
EmcData/PndEmcPoint.cxx
EmcData/PndEmcApdPoint.cxx
EmcData/PndEmcHit.cxx
EmcData/PndEmcApdHit.cxx
EmcData/PndEmcWaveformData.cxx
EmcData/PndEmcWaveform.cxx
EmcData/PndEmcMultiWaveform.cxx
EmcData/PndEmcDigi.cxx
EmcData/PndEmcCluster.cxx
EmcData/PndEmcSharedDigi.cxx
EmcData/PndEmcBump.cxx
EmcData/PndEmcHeader.cxx
EmcData/PndEmcRecoHit.cxx
EmcData/PndEmcClusterProperties.cxx
EmcData/PndEmcAbsClusterProperty.cxx
EmcData/PndEmcClusterDistances.cxx
EmcData/PndEmcClusterEnergySums.cxx
EmcData/PndEmcClusterMoments.cxx
EmcData/PndEmcXClMoments.cxx
EmcDigi/PndEmcHitProducer.cxx
EmcDigi/PndEmcApdHitProducer.cxx
EmcDigi/PndEmcCRRCPulseshape.cxx
EmcDigi/PndEmcCR2RCPulseshape.cxx
EmcDigi/PndEmcExpPulseshape.cxx
EmcDigi/PndEmcAsicPulseshape.cxx
EmcDigi/PndEmcFittedPulseshape.cxx
EmcDigi/PndEmcHitsToWaveform.cxx
#EmcDigi/PndEmcHitsToTmpWaveform.cxx
EmcDigi/PndEmcFullStackedWaveformSimulator.cxx
EmcDigi/PndEmcMultiWaveformSimulator.cxx
EmcDigi/PndEmcWaveformDigitizer.cxx
EmcDigi/PndEmcShapingNoiseAdder.cxx
EmcDigi/PndEmcFWEndcapTimebasedWaveforms.cxx
EmcDigi/PndEmcWaveformToDigi.cxx
EmcDigi/PndEmcWaveformToCalibratedDigi.cxx
EmcDigi/PndEmcMultiWaveformToCalibratedDigi.cxx
#EmcDigi/PndEmcTmpWaveformToDigi.cxx
EmcDigi/PndEmcDigiPar.cxx
EmcDigi/PndEmcDigiNonuniformityPar.cxx
EmcDigi/PndEmcMakeDigi.cxx
EmcDigi/PndEmcAbsPSA.cxx
EmcDigi/PndEmcPSAParabolic.cxx
EmcDigi/PndEmcPSAParabolicBaseline.cxx
EmcDigi/PndEmcPSATrapDigiFilter.cxx
EmcDigi/PndEmcPSAMatchedDigiFilter.cxx
EmcDigi/PndEmcFadcFilter.cxx
EmcDigi/PndEmcFullDigiTask.cxx
EmcDigi/PndEmcWaveformBuffer.cxx
EmcDigi/PndEmcDigiWriteoutBuffer.cxx
EmcDigi/PndEmcDigiRingSorter.cxx
EmcDigi/PndEmcDigiSorterTask.cxx
EmcMC/PndGeoEmc.cxx
EmcMC/PndEmc.cxx
EmcMC/PndEmcReader.cxx
EmcMC/PndGeoEmcApd.cxx
EmcMC/PndEmcApd.cxx
EmcMC/PndEmcGeoPar.cxx
EmcTools/PndEmcTwoCoordIndex.cxx
EmcTools/PndEmcStructure.cxx
EmcTools/PndEmcMapper.cxx
EmcTools/PndEmcXtal.cxx
EmcTools/PndEmcHdrFiller.cxx
#EmcProtoTools/ReadMainzProto60.cxx
#EmcProtoTools/ReadMainzProto60v4.cxx
#EmcProtoTools/ReadMainzProto60v6.cxx
#EmcProtoTools/PndEmcReadProtoData.cxx
#EmcProtoTools/PndEmcReadProto192Data.cxx
#EmcProtoTools/TProtoUnpackEvent.cxx
#EmcProtoTools/TGo4EventElement.cxx
EmcReco/PndEmcMakeCluster.cxx
EmcReco/PndEmcRecoPar.cxx
EmcReco/PndEmcMakeBump.cxx
EmcReco/PndEmc2DLocMaxFinder.cxx
EmcReco/PndEmcExpClusterSplitter.cxx
EmcReco/PndEmcMakeRecoHit.cxx
EmcReco/PndEmcErrorMatrix.cxx
EmcReco/PndEmcErrorMatrixPar.cxx
EmcReco/PndEmcPhiBumpSplitter.cxx
EmcCorr/PndEmcMakeCorr.cxx
EmcCorr/PndEmcCorrection.cxx
EmcCorr/PndEmcClusterCalibrator.cxx
EmcCorr/PndEmcClusterCalibrationPar.cxx
#basetasks/PndGlobalEvtTask.cxx
PndEmcContFact.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGADigitalFilterAnalyser.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGAFilterCF.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGAFilterDelay.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGAFilterLine.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGAFilterMA.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGAFilterMWD.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGAIntegratingAnalyser.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGALinFitter.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGAMLinFitter.cxx
EmcDigi/PndEmcPSAFPGA/PndEmcPSAFPGASampleAnalyser.cxx
EmcDigi/PndEmcAbsCrystalCalibrator.cxx
EmcDigi/PndEmcDummyCrystalCalibrator.cxx
EmcDigi/PndEmcFileCrystalCalibrator.cxx
EmcDigi/PndEmcSimCrystalCalibrator.cxx
EmcTimeBased/PndEmcWaveformWriteoutBuffer.cxx
EmcTimeBased/PndEmcWaveformRingSorter.cxx
EmcTimeBased/PndEmcWaveformSorterTask.cxx
EmcTimeBased/PndEmcFpgaPar.cxx
EmcTimeBased/PndEmcDigiCalibrator.cxx
EmcTimeBased/PndEmcCorrBump.cxx
EmcTimeBased/PndEmcAnalysis.cxx
)

set(LINKDEF  EmcLinkDef.h)
set(LIBRARY_NAME Emc)

set(DEPENDENCIES Base GeoBase ParBase PndData Passive)

GENERATE_LIBRARY()

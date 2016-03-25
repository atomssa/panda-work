# Create a library called "libPid"

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${SYSTEM_INCLUDE_DIRECTORIES}
  ${GEANT3_INCLUDE_DIR}
  ${BASE_INCLUDE_DIRECTORIES}
)

Set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/pndbase/PndStdUtils
  ${CMAKE_SOURCE_DIR}/field
  ${CMAKE_SOURCE_DIR}/trackbase
  ${CMAKE_SOURCE_DIR}/passive
  ${CMAKE_SOURCE_DIR}/PndTools/generalTools
  ${CMAKE_SOURCE_DIR}/emc/EmcMC
  ${CMAKE_SOURCE_DIR}/emc/EmcData
  ${CMAKE_SOURCE_DIR}/emc/EmcReco
  ${CMAKE_SOURCE_DIR}/emc/EmcTools
  ${CMAKE_SOURCE_DIR}/emc/EmcCorr
  ${CMAKE_SOURCE_DIR}/gem
  ${CMAKE_SOURCE_DIR}/genfit
  ${CMAKE_SOURCE_DIR}/pnddata
  ${CMAKE_SOURCE_DIR}/pnddata/SttData
  ${CMAKE_SOURCE_DIR}/pnddata/FtsData
  ${CMAKE_SOURCE_DIR}/pnddata/SciTData
  ${CMAKE_SOURCE_DIR}/pnddata/MdtData
  ${CMAKE_SOURCE_DIR}/pnddata/DrcData
  ${CMAKE_SOURCE_DIR}/pnddata/DskData
  ${CMAKE_SOURCE_DIR}/pnddata/RichData
  ${CMAKE_SOURCE_DIR}/pnddata/SdsData
  ${CMAKE_SOURCE_DIR}/pnddata/GemData
  ${CMAKE_SOURCE_DIR}/pnddata/TrackData
  ${CMAKE_SOURCE_DIR}/pnddata/PidData
  ${CMAKE_SOURCE_DIR}/ftof
  ${CMAKE_SOURCE_DIR}/stt
  ${CMAKE_SOURCE_DIR}/stt/sttreco
  ${CMAKE_SOURCE_DIR}/fts
  ${CMAKE_SOURCE_DIR}/geane
  ${CMAKE_SOURCE_DIR}/trackrep
  ${CMAKE_SOURCE_DIR}/rho
  ${CMAKE_SOURCE_DIR}/rho/RhoBase
  ${CMAKE_SOURCE_DIR}/rho/RhoTools
  ${CMAKE_SOURCE_DIR}/pid
  ${CMAKE_SOURCE_DIR}/pid/PidData
  ${CMAKE_SOURCE_DIR}/pid/PidCorr
  ${CMAKE_SOURCE_DIR}/pid/PidClassifier
  ${CMAKE_SOURCE_DIR}/GenfitTools/recotasks
  ${CMAKE_SOURCE_DIR}/PndTools/MVA
  ${CMAKE_SOURCE_DIR}/PndTools/MVA/KNN
  ${CMAKE_SOURCE_DIR}/PndTools/MVA/LVQ
  ${CMAKE_SOURCE_DIR}/PndTools/MVA/TMVAkd_KNN
  ${CMAKE_SOURCE_DIR}/PndTools/MVA/TMVA_MCL
# ${CMAKE_SOURCE_DIR}/PndTools/MVA/TMVA_PDE_FOAM
  ${CMAKE_SOURCE_DIR}/rich
)

Include_Directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set(SRCS
PidCorr/PndPidCorrPar.cxx
PidCorr/PndPidCorrelator.cxx
PidCorr/PndPidTrackInfo.cxx
PidCorr/PndPidMvdInfo.cxx
PidCorr/PndPidSttInfo.cxx
PidCorr/PndPidGemInfo.cxx
PidCorr/PndPidFtsInfo.cxx
PidCorr/PndPidFscInfo.cxx
PidCorr/PndPidFtofInfo.cxx
PidCorr/PndPidMdtInfo.cxx
PidCorr/PndPidDskInfo.cxx
PidCorr/PndPidDrcInfo.cxx
PidCorr/PndPidRichInfo.cxx
PidCorr/PndPidEmcInfo.cxx
PidCorr/PndPidTofInfo.cxx
PidCorr/PndPidBremCorrector.cxx
PidCorr/PndPidBremCorrectorNT.cxx
PidCorr/PndPidBremCorrectorSaveMCT.cxx
PndPidContFact.cxx

###### Associator
PidClassifier/PndPidIdealAssociatorTask.cxx
PidClassifier/PndPidIdealFwdOnlyAssociatorTask.cxx
PidClassifier/PndPidMvdAssociatorTask.cxx
PidClassifier/PndPidDrcAssociatorTask.cxx
PidClassifier/PndPidDiscAssociatorTask.cxx
PidClassifier/PndPidRichAssociatorTask.cxx
PidClassifier/PndPidMdtHCAssociatorTask.cxx
PidClassifier/PndPidMvaAssociatorTask.cxx
PidClassifier/PndPidEmcAssociatorTask.cxx
PidClassifier/PndPidSttAssociatorTask.cxx
PidClassifier/PndPidEmcBayesAssociatorTask.cxx

###### Parameters
PidClassifier/PndPidMvdPar.cxx
)

set(HEADERS
PidCorr/PndPidCorrPar.h
PidCorr/PndPidCorrelator.h
PidCorr/PndPidBremCorrector.h
PidCorr/PndPidBremCorrectorNT.h
PidCorr/PndPidBremCorrectorSaveMCT.h
PndPidContFact.h

###### Associator Headers
PidClassifier/PndPidIdealAssociatorTask.h
PidClassifier/PndPidIdealFwdOnlyAssociatorTask.h
PidClassifier/PndPidMvdAssociatorTask.h
PidClassifier/PndPidDrcAssociatorTask.h
PidClassifier/PndPidDiscAssociatorTask.h
PidClassifier/PndPidRichAssociatorTask.h
PidClassifier/PndPidMdtHCAssociatorTask.h
PidClassifier/PndPidMvaAssociatorTask.h
PidClassifier/PndPidEmcAssociatorTask.h
PidClassifier/PndPidSttAssociatorTask.h
PidClassifier/PndPidEmcBayesAssociatorTask.h

###### Parameters
PidClassifier/PndPidMvdPar.h
)

set(LINKDEF PidLinkDef.h)
set(LIBRARY_NAME Pid)

set(DEPENDENCIES Base GeoBase ParBase RecoTasks generalTools Mva Emc
Ftof Fts Geane PndData Stt TrkBase Rho)

GENERATE_LIBRARY()

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR}
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_SOURCE_DIR}/base/event
${CMAKE_SOURCE_DIR}/field
${CMAKE_SOURCE_DIR}/trackbase
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/pnddata/PidData
${CMAKE_SOURCE_DIR}/passive
${CMAKE_SOURCE_DIR}/tpc
${CMAKE_SOURCE_DIR}/lhetrack

${CMAKE_SOURCE_DIR}/emc/EmcData
#${CMAKE_SOURCE_DIR}/emc/EmcDigi
#${CMAKE_SOURCE_DIR}/emc/EmcMC
${CMAKE_SOURCE_DIR}/emc/EmcReco
${CMAKE_SOURCE_DIR}/emc/EmcTools
#${CMAKE_SOURCE_DIR}/emc/EmcProtoTools
#${CMAKE_SOURCE_DIR}/emc/EmcCorr
#${CMAKE_SOURCE_DIR}/emc/EmcTimeBased
#${CMAKE_SOURCE_DIR}/emc

${CMAKE_SOURCE_DIR}/fsim
${CMAKE_SOURCE_DIR}/rho
${CMAKE_SOURCE_DIR}/rho/RhoBase
${CMAKE_SOURCE_DIR}/rho/RhoSelector
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools/Fitter
${CMAKE_SOURCE_DIR}/genfit
${CMAKE_SOURCE_DIR}/src/filler
${CMAKE_SOURCE_DIR}/src/AnaTda
${CMAKE_SOURCE_DIR}/src/BremPidReader
)

include_directories( ${INCLUDE_DIRECTORIES})

add_definitions(${CXX11_COMPILER_FLAGS})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set(SRCS BremPidReader.cxx)

set(LINKDEF BremPidReaderLinkDef.h)
set(LIBRARY_NAME BremPidReader)

#set(DEPENDENCIES Base PndData Rho AnalysisTools)

GENERATE_LIBRARY()
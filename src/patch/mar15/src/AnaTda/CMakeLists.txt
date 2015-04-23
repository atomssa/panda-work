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
${CMAKE_SOURCE_DIR}/emc/EmcReco
${CMAKE_SOURCE_DIR}/fsim
${CMAKE_SOURCE_DIR}/rho
${CMAKE_SOURCE_DIR}/rho/RhoBase
${CMAKE_SOURCE_DIR}/rho/RhoSelector
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools/Fitter
${CMAKE_SOURCE_DIR}/genfit
${CMAKE_SOURCE_DIR}/src/filler
${CMAKE_SOURCE_DIR}/src/AnaTda
)

include_directories( ${INCLUDE_DIRECTORIES})

add_definitions(${CXX11_COMPILER_FLAGS})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES} ${FAIRROOT_LIBRARY_DIR})

set(SRCS
AnaTda.cxx
AnaTdav2.cxx
)

set(LINKDEF AnaTdaLinkDef.h)
set(LIBRARY_NAME anatda)

#set(DEPENDENCIES Base PndData Rho AnalysisTools TrkBase)
#set(DEPENDENCIES PndData Rho AnalysisTools)


GENERATE_LIBRARY()
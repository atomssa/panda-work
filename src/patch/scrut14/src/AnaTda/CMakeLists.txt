# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Create a library called "anatda"
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
${CMAKE_SOURCE_DIR}/src/AnaTda
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set(anatda_SRCS
AnaTda.cxx
)

set(anatda_HEADERS
AnaTda.h
)

set(DEPENDENCIES Base PndData Rho AnalysisTools TrkBase)

set(anatda_LINKDEF AnaTdaLinkDef.h)
set(anatda_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/AnaTda_Dict.cxx)

ROOT_GENERATE_DICTIONARY("${anatda_HEADERS}" "${anatda_LINKDEF}" "${anatda_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

set(anatda_SRCS ${anatda_SRCS} ${anatda_DICTIONARY})

add_library(anatda SHARED ${anatda_SRCS})
target_link_libraries(anatda ${ROOT_LIBRARIES})
set_target_properties(anatda PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES})

################ install ###################
install(TARGETS anatda DESTINATION ${CMAKE_BINARY_DIR}/lib)

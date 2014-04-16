# Create a library called "histfillertask" 

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
${BASE_INCLUDE_DIRECTORIES} 
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/emc/EmcData
${CMAKE_SOURCE_DIR}/emc/EmcReco
${CMAKE_SOURCE_DIR}/emc/EmcTools
${CMAKE_SOURCE_DIR}/macro/PhysListHadrG4
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(histfillertask_SRCS
HistFillerTask.cxx
)

set(histfillertask_HEADERS
HistFillerTask.h
)

set(histfillertask_LINKDEF HistFillerLinkDef.h)
set(histfillertask_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/HistFiller_Dict.cxx) 

ROOT_GENERATE_DICTIONARY("${histfillertask_HEADERS}" "${histfillertask_LINKDEF}" "${histfillertask_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

set(histfillertask_SRCS ${histfillertask_SRCS} ${histfillertask_DICTIONARY})

add_library(histfillertask SHARED ${histfillertask_SRCS})
target_link_libraries(histfillertask ${ROOT_LIBRARIES})
set_target_properties(histfillertask PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES})
################ install ###################
install(TARGETS histfillertask DESTINATION ${CMAKE_BINARY_DIR}/lib)

# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Create a library called "threepion" 
set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
${BASE_INCLUDE_DIRECTORIES} 
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/base/event
${CMAKE_SOURCE_DIR}/macro/ThreePions
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(threepions_SRCS
ThreePions.cxx
)

set(threepions_HEADERS
ThreePions.h
)

set(threepions_LINKDEF ThreePionsLinkDef.h)
set(threepions_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/ThreePions_Dict.cxx) 

ROOT_GENERATE_DICTIONARY("${threepions_HEADERS}" "${threepions_LINKDEF}" "${threepions_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

set(threepions_SRCS ${threepions_SRCS} ${threepions_DICTIONARY})


add_library(threepions SHARED ${threepions_SRCS})
target_link_libraries(threepions ${ROOT_LIBRARIES})
set_target_properties(threepions PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES})


################ install ###################
install(TARGETS threepions DESTINATION ${CMAKE_BINARY_DIR}/lib)

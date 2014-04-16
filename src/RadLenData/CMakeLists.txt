# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Create a library called "radlendata" 
set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
${BASE_INCLUDE_DIRECTORIES} 
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/base/event
${CMAKE_SOURCE_DIR}/macro/RadLenData
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(radlendata_SRCS
RadLenData.cxx
)

set(radlendata_HEADERS
RadLenData.h
)

set(radlendata_LINKDEF RadLenDataLinkDef.h)
set(radlendata_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/RadLenData_Dict.cxx) 

ROOT_GENERATE_DICTIONARY("${radlendata_HEADERS}" "${radlendata_LINKDEF}" "${radlendata_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

set(radlendata_SRCS ${radlendata_SRCS} ${radlendata_DICTIONARY})


add_library(radlendata SHARED ${radlendata_SRCS})
target_link_libraries(radlendata ${ROOT_LIBRARIES})
set_target_properties(radlendata PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES})


################ install ###################
install(TARGETS radlendata DESTINATION ${CMAKE_BINARY_DIR}/lib)

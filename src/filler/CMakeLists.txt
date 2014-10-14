include_directories(${ROOT_INCLUDE_DIR})

link_directories(${LINK_DIRECTORIES})

set(SRCS filler.cxx)

set(LIBRARY_NAME filler)

add_definitions(${CXX11_COMPILER_FLAGS})

#set(DEPENDENCIES )

GENERATE_LIBRARY()

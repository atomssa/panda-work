# Create a library called "libEvtGenStandAlone" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
   ${HEPMC_INCLUDE_DIR}
   ${ROOT_INCLUDE_DIR}
   ${BASE_INCLUDE_DIRECTORIES}
   ${EVTGEN_INCLUDE_DIR}
   ${CMAKE_SOURCE_DIR}/pgenerators/EvtGen/EvtGenStandAlone
)

if(EVTGEN_EXTERNAL)
   add_definitions("-DEVTGEN_EXTERNAL")
endif(EVTGEN_EXTERNAL)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
   ${ROOT_LIBRARY_DIR}
   ${HEPMC_LIB_DIR}
   ${FORTRAN_LIBDIR}
)

link_directories( ${LINK_DIRECTORIES})

set(SRCS
   PndEvtGenStandAlone.cxx
)

set(LINKDEF  EvtGenStandAloneLinkDef.h)
set(LIBRARY_NAME EvtGenSA)

set(DEPENDENCIES Base EvtGen EvtGenExternal)

GENERATE_LIBRARY()

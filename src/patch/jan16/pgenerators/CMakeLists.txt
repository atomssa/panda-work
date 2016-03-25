# Create a library called "libPGen" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${SYSTEM_INCLUDE_DIRECTORIES}
  ${CLHEP_INCLUDE_DIR}
  ${GEANT4_INCLUDE_DIR}
  ${PYTHIA8_INCLUDE_DIR}
  ${BASE_INCLUDE_DIRECTORIES}
)

Set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/pgenerators
  ${CMAKE_SOURCE_DIR}/pgenerators/FtfEvtGen
  ${CMAKE_SOURCE_DIR}/pgenerators/FtfEvtGen/include
  ${CMAKE_SOURCE_DIR}/pgenerators/eventFilter
  ${CMAKE_SOURCE_DIR}/pnddata
  ${CMAKE_SOURCE_DIR}/rho
  ${CMAKE_SOURCE_DIR}/rho/RhoBase
  ${CMAKE_SOURCE_DIR}/rho/RhoSelector
  ${CMAKE_SOURCE_DIR}/fairtools
)

Include_Directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
${GEANT4_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set(SRCS
PndRestGas.cxx
PndDpmGenerator.cxx
PndDpmGeneratorMod.cxx
PndFtfGenerator.cxx
PndEvtGenGenerator.cxx
PndDpmDirect.cxx
#PndFtfDirect.cxx
PndBoxGenerator.cxx
PndUrqmdSmmGenerator.cxx
PndJRGenerator.cxx
PndPythia6Generator.cxx
PndPythia8Generator.cxx
PndFlatParticleGenerator.cxx
PndFlukaGenerator.cxx
PndGasGenerator.cxx
PndHypBupGenerator.cxx
PndGiBuuGenerator.cxx
PndGiBuuTrack.cxx
PndGiBuuDecayerAbs.cxx
PndGiBuuJPsiEE.cxx
eventFilter/FairFilteredPrimaryGenerator.cxx
eventFilter/FairEvtFilterParams.cxx
eventFilter/FairEvtFilter.cxx
eventFilter/FairEvtFilterOnSingleParticleCounts.cxx
eventFilter/PndEvtFilter.cxx
eventFilter/PndEvtFilterOnInvMassCounts.cxx
)

set(LINKDEF  PGenLinkDef.h)
set(LIBRARY_NAME PGen)

set(DEPENDENCIES Base DpmEvtGen pythia8)

GENERATE_LIBRARY()



################ FTF DIRECT ###################
#	 	set(SRCS PndFtfDirect.cxx)
#	 	set(LINKDEF  FtfLinkDef.h)
#	 	set(LIBRARY_NAME FtfGenDirect)
#	 	set(DEPENDENCIES G4clhep G4global G4geometry G4physicslists G4materials EG FtfEvtGen)
#	 	#Set(DEFINITIONS G4GLOB_ALLOC_EXPORT=1) # Is it really needed??
#	 	GENERATE_LIBRARY()

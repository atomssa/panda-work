# Create a library called "libPGen" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR}
${PYTHIA8_INCLUDE_DIR}
${SIMPATH}/include
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_SOURCE_DIR}/field
${CMAKE_SOURCE_DIR}/generators
${CMAKE_SOURCE_DIR}/pgenerators
${CMAKE_SOURCE_DIR}/pgenerators/eventFilter
${CMAKE_SOURCE_DIR}/base/sim
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/rho
${CMAKE_SOURCE_DIR}/rho/RhoBase
${CMAKE_SOURCE_DIR}/rho/RhoSelector
${CMAKE_SOURCE_DIR}/fairtools
${CMAKE_SOURCE_DIR}/geobase
${CMAKE_SOURCE_DIR}/parbase
${CMAKE_SOURCE_DIR}/base
${CMAKE_SOURCE_DIR}/base/steer
${CMAKE_SOURCE_DIR}/base/event
${CMAKE_SOURCE_DIR}/base/field
${CMAKE_SOURCE_DIR}/base/sim
${CMAKE_SOURCE_DIR}/base/source
${CMAKE_SOURCE_DIR}/dbase/dbInterface
${CMAKE_SOURCE_DIR}/dbase/dbValidation
${CMAKE_SOURCE_DIR}/dbase/dbUtils
${CMAKE_SOURCE_DIR}/input/db
${CMAKE_SOURCE_DIR}/dbase/dbInput
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set(SRCS
PndRestGas.cxx
PndDpmGenerator.cxx
PndDpmGeneratorMod.cxx
PndFtfGenerator.cxx
PndEvtGenGenerator.cxx
PndDpmDirect.cxx
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

set(DEPENDENCIES Base DpmEvtGen)

GENERATE_LIBRARY()

# Create a library called "libEvtGen" which includes the source files given in
# the array . libEvtGenExternal will be created if enabled by variable

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

set(EVTGEN_PRIVATE_DIR     "Private")
set(EVTGEN_DIR             "R01-03-00")
set(EVTGEN_SRC_DIR         "${EVTGEN_DIR}/src")
set(EVTGEN_EXAMPLE_DIR     "${EVTGEN_DIR}/test")
set(EVTGEN_MODELS_DIR      "${EVTGEN_SRC_DIR}/EvtGenModels")
set(EVTGEN_BASE_DIR        "${EVTGEN_SRC_DIR}/EvtGenBase")
set(EVTGEN_EXTERNAL_DIR    "${EVTGEN_SRC_DIR}/EvtGenExternal")

## important: Private dir has to be BEFORE the original EvtGen in the Path ... just in case we want to use OUR header file

set(INCLUDE_DIRECTORIES
   ${EVTGEN_PRIVATE_DIR}
   ${EVTGEN_DIR}
   ${HEPMC_INCLUDE_DIR}
   ${PYTHIA8_INCLUDE_DIR}
   ${ROOT_INCLUDE_DIR}
   ${PHOTOS_INCLUDE_DIR}
   ${TAUOLA_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
   ${ROOT_LIBRARY_DIR}
   ${HEPMC_LIB_DIR}
   ${PYTHIA8_LIB_DIR}
   ${FORTRAN_LIBDIR}
)

link_directories( ${LINK_DIRECTORIES})

set(EVTGEN_SRCS
   ${EVTGEN_SRC_DIR}/EvtGen.cpp
## PANDA Private
   ${EVTGEN_PRIVATE_DIR}/EvtPpbarGammaGamma.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtPpbarPi0Gamma.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtElectronPositron.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtDeuteronSpectator.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtPpbarEtacGamma.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtLambdaLambdaBar.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtLambdaLambdaBarPol.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtLambdaLambdaBarHE.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtDsDalitz.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtChic1toJpsiGam.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtChic2toJpsiGam.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtDecayNew.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtPPPsiPiPi.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtPPPsiPiPi_matr2_f0.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtPP_EtacGamma.cpp
## overwrite EvtGen ModelReg by Panda one
   ${EVTGEN_PRIVATE_DIR}/EvtModelReg.cpp
##   ${EVTGEN_MODELS_DIR}/EvtModelReg.cpp
   ${EVTGEN_BASE_DIR}/EvtOrthogVector.cpp
   ${EVTGEN_BASE_DIR}/EvtDecayTable.cpp
   ${EVTGEN_BASE_DIR}/EvtRaritaSchwinger.cpp
   ${EVTGEN_BASE_DIR}/EvtDiracSpinor.cpp
   ${EVTGEN_BASE_DIR}/EvtSemiLeptonicBaryonAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtdFunction.cpp
   ${EVTGEN_BASE_DIR}/EvtSecondary.cpp
   ${EVTGEN_BASE_DIR}/EvtdFunctionSingle.cpp
   ${EVTGEN_BASE_DIR}/EvtValError.cpp
   ${EVTGEN_BASE_DIR}/EvtMHelAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtCPUtil.cpp
   ${EVTGEN_BASE_DIR}/EvtVector4R.cpp
   ${EVTGEN_BASE_DIR}/EvtMNode.cpp
   ${EVTGEN_BASE_DIR}/EvtPto3PAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtAbsLineShape.cpp
   ${EVTGEN_BASE_DIR}/EvtDalitzFlatPdf.cpp
   ${EVTGEN_BASE_DIR}/EvtVectorParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtParticleFactory.cpp
   ${EVTGEN_BASE_DIR}/EvtPropGounarisSakurai.cpp
   ${EVTGEN_BASE_DIR}/EvtDiLog.cpp
   ${EVTGEN_BASE_DIR}/EvtCyclic3.cpp
   ${EVTGEN_BASE_DIR}/EvtFlatLineShape.cpp
   ${EVTGEN_BASE_DIR}/EvtPto3PAmpSmpResolution.cpp
   ${EVTGEN_BASE_DIR}/EvtTensorParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtGenKine.cpp
   ${EVTGEN_BASE_DIR}/EvtPDL.cpp
   ${EVTGEN_BASE_DIR}/EvtMParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtSpinDensity.cpp
   ${EVTGEN_BASE_DIR}/EvtRaritaSchwingerParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtConst.cpp
   ${EVTGEN_BASE_DIR}/EvtGammaMatrix.cpp
   ${EVTGEN_BASE_DIR}/EvtBlattWeisskopf.cpp
   ${EVTGEN_BASE_DIR}/EvtSemiLeptonicScalarAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtIncoherentMixing.cpp
   ${EVTGEN_BASE_DIR}/EvtMultiChannelParser.cpp
   ${EVTGEN_BASE_DIR}/EvtVector4C.cpp
   ${EVTGEN_BASE_DIR}/EvtHepMCEvent.cpp
   ${EVTGEN_BASE_DIR}/EvtManyDeltaFuncLineShape.cpp
   ${EVTGEN_BASE_DIR}/EvtPoint1D.cpp
   ${EVTGEN_BASE_DIR}/EvtDecayMode.cpp
   ${EVTGEN_BASE_DIR}/EvtNonresonantAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtDecayAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtDalitzPlot.cpp
   ${EVTGEN_BASE_DIR}/EvtParserXml.cpp
   ${EVTGEN_BASE_DIR}/EvtRandom.cpp
   ${EVTGEN_BASE_DIR}/EvtRadCorr.cpp
   ${EVTGEN_BASE_DIR}/EvtDalitzResPdf.cpp
   ${EVTGEN_BASE_DIR}/EvtTensor4C.cpp
   ${EVTGEN_BASE_DIR}/EvtDiracParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtDalitzPoint.cpp
   ${EVTGEN_BASE_DIR}/EvtMRes.cpp
   ${EVTGEN_BASE_DIR}/EvtPto3PAmpFactory.cpp
   ${EVTGEN_BASE_DIR}/EvtCGCoefSingle.cpp
   ${EVTGEN_BASE_DIR}/EvtRelBreitWignerBarrierFact.cpp
   ${EVTGEN_BASE_DIR}/EvtId.cpp
   ${EVTGEN_BASE_DIR}/EvtTwoBodyKine.cpp
   ${EVTGEN_BASE_DIR}/EvtDecayProb.cpp
   ${EVTGEN_BASE_DIR}/EvtDalitzCoord.cpp
   ${EVTGEN_BASE_DIR}/EvtAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtSemiLeptonicTensorAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtModelAlias.cpp
   ${EVTGEN_BASE_DIR}/EvtParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtAmpSubIndex.cpp
   ${EVTGEN_BASE_DIR}/EvtIntervalFlatPdf.cpp
   ${EVTGEN_BASE_DIR}/EvtDecayBase.cpp
   ${EVTGEN_BASE_DIR}/EvtReport.cpp
   ${EVTGEN_BASE_DIR}/EvtSemiLeptonicAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtLASSAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtTwoBodyVertex.cpp
   ${EVTGEN_BASE_DIR}/EvtRanf.cpp
   ${EVTGEN_BASE_DIR}/Evt3Rank3C.cpp
   ${EVTGEN_BASE_DIR}/EvtHighSpinParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtEvalHelAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtPropFlatte.cpp
   ${EVTGEN_BASE_DIR}/EvtScalarParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtParticleDecayList.cpp
   ${EVTGEN_BASE_DIR}/EvtSpinType.cpp
   ${EVTGEN_BASE_DIR}/EvtModel.cpp
   ${EVTGEN_BASE_DIR}/EvtPropBreitWigner.cpp
   ${EVTGEN_BASE_DIR}/EvtFlatte.cpp
   ${EVTGEN_BASE_DIR}/EvtPropBreitWignerRel.cpp
   ${EVTGEN_BASE_DIR}/EvtSimpleRandomEngine.cpp
   ${EVTGEN_BASE_DIR}/EvtSemiLeptonicVectorAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtMTrivialLS.cpp
   ${EVTGEN_BASE_DIR}/EvtComplex.cpp
   ${EVTGEN_BASE_DIR}/EvtNeutrinoParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtVector3C.cpp
   ${EVTGEN_BASE_DIR}/EvtKine.cpp
   ${EVTGEN_BASE_DIR}/EvtSpinAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtMassAmp.cpp
   ${EVTGEN_BASE_DIR}/EvtMTree.cpp
   ${EVTGEN_BASE_DIR}/EvtAmpIndex.cpp
   ${EVTGEN_BASE_DIR}/EvtStdHep.cpp
   ${EVTGEN_BASE_DIR}/EvtResonance2.cpp
   ${EVTGEN_BASE_DIR}/EvtStringParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtDalitzReso.cpp
   ${EVTGEN_BASE_DIR}/EvtDecayIncoherent.cpp
   ${EVTGEN_BASE_DIR}/EvtPhotonParticle.cpp
   ${EVTGEN_BASE_DIR}/EvtIdSet.cpp
   ${EVTGEN_BASE_DIR}/EvtIntegPdf1D.cpp
   ${EVTGEN_BASE_DIR}/EvtPartProp.cpp
   ${EVTGEN_BASE_DIR}/EvtTensor3C.cpp
   ${EVTGEN_BASE_DIR}/EvtDecayParm.cpp
   ${EVTGEN_BASE_DIR}/EvtSymTable.cpp
   ${EVTGEN_BASE_DIR}/EvtResonance.cpp
   ${EVTGEN_BASE_DIR}/EvtBreitWignerPdf.cpp
   ${EVTGEN_BASE_DIR}/EvtRanFor.cpp
   ${EVTGEN_BASE_DIR}/EvtVector3R.cpp
   ${EVTGEN_BASE_DIR}/EvtParser.cpp
   ${EVTGEN_BASE_DIR}/EvtMBreitWigner.cpp
   ${EVTGEN_BASE_DIR}/EvtParticleDecay.cpp
   ${EVTGEN_BASE_DIR}/EvtExtGeneratorCommandsTable.cpp
   ${EVTGEN_MODELS_DIR}/EvtOmegaDalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtSll.cpp
   ${EVTGEN_MODELS_DIR}/EvtPFermi.cpp
   ${EVTGEN_MODELS_DIR}/EvtBTo4piCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtSSDCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtSTS.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVVNONCPEIGEN.cpp
   ${EVTGEN_MODELS_DIR}/EvtGenericDalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtISGW.cpp
   ${EVTGEN_MODELS_DIR}/EvtHQET.cpp
   ${EVTGEN_MODELS_DIR}/EvtVPHOtoVISR.cpp
   ${EVTGEN_MODELS_DIR}/EvtEtaDalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtPhiDalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgThreeCoeffFcn.cpp
   ${EVTGEN_MODELS_DIR}/EvtDDalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsllUtil.cpp
   ${EVTGEN_MODELS_DIR}/EvtVSSBMixCPT.cpp
   ${EVTGEN_MODELS_DIR}/EvtVub.cpp
   ${EVTGEN_MODELS_DIR}/EvtVubHybrid.cpp
   ${EVTGEN_MODELS_DIR}/EvtbTosllScalarAmp.cpp
   ${EVTGEN_MODELS_DIR}/EvtKKLambdaCFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgTwoCoeffFcn.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgFourCoeffFcn.cpp
   ${EVTGEN_MODELS_DIR}/EvtISGW2.cpp
   ${EVTGEN_MODELS_DIR}/EvtISGW2FF.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgPtrFunction.cpp
   ${EVTGEN_MODELS_DIR}/EvtVVpipi.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgammaFermiUtil.cpp
   ${EVTGEN_MODELS_DIR}/EvtVtoSll.cpp
   ${EVTGEN_MODELS_DIR}/EvtBBScalar.cpp
   ${EVTGEN_MODELS_DIR}/EvtPi0Dalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtPto3P.cpp
   ${EVTGEN_MODELS_DIR}/EvtbTosllVectorAmp.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgammaAliGreub.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoKpiCPiso.cpp
   ${EVTGEN_MODELS_DIR}/EvtVubBLNP.cpp
   ${EVTGEN_MODELS_DIR}/EvtPVVCPLH.cpp
   ${EVTGEN_MODELS_DIR}/EvtY3SToY1SpipiMoxhay.cpp
   ${EVTGEN_MODELS_DIR}/EvtCBTo3piMPP.cpp
   ${EVTGEN_MODELS_DIR}/EvtVubBLNPHybrid.cpp
   ${EVTGEN_MODELS_DIR}/EvtLb2Lll.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVPCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtSTSCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtBToVlnuBallFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtVubNLO.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgFunction.cpp
   ${EVTGEN_MODELS_DIR}/EvtVVP.cpp
   ${EVTGEN_MODELS_DIR}/EvtBcBsStarNPi.cpp
   ${EVTGEN_MODELS_DIR}/EvtSSSCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgamma.cpp
   ${EVTGEN_MODELS_DIR}/EvtVSS.cpp
   ${EVTGEN_MODELS_DIR}/EvtGoityRoberts.cpp
   ${EVTGEN_MODELS_DIR}/EvtBToPlnuBK.cpp
   ${EVTGEN_MODELS_DIR}/EvtVSPPwave.cpp
   ${EVTGEN_MODELS_DIR}/EvtMelikhov.cpp
   ${EVTGEN_MODELS_DIR}/EvtWilsonCoefficients.cpp
   ${EVTGEN_MODELS_DIR}/EvtPhsp.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVVHelCPMix.cpp
   ${EVTGEN_MODELS_DIR}/EvtTVSPwave.cpp
   ${EVTGEN_MODELS_DIR}/EvtSLBKPole.cpp
   ${EVTGEN_MODELS_DIR}/EvtHypNonLepton.cpp
   ${EVTGEN_MODELS_DIR}/EvtBToDDalitzCPK.cpp
   ${EVTGEN_MODELS_DIR}/EvtSSSCPT.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgammaKagan.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsll.cpp
   ${EVTGEN_MODELS_DIR}/EvtBcBsNPi.cpp
   ${EVTGEN_MODELS_DIR}/EvtSLBKPoleFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVSCPLH.cpp
   ${EVTGEN_MODELS_DIR}/EvtSingleParticle.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoKD3P.cpp
   ${EVTGEN_MODELS_DIR}/EvtVll.cpp
   ${EVTGEN_MODELS_DIR}/EvtVPHOtoVISRHi.cpp
   ${EVTGEN_MODELS_DIR}/EvtFlatQ2.cpp
   ${EVTGEN_MODELS_DIR}/EvtMelikhovFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtBToKpipiCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVSCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtSSD_DirectCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVVHelAmp.cpp
   ${EVTGEN_MODELS_DIR}/EvtBaryonPCR.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgammaRootFinder.cpp
   ${EVTGEN_MODELS_DIR}/EvtD0mixDalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtLNuGamma.cpp
   ${EVTGEN_MODELS_DIR}/EvtbTosllBall.cpp
   ${EVTGEN_MODELS_DIR}/EvtTauHadnu.cpp
   ${EVTGEN_MODELS_DIR}/EvtBToPlnuBKFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtSSSCPpng.cpp
   ${EVTGEN_MODELS_DIR}/EvtBTo3piCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgammaFlatEnergy.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVVCP.cpp
   ${EVTGEN_MODELS_DIR}/EvtVVSPwave.cpp
   ${EVTGEN_MODELS_DIR}/EvtKstarnunu.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVVCPLH.cpp
   ${EVTGEN_MODELS_DIR}/EvtKstarstargamma.cpp
   ${EVTGEN_MODELS_DIR}/EvtVVPIPI_WEIGHTED.cpp
   ${EVTGEN_MODELS_DIR}/EvtDMix.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVS.cpp
   ${EVTGEN_MODELS_DIR}/EvtbTosllAliFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtDalitzTable.cpp
   ${EVTGEN_MODELS_DIR}/EvtBcToNPi.cpp
   ${EVTGEN_MODELS_DIR}/EvtbTosllAmp.cpp
   ${EVTGEN_MODELS_DIR}/EvtTaulnunu.cpp
   ${EVTGEN_MODELS_DIR}/EvtBcPsiNPi.cpp
   ${EVTGEN_MODELS_DIR}/EvtHQET2.cpp
   ${EVTGEN_MODELS_DIR}/EvtVubAC.cpp
   ${EVTGEN_MODELS_DIR}/EvtSLN.cpp
   ${EVTGEN_MODELS_DIR}/EvtBto2piCPiso.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgSimpsonIntegrator.cpp
   ${EVTGEN_MODELS_DIR}/EvtLambdaB2LambdaV.cpp
   ${EVTGEN_MODELS_DIR}/EvtBHadronic.cpp
   ${EVTGEN_MODELS_DIR}/EvtVSSMix.cpp
   ${EVTGEN_MODELS_DIR}/EvtCBTo3piP00.cpp
   ${EVTGEN_MODELS_DIR}/EvtTauVectornu.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgAbsFunction.cpp
   ${EVTGEN_MODELS_DIR}/EvtTSS.cpp
   ${EVTGEN_MODELS_DIR}/EvtPropSLPole.cpp
   ${EVTGEN_MODELS_DIR}/EvtSLPole.cpp
   ${EVTGEN_MODELS_DIR}/EvtLambdaP_BarGamma.cpp
   ${EVTGEN_MODELS_DIR}/EvtBToVlnuBall.cpp
   ${EVTGEN_MODELS_DIR}/EvtBaryonPCRFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtHelAmp.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVPHelAmp.cpp
   ${EVTGEN_MODELS_DIR}/EvtbTosllBallFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtSLPoleFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtPartWave.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVSNONCPEIGEN.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgammaAbsModel.cpp
   ${EVTGEN_MODELS_DIR}/EvtVPHOtoV.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsgammaFixedMass.cpp
   ${EVTGEN_MODELS_DIR}/EvtTauScalarnu.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVSCPiso.cpp
   ${EVTGEN_MODELS_DIR}/EvtHQET2FF.cpp
   ${EVTGEN_MODELS_DIR}/EvtVectorIsr.cpp
   ${EVTGEN_MODELS_DIR}/EvtBsquark.cpp
   ${EVTGEN_MODELS_DIR}/EvtItgAbsIntegrator.cpp
   ${EVTGEN_MODELS_DIR}/EvtBtoXsEtap.cpp
   ${EVTGEN_MODELS_DIR}/EvtKKLambdaC.cpp
   ${EVTGEN_MODELS_DIR}/EvtISGWFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtYmSToYnSpipiCLEO.cpp
   ${EVTGEN_MODELS_DIR}/EvtbTosllAli.cpp
   ${EVTGEN_MODELS_DIR}/EvtMultibody.cpp
   ${EVTGEN_MODELS_DIR}/EvtHQETFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtVubdGamma.cpp
   ${EVTGEN_MODELS_DIR}/EvtBTo3pi.F
   ${EVTGEN_MODELS_DIR}/EvtBToKpipi.F
   ${EVTGEN_MODELS_DIR}/EvtBTo3piP00.F
   ${EVTGEN_MODELS_DIR}/EvtLi2Spence.F
   ${EVTGEN_MODELS_DIR}/EvtBTo3piMPP.F
   ${EVTGEN_MODELS_DIR}/EvtBcVNpi.cpp
   ${EVTGEN_MODELS_DIR}/EvtBcVMuNu.cpp
   ${EVTGEN_MODELS_DIR}/EvtXPsiGamma.cpp
   ${EVTGEN_MODELS_DIR}/EvtD0gammaDalitz.cpp
   ${EVTGEN_MODELS_DIR}/EvtTVP.cpp
   ${EVTGEN_MODELS_DIR}/EvtSVP.cpp
   ${EVTGEN_MODELS_DIR}/EvtEta2MuMuGamma.cpp
   ${EVTGEN_MODELS_DIR}/EvtBCVFF.cpp
   ${EVTGEN_MODELS_DIR}/EvtWnPi.cpp
   ${EVTGEN_PRIVATE_DIR}/EvtPpbarJpsiPi0.cpp
   )

#   ${EVTGEN_MODELS_DIR}/EvtResonancePrototype.cpp


set(EVTGEN_EXT_SRCS
   ${EVTGEN_EXTERNAL_DIR}/EvtPythiaEngine.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtExternalGenList.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtPHOTOS.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtPythia6CommandConverter.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtPythia.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtPhotosEngine.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtTauolaEngine.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtExternalGenFactory.cpp
   ${EVTGEN_EXTERNAL_DIR}/EvtTauola.cpp
)

# Fix a problem with g++ on Mac OS X
#set_source_files_properties(${EVTGEN_MODELS_DIR}/EvtPythia.cc PROPERTIES COMPILE_FLAGS "-O0")

if(RULE_CHECKER_FOUND)
  CHECK_RULES("${EVTGEN_SRCS}" "${INCLUDE_DIRECTORIES}" EVTGEN_RULES)
endif(RULE_CHECKER_FOUND)

set(EVTGEN_HEADERS
   EvtGen.hh
)

set(EVTGEN_LINKDEF EvtGenLinkDef.h)
set(EVTGEN_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/EvtGenDict.cxx)

add_library(EvtGen SHARED ${EVTGEN_SRCS})

target_link_libraries(EvtGen ${ROOT_LIBRARIES} ${Fortran_UseLib} ${HEPMCLIBLIST})
set_target_properties(EvtGen PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES})
set_target_properties(EvtGen PROPERTIES LINK_FLAGS "-fPIC")

if(EVTGEN_EXTERNAL)
   add_library(EvtGenExternal SHARED ${EVTGEN_EXT_SRCS})
   target_link_libraries(EvtGenExternal ${ROOT_LIBRARIES} ${PHOTOSLIBLIST} ${TAUOLALIBLIST} ${HEPMCLIBLIST} ${PYTHIALIBLIST} ${Fortran_UseLib} EvtGen)
   set_target_properties(EvtGenExternal PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES})
   set_target_properties(EvtGenExternal PROPERTIES LINK_FLAGS "-fPIC")
endif(EVTGEN_EXTERNAL)

add_executable(evtgenlhc_test1 ${EVTGEN_EXAMPLE_DIR}/evtgenlhc_test1.cc)
set_target_properties(evtgenlhc_test1 PROPERTIES OUTPUT_NAME evtgen_lhc_test1)
set_target_properties(evtgenlhc_test1 PROPERTIES LINK_FLAGS "-fPIC")
target_link_libraries(evtgenlhc_test1 ${ROOT_LIBRARIES} EvtGen ${EVTGENEXTLIBLIST} ${PHOTOSLIBLIST} ${TAUOLALIBLIST} ${HEPMCLIBLIST} ${PYTHIALIBLIST} ${Fortran_UseLib} )


add_executable(evt_dalitz ${EVTGEN_EXAMPLE_DIR}/evt_dalitz.cc)
set_target_properties(evt_dalitz PROPERTIES OUTPUT_NAME evtgen_dalitz)
set_target_properties(evt_dalitz PROPERTIES LINK_FLAGS "-fPIC")
target_link_libraries(evt_dalitz ${ROOT_LIBRARIES} EvtGen ${EVTGENEXTLIBLIST} ${PHOTOSLIBLIST} ${TAUOLALIBLIST} ${HEPMCLIBLIST} ${PYTHIALIBLIST} ${Fortran_UseLib} )


add_executable(example1 ${EVTGEN_EXAMPLE_DIR}/example1.cc)
set_target_properties(example1 PROPERTIES OUTPUT_NAME evtgen_example1)
set_target_properties(example1 PROPERTIES LINK_FLAGS "-fPIC")
target_link_libraries(example1 ${ROOT_LIBRARIES} EvtGen ${EVTGENEXTLIBLIST} ${PHOTOSLIBLIST} ${TAUOLALIBLIST} ${HEPMCLIBLIST} ${PYTHIALIBLIST} ${Fortran_UseLib} )

################ install ###################
install(TARGETS EvtGen DESTINATION ${CMAKE_BINARY_DIR}/lib)

#!/bin/bash

function patch_src_file {
    SRC_FILE=$1
    if [[ -e $PNDROOT_DIR/$SRC_FILE ]]
    then
	if [[ -e $PATCH_DIR/$SRC_FILE ]]
	then
	    echo "####################################"
	    echo "####################################"
	    echo Patching $SRC_FILE:
	    nbkp=$(ls $PNDROOT_DIR/$SRC_FILE.org.* 2>/dev/null | wc -l)
	    diff -b $PNDROOT_DIR/$SRC_FILE $PATCH_DIR/$SRC_FILE
	    mv -vf $PNDROOT_DIR/$SRC_FILE $PNDROOT_DIR/$SRC_FILE.org.$nbkp
	    cp -vf $PATCH_DIR/$SRC_FILE $PNDROOT_DIR/$SRC_FILE
	else
	    echo Patch file $PATCH_DIR/$SRC_FILE not found. Cant patch
	fi
    else
	echo Source file $PNDROOT_DIR/$SRC_FILE not found. Cant patch
    fi
}

function add_src_file {
    SRC_FILE=$1
    SRC_FILE_DIR=$(dirname $SRC_FILE)
    # create parent directory if it doesn't exist
    [ -e $PNDROOT_DIR/$SRC_FILE_DIR ] || mkdir -p $PNDROOT_DIR/$SRC_FILE_DIR
    # for safety if file already exists
    [ -e $PNDROOT_DIR/$SRC_FILE ] && mv $PNDROOT_DIR/$SRC_FILE $PNDROOT_DIR/$SRC_FILE.org
    cp -vf $PATCH_DIR/$SRC_FILE $PNDROOT_DIR/$SRC_FILE
}

export EXT_VER=jul14p3
export PR_VER=26841
export PR_VER_DIR=trunk-$PR_VER

export SOFTDIR=/nfs1/panda/ermias/soft
export SIMPATH=$SOFTDIR/$EXT_VER/install
export FAIRROOTPATH=$SOFTDIR/$EXT_VER/FairRoot/install

echo EXT_VER: $EXT_VER
echo PR_VER: $PR_VER
echo PR_VER_DIR: $PR_VER_DIR
echo SOFTDIR: $SOFTDIR
echo SIMPATH: $SIMPATH
echo FAIRROOTPATH: $FAIRROOTPATH

cd $SOFTDIR/$EXT_VER
mkdir -p pandaroot
cd pandaroot
svn co https://subversion.gsi.de/fairroot/pandaroot/trunk $PR_VER_DIR -r $PR_VER
cd $PR_VER_DIR

# remove "FairRoot" directories because they have already been installed. This won't be necessary in future
# suggested way on wiki (svn propedit followed by svn update) doesn't work in batch environment
rm -vrf base generators geobase parbase cmake geane eventdisplay trackbase fairtools dbase MbsAPI

PNDROOT_DIR=$SOFTDIR/$EXT_VER/pandaroot/$PR_VER_DIR
PATCH_DIR=$SOFTDIR/patch-$PR_VER_DIR

#patch 1: "Brem" tag in PndAnalysis filtering
patch_src_file PndTools/AnalysisTools/PndAnalysis.cxx
patch_src_file PndTools/AnalysisTools/PndAnalysis.h

#patch 2: Comment out emcModule 3 matching condition
patch_src_file pid/PidCorr/PndPidEmcInfo.cxx

#patch 3: Modified Brem Corrector module for analysis ntuples
add_src_file pid/PidCorr/PndPidBremCorrectorNT.cxx
add_src_file pid/PidCorr/PndPidBremCorrectorNT.h
patch_src_file pid/CMakeLists.txt

#patch 4: Make RhoCandList compatible with stl containers
patch_src_file rho/RhoBase/RhoCandList.cxx
patch_src_file rho/RhoBase/RhoCandList.h

#patch 5: TDA generator in EvtPpbarJpsiPi0
add_src_file pgenerators/EvtGen/EvtGen/Private/EvtPpbarJpsiPi0.cpp
add_src_file pgenerators/EvtGen/EvtGen/Private/EvtGenModels/EvtPpbarJpsiPi0.hh
patch_src_file pgenerators/EvtGen/EvtGen/Private/EvtModelReg.cpp
patch_src_file pgenerators/EvtGen/EvtGen/CMakeLists.txt

#patch 6: Modified DPM generator (Decay pion on your own)
add_src_file pgenerators/PndDpmGeneratorMod.cxx
add_src_file pgenerators/PndDpmGeneratorMod.h
patch_src_file pgenerators/CMakeLists.txt

#patch 7: Standalone Event Generator
add_src_file pgenerators/EvtGen/EvtGenStandAlone/CMakeLists.txt
add_src_file pgenerators/EvtGen/EvtGenStandAlone/EvtGenStandAloneLinkDef.h
add_src_file pgenerators/EvtGen/EvtGenStandAlone/PndEvtGenStandAlone.cxx
add_src_file pgenerators/EvtGen/EvtGenStandAlone/PndEvtGenStandAlone.h
patch_src_file pgenerators/EvtGen/CMakeLists.txt

#patch 8: Private analysis modules
patch_src_file CMakeLists.txt
add_src_file src/CMakeLists.txt

#patch 8.1 Tda Analysis module
add_src_file src/AnaTda/AnaTda.cxx
add_src_file src/AnaTda/AnaTda.h
add_src_file src/AnaTda/AnaTdaLinkDef.h
add_src_file src/AnaTda/AnaTdav2.cxx
add_src_file src/AnaTda/AnaTdav2.h
add_src_file src/AnaTda/CMakeLists.txt

#patch 8.2 PhysListHadronic module
add_src_file src/PhysListHadrG4/HistFillerTask.cxx
add_src_file src/PhysListHadrG4/HistFillerLinkDef.h
add_src_file src/PhysListHadrG4/CMakeLists.txt
add_src_file src/PhysListHadrG4/HistFillerTask.h

#patch 8.3 PhysListHadronic module
add_src_file src/filler/filler.cxx
add_src_file src/filler/CMakeLists.txt
add_src_file src/filler/filler.h

mkdir -p build
cd build
cmake -DUSE_DIFFERENT_COMPILER=TRUE ../
make -j 4

cd $SOFTDIR/$EXT_VER
chmod -Rf g+w pandaroot

#!/bin/bash

if [[ $# != 2 ]]; then
    echo "usage: ./patch.pandaroot.sh EXT_VER PR_VER"
    echo "example1(locl): ./compile.pandaroot.sh ext-apr13 oct14"
    echo "example1(grid): ./compile.pandaroot.sh apr13 oct14"
    echo "example2(grid): ./compile.pandaroot.sh jul14p3 trunk-26841"
    exit
else
    EXT_VER=$1
    PR_VER=$2
    echo "Compiling EXT_VER=$EXT_VER and PR_VER=$PR_VER"
fi

SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
UTILS_SCRIPT=$SCRIPT_DIR/utils.sh
if [ ! -e $UTILS_SCRIPT ]; then
    echo "utils.sh not found in the same directory as compile script. quitting"
    exit
fi
. $UTILS_SCRIPT

#input sanitation
if ! check_path $PATCH_DIR; then exit; fi
if ! check_path $EXT_DIR; then exit; fi
if ! check_path $SIMPATH; then exit; fi
if [[ -e  $EXT_DIR/FairRoot ]]; then
    export FAIRROOTPATH=$EXT_DIR/FairRoot/install
    if ! check_path $FAIRROOTPATH; then exit; fi
fi
echo "Pathes check out. Proceeding to compilation"

cd $EXT_DIR

[ -e pandaroot ] || mkdir -p pandaroot
cd pandaroot
if [ ! -e $PR_VER ]; then
    svn co https://subversion.gsi.de/fairroot/pandaroot/release/$PR_VER $PR_VER
else
    echo "$PR_VER directory already exist. Refusing to checkout."
    echo "Will try to compile existing source"
fi

. $SCRIPT_DIR/patch.pandaroot.sh -n $EXT_VER $PR_VER

cd $PNDROOT_DIR
if [[ $PR_VER == "oct14" ]]; then
    rm -vrf base generators geobase parbase cmake geane eventdisplay trackbase fairtools dbase MbsAPI
fi

echo "SIMPATH= $SIMPATH"
echo "FAIRROOTPATH= $FAIRROOTPATH"

[ -e build ] && mv -v build build.bkp.$(ls -d build.bkp.* 2>/dev/null | wc -l)
mkdir -p build
cd build
cmake -DUSE_DIFFERENT_COMPILER=TRUE ../
make -j $NCORE

if [[ $HN != "ipnphen01" ]]; then
    cd $SOFT_DIR/$PR_VER
    chmod -Rf g+w pandaroot
fi

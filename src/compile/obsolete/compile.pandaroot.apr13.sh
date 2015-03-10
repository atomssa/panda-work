#!/bin/bash

export HN=$(hostname)
export VER=apr13

if [[ $HN == "ipnphen01" ]]; then
    export SOFTDIR=/vol0/panda/svn
    export EXTDIR=$SOFTDIR/ext-apr13
    export SIMPATH=$SOFTDIR/ext-apr13
    export PATCH_DIR=/vol0/panda/work/src/patch-$VER
    export NCORE=16
elif [[ $HN == "ipngrid02" ]]; then
    export SOFTDIR=/nfs1/panda/ermias/soft
    export EXTDIR=$SOFTDIR/apr13
    export SIMPATH=$EXTDIR/install
    export PATCH_DIR=$SOFTDIR/patch-$VER
    export NCORE=4
else
    echo "Unrecongnized host name: " $HN
    exit
fi

cd $EXTDIR

[ -e pandaroot ] || mkdir -p pandaroot
cd pandaroot
svn co https://subversion.gsi.de/fairroot/pandaroot/release/$VER $VER
cd $VER

#PNDROOT_DIR=$EXTDIR/pandaroot/$VER

mkdir -p build
cd build
cmake ../
make -j $NCORE

if [[ $HN == "ipngrid02" ]]; then
    cd $SOFTDIR/$VER
    chmod -Rf g+w pandaroot
fi

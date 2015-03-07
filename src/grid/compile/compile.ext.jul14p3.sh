#!/bin/bash

SOFTDIR=/nfs1/panda/ermias/soft
cd $SOFTDIR
echo PWD: $PWD
ls -altr

rm -vrf jul14p3

# FairSoft
echo installing jul14p3 externals release version
git clone https://github.com/FairRootGroup/FairSoft jul14p3
cd $SOFTDIR/jul14p3
git checkout -b jul14p3 jul14p3
cp -vf $SOFTDIR/configure-jul14p3.sh ./configure.sh
. ./configure.sh grid

cd $SOFTDIR/jul14p3/install
export SIMPATH=$PWD

# FairRoot
cd $SOFTDIR/jul14p3
git clone https://github.com/FairRootGroup/FairRoot.git
_FAIRROOTDIR=$SOFTDIR/jul14p3/FairRoot
cd $_FAIRROOTDIR
git checkout -b v-14.11 v-14.11
mkdir $_FAIRROOTDIR/build
cd $_FAIRROOTDIR/build
cmake -DCMAKE_INSTALL_PREFIX=$_FAIRROOTDIR/install -DUSE_DIFFERENT_COMPILER=TRUE $_FAIRROOTDIR
make -j 4
make install

# regain permission
cd $SOFTDIR
chmod -Rf g+w jul14p3

#!/bin/bash

HN=$(hostname)

if [[ $HN == "ipnphen01" ]];then
    SOFT_DIR=/vol0/panda/svn
else
    SOFT_DIR=/nfs1/panda/ermias/soft
fi
cd $SOFT_DIR
echo PWD: $PWD
ls -altr

[ -e jul14p3 ] && mv -v build jul14p3.bkp.$(ls -d jul14p3.bkp.* 2>/dev/null | wc -l)

# FairSoft
echo installing jul14p3 externals release version
git clone https://github.com/FairRootGroup/FairSoft jul14p3
cd $SOFT_DIR/jul14p3
git checkout -b jul14p3 jul14p3

if [[ $HN == "ipnphen01" ]]; then
    . ./configure.sh <<EOF
1
1
2
1
$SOFT_DIR/jul14p3/install
2
EOF
else
    sed -i.org.0 s/geant4_download_install_data_automatic=no/geant4_download_install_data_automatic=yes/ configure.sh
    sed -i.org.1 s/geant4_install_data_from_dir=yes/geant4_install_data_from_dir=no/ configure.sh
    sed -i.org.2 s/'SIMPATH_INSTALL=$PWD\/installation'/'export SIMPATH_INSTALL=$PWD\/install'/ configure.sh
    . ./configure.sh grid
fi
export SIMPATH=$SOFT_DIR/jul14p3/install

# FairRoot
cd $SOFT_DIR/jul14p3
git clone https://github.com/FairRootGroup/FairRoot.git
FAIRROOT_DIR=$SOFT_DIR/jul14p3/FairRoot
cd $FAIRROOT_DIR
git checkout -b v-14.11 v-14.11
mkdir $FAIRROOT_DIR/build
cd $FAIRROOT_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=$FAIRROOT_DIR/install -DUSE_DIFFERENT_COMPILER=TRUE $FAIRROOT_DIR
make -j 4
make install

if [[ $HN != "ipnphen01" ]]; then
   # regain permission
   cd $SOFT_DIR
   chmod -Rf g+w jul14p3
fi

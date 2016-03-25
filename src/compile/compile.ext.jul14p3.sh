#!/bin/bash

HN=$(hostname)

if [[ $HN == "ipnphen01" ]];then
    SOFT_DIR=/vol0/panda/svn
elif [[ $HN == "rasalula" ]];then
    SOFT_DIR=/Users/tujuba/panda/svn
elif [[ $HN == "ipngrid01.in2p3.fr" ]]; then
    SOFT_DIR=/nfs1/panda/ermias/soft_ipngrid01
else
    SOFT_DIR=/nfs1/panda/ermias/soft
fi
cd $SOFT_DIR
echo PWD: $PWD
ls -altr

[ -e jul14p3 ] && mv -v jul14p3 jul14p3.bkp.$(ls -d jul14p3.bkp.* 2>/dev/null | wc -l)

# FairSoft
echo installing jul14p3 externals release version
git clone https://github.com/FairRootGroup/FairSoft jul14p3
cd $SOFT_DIR/jul14p3
git checkout -b jul14p3 jul14p3

# macbook can't compile early versions of root
if [[ $HN == "rasalula" ]]; then
    cd $SOFT_DIR/jul14p3/scripts
    sed -i.bkp.0 s/'ROOTVERSION=87c1d2a0e691781d2ac93a4f05afee1436feaf4a'/'ROOTVERSION=v5-34-26'/ package_versions.sh
    cd $SOFT_DIR/jul14p3
fi

if [[ $HN == "ipnphen01" || $HN == "rasalula" || $HN == "ipngrid01.in2p3.fr" ]]; then
    . ./configure.sh <<EOF
1
1
2
1
$SOFT_DIR/jul14p3/install
2
EOF
else
    sed -i.bkp.0 s/geant4_download_install_data_automatic=no/geant4_download_install_data_automatic=yes/ configure.sh
    sed -i.bkp.1 s/geant4_install_data_from_dir=yes/geant4_install_data_from_dir=no/ configure.sh
    sed -i.bkp.2 s/'SIMPATH_INSTALL=$PWD\/installation'/'export SIMPATH_INSTALL=$PWD\/install'/ configure.sh
    . ./configure.sh grid
fi
export SIMPATH=$SOFT_DIR/jul14p3/install

# FairRoot
cd $SOFT_DIR/jul14p3
git clone https://github.com/FairRootGroup/FairRoot.git
FAIRROOT_DIR=$SOFT_DIR/jul14p3/FairRoot
cd $FAIRROOT_DIR
git checkout -b v-14.11 v-14.11

# CUDA and PROTObuf not compileable on OSX, and not required
if [[ $HN == "rasalula" ]]; then
    sed -i.bkp.0 s/'find_package(CUDA)'/'#find_package(CUDA)'/ CMakeLists.txt
    sed -i.bkp.0 s/'find_package(Protobuf)'/'#find_package(Protobuf)'/ CMakeLists.txt
fi

mkdir $FAIRROOT_DIR/build
cd $FAIRROOT_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=$FAIRROOT_DIR/install -DUSE_DIFFERENT_COMPILER=TRUE $FAIRROOT_DIR
make -j 4
make install

if [[ $HN != "ipnphen01" && $HN != "rasalula" && $HN != "ipngrid01.in2p3.fr" ]]; then
    # regain permission
    cd $SOFT_DIR
    chmod -Rf g+w jul14p3
fi

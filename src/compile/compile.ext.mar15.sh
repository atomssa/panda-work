#!/bin/bash

HN=$(hostname)

if [[ $HN == "ipnphen01" ]];then
    SOFT_DIR=/vol0/panda/svn
elif [[ $HN == "rasalula" ]];then
    SOFT_DIR=/Users/tujuba/panda/svn
else
    SOFT_DIR=/nfs1/panda/ermias/soft
fi
cd $SOFT_DIR
echo PWD: $PWD
ls -altr

VFS=mar15
VFR=v-15.03

[ -e ${VFS} ] && mv -v build ${VFS}.bkp.$(ls -d ${VFS}.bkp.* 2>/dev/null | wc -l)

# FairSoft
echo installing ${VFS} externals release version
git clone https://github.com/FairRootGroup/FairSoft ${VFS}
cd $SOFT_DIR/${VFS}
git checkout -b ${VFS} ${VFS}

# macbook can't compile early versions of root
#if [[ $HN == "rasalula" ]]; then
#    cd $SOFT_DIR/${VFS}/scripts
#    sed -i.bkp.0 s/'ROOTVERSION=87c1d2a0e691781d2ac93a4f05afee1436feaf4a'/'ROOTVERSION=v5-34-26'/ package_versions.sh
#    cd $SOFT_DIR/${VFS}
#fi

if [[ $HN == "rasalula" ]]; then
    . ./configure.sh <<EOF
5
1
1
2
1
$SOFT_DIR/${VFS}/install
2
EOF
elif [[ $HN == "ipnphen01" ]]; then
    . ./configure.sh <<EOF
1
1
1
2
1
$SOFT_DIR/${VFS}/install
2
EOF
else
    sed -i.bkp.0 s/geant4_download_install_data_automatic=no/geant4_download_install_data_automatic=yes/ configure.sh
    sed -i.bkp.1 s/geant4_install_data_from_dir=yes/geant4_install_data_from_dir=no/ configure.sh
    sed -i.bkp.2 s/'SIMPATH_INSTALL=$PWD\/installation'/'export SIMPATH_INSTALL=$PWD\/install'/ configure.sh
    . ./configure.sh grid
fi
export SIMPATH=$SOFT_DIR/${VFS}/install

# FairRoot
cd $SOFT_DIR/${VFS}
git clone https://github.com/FairRootGroup/FairRoot.git
FAIRROOT_DIR=$SOFT_DIR/${VFS}/FairRoot
cd $FAIRROOT_DIR
git checkout -b $VFR $VFR

## CUDA and PROTObuf not compileable on OSX, and not required
#if [[ $HN == "rasalula" ]]; then
#    sed -i.bkp.0 s/'find_package(CUDA)'/'#find_package(CUDA)'/ CMakeLists.txt
#    sed -i.bkp.0 s/'find_package(Protobuf)'/'#find_package(Protobuf)'/ CMakeLists.txt
#fi

mkdir $FAIRROOT_DIR/build
cd $FAIRROOT_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=$FAIRROOT_DIR/install -DUSE_DIFFERENT_COMPILER=TRUE $FAIRROOT_DIR
make -j 4
make install

if [[ $HN != "ipnphen01" && $HN != "rasalula" ]]; then
    # regain permission
    cd $SOFT_DIR
    chmod -Rf g+w ${VFS}
fi

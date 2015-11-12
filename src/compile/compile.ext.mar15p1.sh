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

VFS=mar15p1
VFR=v-15.03

[ -e ${VFS} ] && mv -v build ${VFS}.bkp.$(ls -d ${VFS}.bkp.* 2>/dev/null | wc -l)

# FairSoft
echo installing ${VFS} externals release version
git clone https://github.com/FairRootGroup/FairSoft ${VFS}
cd $SOFT_DIR/${VFS}
git checkout -b ${VFS} ${VFS}

if [[ $HN == "rasalula" ]]; then
    # on new mac-osx (El Capitan and more) openssl is not shipped with OS
    # no option but to use the homebrew (still bottled) version
    export OPENSSL_ROOT_DIR=/usr/local/opt/openssl
    ./configure.sh <<EOF
5
1
1
2
1
$SOFT_DIR/${VFS}/install
2
EOF
elif [[ $HN == "ipnphen01" || $HN == "ipngrid01.in2p3.fr" ]]; then
    ./configure.sh <<EOF
1
1
1
2
1
$SOFT_DIR/${VFS}/install
2
EOF
else
    ./configure.sh <<EOF
1
1
1
2
1
$SOFT_DIR/${VFS}/install
2
EOF
fi
export SIMPATH=$SOFT_DIR/${VFS}/install

# FairRoot
cd $SOFT_DIR/${VFS}
git clone https://github.com/FairRootGroup/FairRoot.git
FAIRROOT_DIR=$SOFT_DIR/${VFS}/FairRoot
cd $FAIRROOT_DIR
git checkout -b $VFR $VFR

mkdir $FAIRROOT_DIR/build
cd $FAIRROOT_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=$FAIRROOT_DIR/install -DUSE_DIFFERENT_COMPILER=TRUE $FAIRROOT_DIR
make -j 4
make install

if [[ $HN != "ipnphen01" && $HN != "rasalula" && $HN != "ipngrid01.in2p3.fr" ]]; then
    # regain permission
    cd $SOFT_DIR
    chmod -Rf g+w ${VFS}
fi

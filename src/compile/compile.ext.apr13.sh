#!/bin/bash

HN=$(hostname)

if [[ $HN == "ipnphen01" ]]; then
    SOFT_DIR=/vol0/panda/svn
else
    SOFT_DIR=/nfs1/panda/ermias/soft
fi

cd $SOFT_DIR
echo PWD: $PWD
ls -altr

[ -e apr13 ] && mv -v build apr13.bkp.$(ls -d apr13.bkp.* 2>/dev/null | wc -l)

echo installing apr13 externals release version
svn co https://subversion.gsi.de/fairroot/fairsoft/release/apr13

echo root is missing for apr13 release download it manually...
cd $SOFT_DIR/apr13/tools
svn co https://root.cern.ch/svn/root/tags/v5-34-05 root

cd $SOFT_DIR/apr13
echo PWD: $PWD
ls -altr

if [[ $HN == "ipnphen01" ]]; then
    . ./configure.sh <<EOF
1
1
2
$SOFT_DIR/apr13/install
2
EOF
else
# for grid install, we want the geant4 data archives to be
# downloaded and installed automatically.
    sed -i.org s/geant4_download_install_data_automatic=no/geant4_download_install_data_automatic=yes/ configure.sh
    . ./configure.sh grid
fi

if [[ $HN != "ipnphen01" ]]; then
    cd /nfs1/panda/ermias/soft
    chmod -Rf g+w *
fi

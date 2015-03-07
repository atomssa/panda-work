#!/bin/bash

cd /nfs1/panda/ermias/soft.fail1
chmod -Rf g+w *

cd /nfs1/panda/ermias/soft.fail2
chmod -Rf g+w *

cd /nfs1/panda/ermias/soft

echo PWD: $PWD
ls -altr

echo installing apr13 externals release version
svn co https://subversion.gsi.de/fairroot/fairsoft/release/apr13

echo root is missing download it manually...
cd /nfs1/panda/ermias/soft/apr13/tools
svn co https://root.cern.ch/svn/root/tags/v5-34-05 root

cd /nfs1/panda/ermias/soft/apr13
echo PWD: $PWD
ls -altr

mv configure.sh configure.sh.org
cp -f ../configure.sh .

. ./configure.sh grid

cd /nfs1/panda/ermias/soft
chmod -Rf g+w *

#echo installing oct14 pandaroot release version
#export SIMPATH=$PWD
#svn co https://subversion.gsi.de/fairroot/pandaroot/oct14 pandaroot_oct14
#cd pandaroot_oct14

## Patch oct14
# git checkout ... somepatchfile.patch
# patch -i somepatchfile.patch

## launch build
#mkdir build
#cd build
#cmake ../
#make -j 4

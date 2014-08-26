#!/bin/bash

# This macro is needed to start the Root macros used for automatic testing
# from inside CMake using the add_test functionality. Since the tests 
# starts with a fresh environment on has to set first the correct environment
# needed to run FairRoot.
# Also parameters defined in add_test will be converted in the correct format
# to be passed to root.

export SIMPATH=/Users/tujuba/panda/ext

# Setup the needed environment
prepend LD_LIBRARY_PATH '/Users/tujuba/panda/ext/lib'
prepend LD_LIBRARY_PATH '/Users/tujuba/panda/ext/lib/root'
prepend LD_LIBRARY_PATH '/Users/tujuba/panda/root/build/lib'
export LD_LIBRARY_PATH=/Users/tujuba/panda/root/build/lib:$LD_LIBRARY_PATH #No idea why this has to be put in twice

prepend DYLD_LIBRARY_PATH '/Users/tujuba/panda/ext/lib'
prepend DYLD_LIBRARY_PATH '/Users/tujuba/panda/ext/lib/root'
prepend DYLD_LIBRARY_PATH '/Users/tujuba/panda/root/build/lib'
export DYLD_LIBRARY_PATH=/Users/tujuba/panda/root/build/lib:$DYLD_LIBRARY_PATH #No idea why this has to be put in twice

prepend PATH '/Users/tujuba/panda/ext/bin'

export ROOTSYS=/Users/tujuba/panda/ext/tools/root
export ROOTEXE=/Users/tujuba/panda/ext/bin/root.exe
export VMCWORKDIR=/Users/tujuba/panda/root
export GEANT4VMC_MACRO_DIR=/Users/tujuba/panda/ext/transport/macro
export USE_VGM=1
export BOOST=
export G3SYS="/Users/tujuba/panda/ext/share/geant3"

if [ -e /Users/tujuba/panda/ext/bin/geant4.sh ]; then
    . /Users/tujuba/panda/ext/bin/geant4.sh
fi

if [ -e /Users/tujuba/panda/ext/bin/env.sh ]; then
  . /Users/tujuba/panda/ext/bin/env.sh
fi

#!/bin/bash

# This macro is needed to start the Root macros used for automatic testing
# from inside CMake using the add_test functionality. Since the tests 
# starts with a fresh environment on has to set first the correct environment
# needed to run FairRoot.
# Also parameters defined in add_test will be converted in the correct format
# to be passed to root.

export SIMPATH=/Users/tujuba/panda/apr13

# Setup the needed environment
prepend LD_LIBRARY_PATH '/Users/tujuba/panda/apr13/lib'
prepend LD_LIBRARY_PATH '/Users/tujuba/panda/apr13/lib/root'
prepend LD_LIBRARY_PATH '/Users/tujuba/panda/pandaroot/build/lib'
export LD_LIBRARY_PATH=/Users/tujuba/panda/pandaroot/build/lib:$LD_LIBRARY_PATH #No idea why this has to be put in twice

prepend DYLD_LIBRARY_PATH '/Users/tujuba/panda/apr13/lib'
prepend DYLD_LIBRARY_PATH '/Users/tujuba/panda/apr13/lib/root'
prepend DYLD_LIBRARY_PATH '/Users/tujuba/panda/pandaroot/build/lib'
export DYLD_LIBRARY_PATH=/Users/tujuba/panda/pandaroot/build/lib:$DYLD_LIBRARY_PATH #No idea why this has to be put in twice

prepend PATH '/Users/tujuba/panda/apr13/bin'

export ROOTSYS=/Users/tujuba/panda/apr13/tools/root
export ROOTEXE=/Users/tujuba/panda/apr13/bin/root.exe
export VMCWORKDIR=/Users/tujuba/panda/pandaroot
export GEANT4VMC_MACRO_DIR=/Users/tujuba/panda/apr13/transport/macro
export USE_VGM=1
export BOOST=
export G3SYS="/Users/tujuba/panda/apr13/share/geant3"

if [ -e /Users/tujuba/panda/apr13/bin/geant4.sh ]; then
    . /Users/tujuba/panda/apr13/bin/geant4.sh
fi

if [ -e /Users/tujuba/panda/apr13/bin/env.sh ]; then
  . /Users/tujuba/panda/apr13/bin/env.sh
fi

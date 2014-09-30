#!/bin/bash

# usage : fill x y z (x=1 => Sig, y=0 => no eff, y=1 => eff, z=0 => unfilt, z=1 => filt)
./fill 1 0 0 # unfiltered

./fill 1 0 1 # filtered with no eff

./fill 1 1 1 # filtered with eff

#!/bin/bash

if [[ $1 == 0 ]]
then
    charge=esim_trunk
elif [[ $1 == 1 ]]
then
     charge=psim_trunk
else
    echo "Wrong arg. 0=>esim_trunk 1=>psim_trunk"
    exit 1
fi

. /vol0/panda/svn/fairsoft_jul14p3/pandaroot-t26841/build/config.sh

for i in `seq 0 49`; do
    echo $i
    indir=/vol0/panda/work/jpsi_pi0/grid.out/$charge/runall.$i.out/
    root -b <<EOF
.x /vol0/panda/work/jpsi_pi0/debug/mcfix/test.C("$indir")
.q
EOF
done

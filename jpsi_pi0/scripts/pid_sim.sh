#!/bin/bash

pnd scrut14

for i in 1 2 3 4 5;
do
    cd /vol0/panda/work/jpsi_pi0/.batch/tmp_${1}
    root -b <<EOF
.x pid_complete.C
.q
EOF
    echo "done with id " $i
done

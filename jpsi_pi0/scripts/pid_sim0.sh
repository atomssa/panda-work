#!/bin/bash

pnd scrut14

for i in 6 7 8 9 10 11 12 13 14 15 16 17 18;
do
    cd /vol0/panda/work/jpsi_pi0/.batch/tmp_${i}
    root -b <<EOF
.x pid_complete.C
.q
EOF
    mv pid_complete.root /vol0/panda/work/jpsi_pi0/output/pid/pbar_p_pip_pim_pi0_${i}.root
    echo "done with id " $i
done

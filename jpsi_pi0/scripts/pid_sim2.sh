#!/bin/bash

pnd scrut14

for i in 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130;
do
    cd /vol0/panda/work/jpsi_pi0/.batch/tmp_${i}
    root -b <<EOF
.x pid_complete.C
.q
EOF
    mv pid_complete.root /vol0/panda/work/jpsi_pi0/output/pid/pbar_p_jpsi_pi0_${i}.root
    echo "done with id " $i
done

#!/bin/bash

batch=$1
echo $batch

root -b <<EOF >& output/log${batch}.log
.x sim_complete_prod.C($batch)
.q
EOF


#

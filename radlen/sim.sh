#!/bin/bash

batch=$1
echo $batch

root -b <<EOF >& output/log${batch}_g4.log
.x sim_complete_prod.C($batch, "TGeant4")
.q
EOF


root -b <<EOF >& output/log${batch}_g3.log
.x sim_complete_prod.C($batch, "TGeant3")
.q
EOF


#

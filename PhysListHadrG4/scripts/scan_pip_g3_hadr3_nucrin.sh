#!/bin/bash

PHYS_LIST=G3_HADR3_NUCRIN

NEVT=50000

source /vol0/panda/Example_Env_Setup.sh

echo VMCWORKDIR: $VMCWORKDIR
echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH
echo PhysList: $PHYS_LIST

CONFIG_DIR=/vol0/panda/panda-apr13/pandaroot/macro/PhysListHadrG4/gconfig/$PHYS_LIST
echo CONFIG_DIR: $CONFIG_DIR


MOM=1.0
CHARGE=+1

SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`

echo Processing $PHYS_LIST $MOM GeV EMCal Only for $SPECIES

OUT_DAT=output/data_${PHYS_LIST}_MOM_${MOM}_${SPECIES}.root
OUT_PAR=output/simpar_${PHYS_LIST}_MOM_${MOM}_${SPECIES}.root
LOG=logs/${PHYS_LIST}_MOM_${MOM}_${SPECIES}.log
echo OUT_DAT: $OUT_DAT
echo OUT_PAR: $OUT_PAR
echo LOG: $LOG

root -b <<EOF >& $LOG 
.x emc_complete.C($NEVT,$MOM,$CHARGE,"$PHYS_LIST",0,"$OUT_DAT","$OUT_PAR")
.q
EOF
echo Done processing $PHYS_LIST $MOM GeV EMCal Only for $SPECIES



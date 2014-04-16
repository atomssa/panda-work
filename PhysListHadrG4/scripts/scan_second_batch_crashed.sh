#!/bin/bash

PHYS_LIST=$1

NEVT=50000

source /vol0/panda/Example_Env_Setup.sh

echo VMCWORKDIR: $VMCWORKDIR
echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH
echo PhysList: $PHYS_LIST

CONFIG_DIR=/vol0/panda/panda-apr13/pandaroot/macro/PhysListHadrG4/gconfig/$PHYS_LIST
echo CONFIG_DIR: $CONFIG_DIR

for MOM in 3.0 5.0;
do
    for CHARGE in +1 -1;
    do
	SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`

	echo Processing $PHYS_LIST $MOM GeV EMCal Only for $SPECIES

	OUT_DAT=output/data_${PHYS_LIST}_MOM_${MOM}_${SPECIES}.root
	OUT_PAR=output/simpar_${PHYS_LIST}_MOM_${MOM}_${SPECIES}.root
	LOG=logs/${PHYS_LIST}_MOM_${MOM}_${SPECIES}.log
	echo OUT_DAT: $OUT_DAT
	echo OUT_PAR: $OUT_PAR
	echo LOG: $LOG

	if [[ "$MOM" == "3.0" ]];
	then
	    if [[ "$CHARGE" == "+1" ]];
	    then 
		echo PHYS_LIST $PHYS_LIST MOM $MOM CHARGE $CHARGE already processed. not redoing
		continue
	    fi
	    if [[ "$PHYS_LIST" != "G3_HADR4_FLUKA" ]] && [[ "$PHYS_LIST" != "G3_HADR5_MICAP" ]] && [[ "$PHYS_LIST" != "G3_HADR6_GCALOR" ]] ;
	    then
		echo PHYS_LIST $PHYS_LIST MOM $MOM CHARGE $CHARGE already processed. not redoing
		continue
	    fi
	fi
	echo PHYS_LIST $PHYS_LIST MOM $MOM CHARGE $CHARGE crashed. redoing

	root -b <<EOF >& $LOG 
.x emc_complete.C($NEVT,$MOM,$CHARGE,"$PHYS_LIST",0,"$OUT_DAT","$OUT_PAR")
.q
EOF
	echo Done processing $PHYS_LIST $MOM GeV EMCal Only for $SPECIES
    done
done

#MOM=1.5
#for CHARGE in +1 -1;
#do
#    SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`
#
#    echo Processing $PHYS_LIST $MOM GeV Full PANDA for $SPECIES
#
#    OUT_DAT=output/data_FULLPANDA_${PHYS_LIST}_MOM_${MOM}_${SPECIES}.root
#    OUT_PAR=output/simpar_FULLPANDA_${PHYS_LIST}_MOM_${MOM}_${SPECIES}.root
#    LOG=logs/FULLPANDA_${PHYS_LIST}_MOM_${MOM}_${SPECIES}.log
#    echo OUT_DAT: $OUT_DAT
#    echo OUT_PAR: $OUT_PAR
#    echo LOG: $LOG
#
#    root -b <<EOF >& $LOG 
#.x emc_complete.C($NEVT,$MOM,$CHARGE,"$PHYS_LIST",1,"$OUT_DAT","$OUT_PAR")
#.q
#EOF
#    echo Done processing $PHYS_LIST $MOM GeV Full PANDA for $SPECIES
#done



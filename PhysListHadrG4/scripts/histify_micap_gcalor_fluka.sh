#!/bin/bash

echo VMCWORKDIR=$VMCWORKDIR
echo ROOTSYS=$ROOTSYS

function run_root_macro {

    FILETAG=$1
    LOGFILE=logs/hist_${FILETAG}.log

    root -b <<EOF >& $LOGFILE 
.x hist_filler.C("$FILETAG")
.q
EOF

}

for BATCH in G3_HADR6_GCALOR G3_HADR5_MICAP G3_HADR4_FLUKA G3_HADR1_GEANH
do
    for MOM in 0.5 0.8 1.0 1.5 3.0 5.0;
    do
	for CHARGE in +1 -1;
	do
	    
	    SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`
	    
	    FILETAG=${BATCH}_MOM_${MOM}_${SPECIES}
	    
	    echo FILETAG=$FILETAG

	    if [ -e output/data_${FILETAG}.root ] 
	    then
		echo "Running ROOT for " $FILETAG
		run_root_macro $FILETAG
	    fi

	    #if [ $MOM == "1.0" ];
	    #then
	    #	FILETAG=FULLPANDA_${FILETAG}
	    #	
	    #	if [ -e output/data_${FILETAG}.root ]
	    #	then 
	    #	    echo "Running ROOT for " $FILETAG
	    #	    run_root_macro $FILETAG
	    #	fi
	    #fi

	done
    done
done


#BATCH=QGSP_BERT_HP
#MOM=1.0
#for CHARGE in +1 -1;
#do
#    
#    SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`
#    
#    FILETAG=${BATCH}_MOM_${MOM}_${SPECIES}
#    
#    echo FILETAG=$FILETAG
#
#    if [ -e output/data_${FILETAG}.root ] 
#    then
#	echo "Running ROOT for " $FILETAG
#	run_root_macro $FILETAG
#    fi
#
#    if [ $MOM == "1.0" ];
#    then
#	FILETAG=FULLPANDA_${FILETAG}
#	
#	if [ -e output/data_${FILETAG}.root ]
#	then 
#	    echo "Running ROOT for " $FILETAG
#	    run_root_macro $FILETAG
#	fi
#    fi
#
#done





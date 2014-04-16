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

MOM=1.0
for BATCH in QGSP_BIC_EMV QGS_BIC_EMV FTFP_BERT_EMV FTF_BIC_EMV LHEP QGSP_INCLXX_EMV QGSP_BERT_HP FTFP_BERT_TRV_EMV
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
	#    FILETAG=FULLPANDA_${FILETAG}
	#    
	#    if [ -e output/data_${FILETAG}.root ]
	#    then 
	#	echo "Running ROOT for " $FILETAG
	#	run_root_macro $FILETAG
	#    fi
	#fi

    done
done


#BATCH=$1
#for MOM in 0.5 0.8 1.0 1.5 3.0 5.0;
#do
#    for CHARGE in +1 -1;
#    do
#	
#	SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`
#	
#	FILETAG=${BATCH}_MOM_${MOM}_${SPECIES}
#	
#	echo FILETAG=$FILETAG
#
#	if [ -e output/data_${FILETAG}.root ] 
#	then
#	    echo "Running ROOT for " $FILETAG
#	    run_root_macro $FILETAG
#	fi
#
#	#if [ $MOM == "1.0" ];
#	#then
#	#    FILETAG=FULLPANDA_${FILETAG}
#	#    
#	#    if [ -e output/data_${FILETAG}.root ]
#	#    then 
#	#	echo "Running ROOT for " $FILETAG
#	#	run_root_macro $FILETAG
#	#    fi
#	#fi
#
#    done
#done

#for BATCH in G3_HADR3_NUCRIN QGSP_BERT_EMV QGSP_BIC_EMV QGSP_BERT_EMV_OPTICAL QGSP_BERT_OPTICAL QGSC_BERT_EMV QGSP_BERT_HP
#do
#    for MOM in 0.5 0.8 1.0 1.5 3.0 5.0;
#    do
#	for CHARGE in +1 -1;
#	do
#	    
#	    SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`
#	    
#	    FILETAG=${BATCH}_MOM_${MOM}_${SPECIES}
#	    
#	    echo FILETAG=$FILETAG
#
#	    if [ -e output/data_${FILETAG}.root ] 
#	    then
#		echo "Running ROOT for " $FILETAG
#		#run_root_macro $FILETAG
#	    fi
#
#	    if [ $MOM == "1.0" ];
#	    then
#		FILETAG=FULLPANDA_${FILETAG}
#		
#		if [ -e output/data_${FILETAG}.root ]
#		then 
#		    echo "Running ROOT for " $FILETAG
#		    run_root_macro $FILETAG
#		fi
#	    fi
#
#	done
#    done
#done


#for BATCH in QGSP_BERT_HP
#do
#    for MOM in 0.5 0.8 1.0 1.5 3.0 5.0;
#    do
#	for CHARGE in +1 -1;
#	do
#	    
#	    SPECIES=`if (( $CHARGE > 0 )); then echo pip; else echo pim; fi`
#	    
#	    FILETAG=${BATCH}_MOM_${MOM}_${SPECIES}
#	    
#	    echo FILETAG=$FILETAG
#
#	    if [ -e output/data_${FILETAG}.root ] 
#	    then
#		echo "Running ROOT for " $FILETAG
#		run_root_macro $FILETAG
#	    fi
#
#	    #if [ $MOM == "1.0" ];
#	    #then
#	    #	FILETAG=FULLPANDA_${FILETAG}
#	    #	
#	    #	if [ -e output/data_${FILETAG}.root ]
#	    #	then 
#	    #	    echo "Running ROOT for " $FILETAG
#	    #	    run_root_macro $FILETAG
#	    #	fi
#	    #fi
#
#	done
#    done
#done
#

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





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

#0
#data_G3_HADR1_GEANH :        12
#data_G3_HADR3_NUCRIN :        12   +    data_FULLPANDA_G3_HADR3_NUCRIN :         2
#data_G3_HADR4_FLUKA :        12
#data_G3_HADR5_MICAP :        12
#
#1
#data_G3_HADR6_GCALOR :        12
#data_FTFP_BERT_EMV :        12
#data_FTFP_BERT_TRV_EMV :        12
#data_QGSC_BERT_EMV :        12   + data_FULLPANDA_QGSC_BERT_EMV :         2
#
#2
#data_LHEP :        12
#data_FTF_BIC_EMV :        12
#data_QGSP_BERT_EMV :        12   + data_FULLPANDA_QGSP_BERT_EMV :         2
#data_QGSP_BERT_HP :        12    + data_FULLPANDA_QGSP_BERT_HP :         2
#
#3
#data_QGSP_BIC_EMV :        12    + data_FULLPANDA_QGSP_BIC_EMV :         2
#data_QGSP_EMV :        12
#data_QGSP_FTFP_BERT_EMV :        12
#data_QGSP_INCLXX_EMV :        12
#data_QGS_BIC_EMV :        12

for BATCH in G3_HADR6_GCALOR FTFP_BERT_EMV FTFP_BERT_TRV_EMV QGSC_BERT_EMV
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

	    if [ $MOM == "1.0" ];
	    then
		FILETAG=FULLPANDA_${FILETAG}
		
		if [ -e output/data_${FILETAG}.root ]
		then 
		    echo "Running ROOT for " $FILETAG
		    run_root_macro $FILETAG
		fi
	    fi

	done
    done
done

#

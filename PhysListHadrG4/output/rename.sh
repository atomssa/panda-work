#!/bin/bash

BATCH=$1

#mv -vi simpar_${BATCH}_MOM_1.0_pim.root simpar_FULLPANDA_${BATCH}_MOM_1.0_pim.root
#mv -vi data_${BATCH}_MOM_1.0_pim.root data_FULLPANDA_${BATCH}_MOM_1.0_pim.root
#mv -vi FairRunInfo_data_${BATCH}_MOM_1.0_pim.root FairRunInfo_data_FULLPANDA_${BATCH}_MOM_1.0_pim.root

mv -vi _data_${BATCH}_MOM_1.0_pim.root data_${BATCH}_MOM_1.0_pim.root 
mv -vi _simpar_${BATCH}_MOM_1.0_pim.root simpar_${BATCH}_MOM_1.0_pim.root 
mv -vi _FairRunInfo_data_${BATCH}_MOM_1.0_pim.root FairRunInfo_data_${BATCH}_MOM_1.0_pim.root 

mv -vi data_${BATCH}_MOM_1.0_pip.root data_FULLPANDA_${BATCH}_MOM_1.0_pip.root 
mv -vi simpar_${BATCH}_MOM_1.0_pip.root simpar_FULLPANDA_${BATCH}_MOM_1.0_pip.root 
mv -vi FairRunInfo_data_${BATCH}_MOM_1.0_pip.root FairRunInfo_data_FULLPANDA_${BATCH}_MOM_1.0_pip.root 

mv -vi _data_${BATCH}_MOM_1.0_pip.root data_${BATCH}_MOM_1.0_pip.root 
mv -vi _simpar_${BATCH}_MOM_1.0_pip.root simpar_${BATCH}_MOM_1.0_pip.root 
mv -vi _FairRunInfo_data_${BATCH}_MOM_1.0_pip.root FairRunInfo_data_${BATCH}_MOM_1.0_pip.root 

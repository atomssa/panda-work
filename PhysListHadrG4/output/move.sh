#!/bin/bash

BATCH=$1

#mv -v data_${BATCH}_MOM_1.0_pip.root _data_${BATCH}_MOM_1.0_pip.root
mv -v data_${BATCH}_MOM_1.0_pim.root _data_${BATCH}_MOM_1.0_pim.root
#mv -v simpar_${BATCH}_MOM_1.0_pip.root _simpar_${BATCH}_MOM_1.0_pip.root
mv -v simpar_${BATCH}_MOM_1.0_pim.root _simpar_${BATCH}_MOM_1.0_pim.root
#mv -v FairRunInfo_data_${BATCH}_MOM_1.0_pip.root _FairRunInfo_data_${BATCH}_MOM_1.0_pip.root
mv -v FairRunInfo_data_${BATCH}_MOM_1.0_pim.root _FairRunInfo_data_${BATCH}_MOM_1.0_pim.root


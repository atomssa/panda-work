#!/usr/bin/python

import sys
import my_utils
import full_sim
import evtgen

for i in range(1000,1100):

    my_utils.uniq_id = i
    my_utils.dbg_msg( "run_gen_series: doing uniq id %d" % my_utils.uniq_id)

    if ( my_utils.uniq_id % 3 == 0):
        my_utils.dpm_pbar_lab_mom = 5.513
    elif ( my_utils.uniq_id % 3 == 1):
        my_utils.dpm_pbar_lab_mom = 10.0
    else:
        my_utils.dpm_pbar_lab_mom = 20.0
    
    #delete the original sim file, since it's a disk space killer
    evtgen.generate()



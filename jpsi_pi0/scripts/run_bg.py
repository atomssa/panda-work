#!/usr/bin/python

import sys
import my_utils
import full_sim

my_utils.sim_type = my_utils.sim_bg
#my_utils.test_run = False

plab_values = [5.513, 8, 12]
#nevt = [1000, 1000, 1000]
nevt = [2000000, 20000000, 100000000]
start=int(sys.argv[1])
number=int(sys.argv[2])
series = range(start, start+number)

for iplab in range(0, len(plab_values)):
    for i in series:
        my_utils.uniq_id = i
        my_utils.pbar_lab_mom = plab_values[iplab]
        my_utils.dpm_nevt_per_file = nevt[iplab]
        my_utils.dbg_msg("run_series: doing uniq id %d" % my_utils.uniq_id)
        full_sim.run()

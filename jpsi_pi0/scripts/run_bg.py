#!/usr/bin/python

import sys
import my_utils
import full_sim

my_utils.sim_type = my_utils.sim_bg
#my_utils.test_run = False
#my_utils.delete_unfiltered_dpm = False

iplab=int(sys.argv[1]) # 0=>(p=5.513), 1=>(p=8.0), 2=>(p=12.0)
if iplab < 0 or iplab > 2:
    print("1st argument should be between 1 and 3. Supplied = %d" % iplab)
    exit(0)

#nevt = [10000, 40000, 100000]
nevt = [2000000, 8000000, 40000000]
start=int(sys.argv[2])
number=int(sys.argv[3])
series = range(start, start+number)

#for iplab in range(0, len(my_utils.plab_values)):
for i in series:
    my_utils.uniq_id = i
    my_utils.pbar_lab_mom = my_utils.plab_values[iplab]
    my_utils.dpm_nevt_per_file = nevt[iplab]
    my_utils.dbg_msg("run_series: doing uniq id %d" % my_utils.uniq_id)
    full_sim.run()

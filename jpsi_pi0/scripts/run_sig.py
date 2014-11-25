#!/usr/bin/python

import sys
import my_utils
import full_sim

my_utils.sim_type = my_utils.sim_sig
#my_utils.test_run = True;

_sim_type=int(sys.argv[1]) # 1=>(p=5.513), 2=>(p=8.0), 3=>(p=12.0)
if _sim_type < 1 or _sim_type > 3:
    print("1st argument should be between 1 and 3. Supplied = %d" % _sim_type)
start=int(sys.argv[2])
number=int(sys.argv[3])
series = range(start, start+number)
for i in series:
    my_utils.uniq_id = i
    my_utils.dbg_msg("run_series: doing uniq id %d" % my_utils.uniq_id)
    my_utils.sim_type = _sim_type
    my_utils.pbar_lab_mom = my_utils.plab_values[_sim_type-1]
    full_sim.run()

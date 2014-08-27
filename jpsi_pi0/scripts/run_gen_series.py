#!/usr/bin/python

import sys
import my_utils
import full_sim
import evtgen

#for i in range(1000,1100):
#    my_utils.uniq_id = i
#    my_utils.dbg_msg( "run_gen_series: doing uniq id %d" % my_utils.uniq_id)
#    if ( my_utils.uniq_id % 3 == 0):
#        my_utils.dpm_pbar_lab_mom = 5.513
#    elif ( my_utils.uniq_id % 3 == 1):
#        my_utils.dpm_pbar_lab_mom = 10.0
#    else:
#        my_utils.dpm_pbar_lab_mom = 20.0
#    evtgen.generate()

#my_utils.dpm_pbar_lab_mom = 5.513
#my_utils.test_run = False
#s1 = range(1000, 1100, 3)
#s2 = range(1001, 1100, 3)
#s = []
#for idx in range(0,len(s1)-1):
#    s.append(s1[idx])
#    if (idx < len(s2)):
#        s.append(s2[idx])
#print s
#for i in s:
#    my_utils.uniq_id = i
#    my_utils.dbg_msg( "run_gen_series: doing uniq id %d" % my_utils.uniq_id)
#    evtgen.generate()

plabs = [8.0, 14.0]
nevts = [10000000, 50000000]
plab_index = int(sys.argv[1])
my_utils.sim_type = my_utils.sim_bg
my_utils.dpm_pbar_lab_mom = plabs[plab_index]
my_utils.dpm_nevt_per_file = nevts[plab_index]
my_utils.test_run = False

start=int(sys.argv[2])
number=int(sys.argv[3])
series = range(start, start+number)
for i in series:
    my_utils.uniq_id = i
    my_utils.dbg_msg( "run_gen_series: doing uniq id %d at energy %f" % (my_utils.uniq_id,my_utils.dpm_pbar_lab_mom))
    evtgen.generate()

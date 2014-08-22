#!/usr/bin/python

import sys
import my_utils
import full_sim

#start=int(sys.argv[1])
#number=int(sys.argv[2])
#series = range(start, start+number)
#for i in series:
#    my_utils.uniq_id = i
#    my_utils.dbg_msg("run_series: doing uniq id %d" % my_utils.uniq_id)
#    full_sim.run()

series = [0, 3]
for i in range(9, 21):
    series2.append(i)
print "series: %s" % series

for i in series:
    my_utils.uniq_id = i
    my_utils.dbg_msg("run_series: doing uniq id %d" % my_utils.uniq_id)
    full_sim.run()

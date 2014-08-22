#!/usr/bin/python

import sys
import my_utils
import tmp_full_sim

series = [0, 3]
print "series: %s" % series

for i in series:
    my_utils.uniq_id = i
    my_utils.dbg_msg("tmp_run_series: doing uniq id %d" % my_utils.uniq_id)
    tmp_full_sim.run()

series2= []
for i in range(9, 21):
    series2.append(i)

print "series2: %s" % series2

for i in series2:
    my_utils.uniq_id = i
    my_utils.dbg_msg("tmp_run_series: doing uniq id %d" % my_utils.uniq_id)
    tmp_full_sim.run()



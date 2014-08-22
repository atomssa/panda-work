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

series = []

arg=int(sys.argv[1])

if arg == 1:
    series.append(0)
    series.append(3)
    series.append(11)
    for i in range(14, 30):
        series.append(i)
if arg == 2:
    series.append(31)
    series.append(34)
    series.append(35)
    series.append(38)
    series.append(39)    
    series.append(41)
    series.append(44)
    for i in range(45, 61):
        series.append(i)
        
print "series: %s" % series

my_utils.test_run = False

for i in series:
    my_utils.uniq_id = i
    my_utils.dbg_msg("run_series: doing uniq id %d" % my_utils.uniq_id)
    full_sim.run()

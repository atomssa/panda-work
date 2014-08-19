#!/usr/bin/python

import sys
import full_sim

def run_series(uniq_id, start):
    series = range(start, start+3)
    for i in series:
        uniq_id = i
        print "run_series: doing uniq id %d" % uniq_id
        full_sim.run(uniq_id)

if __name__ == "__main__":
    start=int(sys.argv[1])
    uniq_id = 0
    run_series(uniq_id, start)

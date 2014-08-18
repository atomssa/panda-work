#!/usr/bin/python

import subprocess
import sys,os
import utils

nbatch = int(4)
batch_id = int(sys.argv[1])

def execute_root(args):
    subprocess.Popen()

def execute_step(step, uniq_id):
    log_dir = prep_dir("output", step)
    out_dir = prep_dir("log", step)
    log_file = "%s/%s.log" % (log_dir, file_tag)
    out_file = "%s/%s.root" % (out_dir, file_tag)
    batch_dir = prep_batch_dir(file_id)
    os.chdir(batch_dir)
    excute_root(...args)
    os.rename("%s_complete.root" % step, out_file)
    
for job in range(0,3):
    uniq_id = job*nbatch + batch
    file_tag = "%s_%d" % (utils.proc_tag, uniq_id)
    print "launching simulation for uniq_id = %d, file_tag= %s" % (uniq_id, file_tag)

    # input file for simulation
    dpm_file "%s/output/evt/%s.root" % (utils.base_dir, file_tag)
    batch_dir = utils.prep_batch_dir( file_id )

    for each step in util.steps:
        execute_step(step, uniq_id)

    out_dir = prep_dir("output", "simpar")
    os.rename("simparams.root", )

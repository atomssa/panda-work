#!/usr/bin/python

import subprocess
import sys,os
import my_utils
import evtgen

def execute_root(step, log, prev_out, prev_step, args=None):
    print "========= executing root with args %s " % args
    os.chdir(my_utils.prep_batch_dir( uniq_id ))
    my_utils.soft_link(prev_out, "%s_complete.root" % prev_step)
    proc = subprocess.Popen("root -b -q %s" % args,shell=True, stdout=log, stderr=log)
    proc.wait()
    
def execute_step(step, uniq_id, prev_out, prev_step):
    print "====== executing step %s with uniq id %d " % (step,uniq_id)
    (out_file,log_file,dummy) = my_utils.file_names(uniq_id, step)
    print step
    print out_file
    print log_file    
    if not os.path.exists(out_file):
        if (step != "sim"):
            root_args =  "%s/%s_complete.C" % (my_utils.macro_dir,step)
            execute_root(step, open(log_file,"w"), prev_out, prev_step, root_args)
        else:
            root_args = "\'%s/%s_complete.C(\"%s\")\'" % (my_utils.macro_dir, step, prev_out)
            execute_root(step, open(log_file,"w"), prev_out, prev_step, root_args)
        os.rename("%s_complete.root" % step, out_file)
    return (step,out_file)

def run(uniq_id):
    print "==== launching simulation uniq_id %d ===" % uniq_id
    prev_step = "filt"
    prev_out = evtgen.generate(uniq_id)
    for step in my_utils.simulation_steps:
        (prev_step,prev_out) = execute_step(step, uniq_id, prev_out, prev_step)
    (out_file,dummy,dummy) = my_utils.file_names(uniq_id, "par")
    os.rename("simparams.root", out_file)

if __name__ == "__main__":
    uniq_id = int(sys.argv[1])
    run(uniq_id)

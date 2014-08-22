#!/usr/bin/python

import subprocess
import sys,os
import my_utils
import evtgen

def test_exec(step,log,prev_step,prev_out,args):
    out = "%s_complete.root" % step
    my_utils.cd_to_batch_dir()
    my_utils.dbg_msg("writing to %s"%out)
    checklist = []
    checklist.append("prev. step = %s\n"%prev_step)
    checklist.append("prev. file %s exsists? : %s\n"%(prev_out, os.path.exists(prev_out)))
    checklist.append("plab= %f\n"%my_utils.dpm_pbar_lab_mom)
    checklist.append("args= %s\n"%args)
    checklist.append("output from current step= %s\n"%out)
    checklist.append("current dir content= %s\n"%os.listdir(os.getcwd()))
    with open(out,"w") as fout:
        for check in checklist: fout.write(check)
    with open('simparams.root','w') as fout:
        fout.write("mehmehmeh")
    
def execute_root(step, log, prev_out, prev_step, args=None):
    my_utils.dbg_msg("========= executing root with args %s " % args)
    my_utils.cd_to_batch_dir()
    my_utils.soft_link(prev_out, "%s_complete.root" % prev_step)
    proc = subprocess.Popen("root -b -q %s" % args,shell=True, stdout=log, stderr=log)
    proc.wait()
    #test_exec(step,log,prev_step,prev_out,args)
    #os.unlink("%s_complete.root"%prev_step)
    
def execute_step(step, prev_out, prev_step):
    my_utils.dbg_msg( "====== executing step %s with uniq id %d " % (step,my_utils.uniq_id) )
    (out_file,log_file,dummy) = my_utils.file_names(step)
    if (step != "sim"):
        root_args =  "%s/%s_complete.C" % (my_utils.macro_dir,step)
        with open(log_file,"w") as log:
            execute_root(step, log, prev_out, prev_step, root_args)
    else:
        root_args = "\'%s/%s_complete.C(\"%s\")\'" % (my_utils.macro_dir, step, prev_out)
        with open(log_file,"w") as log:
            execute_root(step, log, prev_out, prev_step, root_args)
    my_utils.move_file("%s_complete.root" % step, out_file)
    return (step,out_file)

def run():
    if (my_utils.uniq_id==None):
        my_utils.dbg_msg( "my_utils.uniq_id should be set before calling this function")
        sys.exit(-1)
    else:
        my_utils.dbg_msg( "==== launching simulation uniq_id %d ===" % my_utils.uniq_id)

    prev_step = "filt"
    prev_out = evtgen.generate()

    for step in my_utils.simulation_steps:
        (prev_step,prev_out) = execute_step(step, prev_out, prev_step)

    (out_file,dummy,dummy) = my_utils.file_names("par")
    my_utils.move_file("simparams.root", out_file)

if __name__ == "__main__":
    my_utils.uniq_id = int(sys.argv[1])
    run()

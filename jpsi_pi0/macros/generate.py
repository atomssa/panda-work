#!/usr/bin/python

import subprocess;
from random import randint
import sys
import utils

entry_chan = "pbar_p"
exit_chan = "pip_pim_pi0"

nevt_per_file = 2000000

pbar_lab_mom = 5.513

proc_selection = 0 # inelastic only

part_status_mod = "0 0"

nbatch = int(4)
batch_id = int(sys.argv[1])

for job in range(0,3):

    file_id = job*nbatch+ batch_id

    print "launching process for file_id = %s" % file_id

    base_dir = "/vol0/panda/work/jpsi_pi0/"
    batch_dir = ".tmp_%d" %file_id
    dpm_out = "Background-micro.root"
    proc_tag = "%s_%s" %  (entry_chan, exit_chan)
    
    out_dir = "%s/output/evt" % base_dir
    log_dir = "%s/log/evt" % base_dir

    utils.mkdir_if_absent(out_dir)
    utils.mkdir_if_absent(log_dir)
    
    out_file = "%s/%s_%d.root" % (out_dir, proc_tag, file_id)
    log_file = "%s/%s_%d.log" % (log_dir, proc_tag, file_id)
    in_file =  "%s/%s_%d.input" % (log_dir, proc_tag, file_id)
    

    seed = randint(0, 1000000)
    sub_com = [];
    if os.path.exists(batch_dir):
        os.rmtree(batch_dir)
    input_file = open("%s/input" % batch_dir)
    input_file.writelines([seed,pbar_lab_mom,proc_selection, part_status_mod, nevt_per_file])
    print input_file.read()

    os.chdir(batch_dir)

    sub_com.append("DPMGen < %s > %s 2>&1" % (in_file, log_file))

    
    sub_com.append("mv %s %s" %( dpm_out, out_file))
    sub_com.append("cd ../; rmdir %s" % batch_dir)
    cmd = ';'.join(sub_com)
    #print "cmd = %s"%cmd 
    
    proc = subprocess.Popen(cmd, shell=True);
    proc.wait()

    print "...... done lauching for %s " % file_id

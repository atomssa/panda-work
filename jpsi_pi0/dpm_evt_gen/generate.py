#!/usr/bin/python

import subprocess;
from random import randint
import sys

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

    batch_dir = ".tmp_%d" %file_id
    dpm_out = "Background-micro.root"
    proc_tag = "%s_%s" %  (entry_chan, exit_chan)
    out_dir = "/vol0/panda/work/jpsi_pi0/output"
    log_dir = "/vol0/panda/work/jpsi_pi0/log"
    out_file = "%s/%s_%d.root" % (out_dir, proc_tag, file_id)
    log_file = "%s/%s_%d.log" % (log_dir, proc_tag, file_id)
    in_file =  "%s/%s_%d.input" % (log_dir, proc_tag, file_id)

    seed = randint(0, 1000000)
    sub_com = [];
    sub_com.append("[ -e %s ] || mkdir %s; cd %s; rm -vrf *" % (batch_dir, batch_dir, batch_dir))
    sub_com.append("[ -e %s ] && rm -vf %s " % (log_file, log_file))
    sub_com.append("[ -e %s ] && rm -vf %s " % (in_file, in_file))
    sub_com.append("touch %s" % in_file)
    sub_com.append("echo %d >> %s" % (seed, in_file))
    sub_com.append("echo %f >> %s" % (pbar_lab_mom, in_file))
    sub_com.append("echo %d >> %s" % (proc_selection, in_file))
    sub_com.append("echo %s >> %s" % (part_status_mod, in_file))
    sub_com.append("echo %d >> %s" % (nevt_per_file, in_file))
    sub_com.append("more %s" % in_file)
    sub_com.append("DPMGen < %s > %s 2>&1" % (in_file, log_file))
    sub_com.append("mv %s %s" %( dpm_out, out_file))
    sub_com.append("cd ../; rmdir %s" % batch_dir)
    cmd = ';'.join(sub_com)
    #print "cmd = %s"%cmd 
    
    proc = subprocess.Popen(cmd, shell=True);
    proc.wait()

    print "...... done lauching for %s " % file_id

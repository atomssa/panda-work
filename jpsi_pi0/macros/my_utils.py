#!/usr/bin/python

import os
from random import randint

uniq_id = 0

base_dir = "/vol0/panda/work/jpsi_pi0"
macro_dir = "%s/macros" % base_dir

entry_chan = "pbar_p"
exit_chan = "pip_pim_pi0"
proc_tag = "%s_%s" % (entry_chan, exit_chan)

dpm_nevt_per_file = 2000000
dpm_pbar_lab_mom = 5.513
dpm_proc_selection = 0
dpm_part_status_mod = "0 0"
dpm_default_out = "Background-micro.root"

simulation_steps = [ "sim", "digi", "reco", "pid" ]

def get_rnd_seed():
    return randint(0, 1000000)

def file_tag(uniq_id):
    return "%s_%d" % (proc_tag, uniq_id)

def file_names(uniq_id,step):
    out_dir = prep_dir("output", step)
    log_dir = prep_dir("log", step)
    tmp = file_tag(uniq_id)
    return("%s/%s.root" % (out_dir, file_tag(uniq_id)),
           "%s/%s.log" % (log_dir, file_tag(uniq_id)),
           "%s/%s.in" % (log_dir, file_tag(uniq_id)))

def soft_link(src,dest):
    print "ln -sf %s %s" % (src,dest)
    if not os.path.exists(dest):
        os.symlink(src,dest)
    else: # it may be stale
        os.unlink(dest)
        os.symlink(src,dest)

def prep_dir(dir_type, step):
    the_dir = "%s/%s/%s" % (base_dir,dir_type,step)
    _mkdir_if_absent(the_dir)
    return the_dir

def prep_batch_dir(identifier):
    the_dir = "%s/.batch/tmp_%d" % (base_dir, identifier)
    _mkdir_if_absent(the_dir)
    #_empty_dir_contents(the_dir)
    return the_dir

# this should be called with extreme care!!!
# probably only for temporary batch directories    
# even then maybe not :o
#def _empty_dir_contents(_dir):
#    file_list = os.listdir(_dir)
#    for file_name in file_list:
#        os.remove(_dir+"/"+file_name)

def _mkdir_if_absent(the_dir):
    if not os.path.exists(the_dir):
        print "created dir: %s" % the_dir
        os.makedirs(the_dir)
    #else:
    #    print "[info] dir %s alrady exists... " % the_dir

    

    

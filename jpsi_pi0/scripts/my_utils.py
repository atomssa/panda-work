#!/usr/bin/python

import os,sys,shutil,datetime
from random import randint

uniq_id = None

debug = True

base_dir = "/vol0/panda/work/jpsi_pi0"
macro_dir = "%s/macros" % base_dir
batch_dir = None

entry_chan = "pbar_p"
exit_chan = "pip_pim_pi0"
proc_tag = "%s_%s" % (entry_chan, exit_chan)

dpm_nevt_per_file = 2000000
dpm_pbar_lab_mom = 5.513
dpm_proc_selection = 0
dpm_part_status_mod = "0 0"
dpm_default_out = "Background-micro.root"
filter_default_out = "Background-nano.root"

delete_unfiltered_dpm = True

simulation_steps = [ "sim", "digi", "reco", "pid" ]

def dbg_msg(msg):
    if (debug):
        print msg

def get_rnd_seed():
    return randint(0, 1000000)

def file_tag():
    if uniq_id != None:
        return "%s_%d" % (proc_tag, uniq_id)
    else:
        dbg_msg("file_tag called with uniq_id not set yet. calling sys.exit(-1)")
        sys.exit(-1)

def file_names(step):
    out_dir = "%s/output/%s" % (base_dir,step)
    log_dir = "%s/log/%s" % (base_dir,step)

    # small time logging logic
    timer_file = "%s/log/timer/%s.log"%(base_dir,file_tag())
    if not os.path.exists(os.path.dirname(timer_file)):
        os.makedirs(os.path.dirname(timer_file))
    with open(timer_file,'a') as timer:
        timer.write("Start(%s): %s\n"%(step, datetime.datetime.now()))
        timer.flush()

    return("%s/%s.root" % (out_dir, file_tag()),
           "%s/%s.log" % (log_dir, file_tag()),
           "%s/%s.in" % (log_dir, file_tag()))
    
def soft_link(src,dest):
    dbg_msg("Executing ln -sf %s %s" % (src,dest))
    dbg_msg("soft_link: listing cwd %s" % os.getcwd())
    dbg_msg("%s" % os.listdir(os.getcwd()))
    if os.path.exists(src):
        if not os.path.exists(dest): # dest path doesn't exist
            if not os.path.dirname(dest) != '' and os.path.exists(os.path.dirname(dest)): #create destination's parent dir if it doesnt exist and is not current dir
                dbg_msg("[W] destination dir %s doesnt exist, creating it" % dest )
                os.makedirs(os.path.dirname(dest))
            if os.path.islink(dest): # if the path is a link and link is broken, os.path.exists() returns false, so remove it
                dbg_msg("unlinking %s"%dest)
                os.unlink(dest)
            os.symlink(src,dest)
        else: # destination exists
            if os.path.isdir(dest): # dest path is directory
                os.symlink(src,"%s/%s"%(dest,os.path.basename(src)))
            else: 
                dbg_msg("[W] overwriting destination file: %s" % dest )
                os.unlink(dest)
                os.symlink(src,dest)
    else:
        dbg_msg("soft_link failure: source file %s doesn't exist" % src)

def copy_file(src,dest):
    dbg_msg("copy_file: listing cwd %s" % os.getcwd())
    dbg_msg("Executing copy -f %s %s" % (src,dest))
    dbg_msg("%s" % os.listdir(os.getcwd()))
    if os.path.exists(src): # source file exists
        if not os.path.exists(dest): # dest path doesnt exist
            if not os.path.dirname(dest) != '' and os.path.exists(os.path.dirname(dest)): #create destination's parent dir if it doesnt exist and is not current dir
                dbg_msg("[W] destination dir %s doesnt exist, creating it" % dest )
                os.makedirs(os.path.dirname(dest))
            shutil.copy2(src,dest)
        else: # dest path exists too
            if os.path.isdir(dest): # dest path is directory
                dbg_msg("dest= %s/%s"%(dest,os.path.basename(src)))
                shutil.copy2(src,"%s/%s"%(dest,os.path.basename(src)))
            else: # overwrite
                dbg_msg("[W] overwriting destination file: %s" % dest )
                shutil.copy2(src,dest)
    else:
        dbg_msg("copy_file failure: source file %s doesn't exist" % src)

def move_file(src,dest):
    dbg_msg("Executing move -f %s %s" % (src,dest))
    dbg_msg("move_file: listing cwd %s" % os.getcwd())
    dbg_msg("%s" % os.listdir(os.getcwd()))
    if os.path.exists(src): # source file exists
        if not os.path.exists(dest): # dest path doesnt exist
            if not os.path.dirname(dest) != '' and os.path.exists(os.path.dirname(dest)): #create destination's parent dir if it doesnt exist and is not current dir
                dbg_msg("[W] destination dir %s doesnt exist, creating it" % dest )
                os.makedirs(os.path.dirname(dest))
            if os.path.islink(dest):
                os.unlink(dest)
            os.rename(src,dest)
        else: # dest path exists too
            if os.path.isdir(dest): # dest path is directory
                dbg_msg("dest= %s/%s"%(dest,os.path.basename(src)))
                os.rename(src,"%s/%s"%(dest,os.path.basename(src)))
            else: # overwrite
                dbg_msg("[W] overwriting destination file: %s" % dest )
                os.rename(src,dest)
    else:
        dbg_msg("move_file failure: source file %s doesn't exist" % src)
        
def copy_file_to_batch_dir(src):
    if uniq_id != None:
        batch_dir = "%s/.batch/tmp_%d" % (base_dir, uniq_id)        
        copy_file(src,batch_dir)
    else:
        dbg_msg("cd_to_batch_dir called with uniq_id not set yet. calling sys.exit(-1)")
        sys.exit(-1)

def cd_to_batch_dir():
    if uniq_id != None:
        batch_dir = "%s/.batch/tmp_%d" % (base_dir, uniq_id)
        if not os.path.exists(batch_dir):
            os.makedirs(batch_dir)
        os.chdir(batch_dir)
    else:
        dbg_msg("cd_to_batch_dir called with uniq_id not set yet. calling sys.exit(-1)")
        sys.exit(-1)


#def _mkdir_if_absent(the_dir):
#    dbg_msg("[I] creating dir %s if it doesnt exis" % the_dir)
#    if not os.path.exists(the_dir):
#        os.makedirs(the_dir)
#    else:
#        dbg_msg("[I] dir %s alrady exists... " % the_dir)

# this should be called with extreme care!!!
# probably only for temporary batch directories    
# even then maybe not :o
#def _empty_dir_contents(_dir):
#    file_list = os.listdir(_dir)
#    for file_name in file_list:
#        os.remove(_dir+"/"+file_name)

    

    

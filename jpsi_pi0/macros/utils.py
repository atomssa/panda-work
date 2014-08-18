import os

base_dir = "/vol0/panda/work/jpsi_pi0"
entry_chan = "pbar_p"
exit_chan = "pip_pim_pi0"
proc_tag = "%s_%s" % (entry_chan, exit_chan)

steps = [ "filt", "sim", "dig", "rec", "pid"]

def mkdir_if_absent(_dir):
    if not os.path.exists(_dir):
        print "created dir: %s" % _dir
        os.path.makedirs(_dir)
    else:
        print "[info] dir %s alrady exists... " % _dir

def prep_dir(dir_type, step):
    _dir = "%s/%s/%s" % (base_dir,dir_type,step)
    mkdir_if_absent(out_dir)
    return _dir

# this should be called with extreme care!!!
# probably only for temporary batch directories    
# even then maybe not :o
def empty_dir_contents(_dir):
    fileList = os.listdir(_dir)
    for file_name in file_list:
        os.remove(_dir+"/"+file_name)
    
def prep_batch_dir(base, identifier):
    _dir = "%s/.batch/tmp_%d" % (base, identifier)
    utils.mkdir_if_absent(_dir)
    utils.empty_dir_contents(_dir)
    return _dir
    

    

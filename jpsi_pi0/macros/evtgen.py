#!/usr/bin/python

import shutil,sys, os, subprocess
import my_utils

def execute(uniq_id,cmd,log,batch_input):
    os.chdir(my_utils.prep_batch_dir(uniq_id))
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log, stdin=batch_input);
    proc.wait()
    
def write_lines(file_handle, lines):
    for line in lines:
        file_handle.write(line + " \n")
    
def dpm_in(uniq_id,in_file_name):
    lines = [ str(my_utils.get_rnd_seed()), str(my_utils.dpm_pbar_lab_mom),
              str(my_utils.dpm_proc_selection), my_utils.dpm_part_status_mod, str(my_utils.dpm_nevt_per_file)]
    write_lines(open(in_file_name,"w"),lines)
    return in_file_name

def run_dpm(uniq_id):
    print "========= running dpm for uniq_id = %d " % uniq_id
    (out_file,log_file,in_file) = my_utils.file_names(uniq_id,"evt")
    if not os.path.exists(out_file):
        execute(uniq_id, "DPMGen", open(log_file,'w'), open(dpm_in(uniq_id,in_file),"r"))
        os.rename(my_utils.dpm_default_out, out_file)
    return out_file
    
def root_in(in_file_name,dpm_file_name,filt_file_name):
    shutil.copy2("%s/dpm_filter.C" % my_utils.macro_dir, my_utils.prep_batch_dir(uniq_id))
    lines = [".L %s/dpm_filter.C++" % my_utils.macro_dir, "dpm_filter(\"%s\",\"%s\")" % (dpm_file_name,filt_file_name), ".q"]
    write_lines(open(in_file_name,"w"),lines)
    return in_file_name

def filter_dpm(uniq_id, dpm_out):
    print "========= filtering for uniq_id = %d " % uniq_id
    (out_file,log_file,in_file) = my_utils.file_names(uniq_id,"filt")
    if not os.path.exists(out_file):
        execute(uniq_id, "root -b", open(log_file,'w'), open(root_in(in_file, dpm_out, out_file),"r"))
    return out_file
        
def generate(uniq_id):
    print "====== launching event generation for uniq_id = %d" % uniq_id
    #dpm_out = run_dpm(uniq_id)
    #filt_out = filter_dpm(uniq_id, dpm_out)
    #return filt_out
    
if __name__ == "__main__":
    generate(int(sys.argv[1]))

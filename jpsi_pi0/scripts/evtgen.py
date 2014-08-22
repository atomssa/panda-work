#!/usr/bin/python

import os, subprocess
import my_utils

def test_exec(log,inpipe):
    my_utils.cd_to_batch_dir()
    first_input_line = inpipe.readline().strip()
    ref = ".L %s/dpm_filter.C++" % my_utils.macro_dir
    my_utils.dbg_msg("first input line= %s" % first_input_line)
    my_utils.dbg_msg("reference= %s" % ref)
    if ( first_input_line == ref):
        my_utils.dbg_msg("its filter")
        out = my_utils.filter_default_out
        proc = subprocess.Popen(
            "cd blah; ls>%s;echo %s>> %s; echo %f>>%s"%(out,inpipe.readlines(), out, my_utils.dpm_pbar_lab_mom ,out),
            shell=True, stdout=log, stderr=log);
    else:
        my_utils.dbg_msg("its dpm")
        out = my_utils.dpm_default_out
        proc = subprocess.Popen(
            "cd blah; ls>%s;echo %s>> %s;echo plab %f>>%s"%(out,inpipe.readlines(), out, my_utils.dpm_pbar_lab_mom, out),
            shell=True, stdout=log, stderr=log);
    proc.wait()
        
def execute(cmd,log,inpipe):
    my_utils.cd_to_batch_dir()
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log, stdin=inpipe);
    proc.wait()
    #test_exec(log,inpipe)
    
def write_lines(file_handle, lines):
    for line in lines:
        file_handle.write(line + " \n")
    
def dpm_in(in_file):
    lines = [ str(my_utils.get_rnd_seed()), str(my_utils.dpm_pbar_lab_mom),
              str(my_utils.dpm_proc_selection), my_utils.dpm_part_status_mod, str(my_utils.dpm_nevt_per_file)]
    with open(in_file,"w") as inpipe:
        write_lines(inpipe,lines)
    return in_file

def run_dpm():
    my_utils.dbg_msg("========= running dpm for uniq_id = %d " % my_utils.uniq_id)
    (out_file,log_file,in_file) = my_utils.file_names("evt")
    with open(log_file,'w') as log, open(dpm_in(in_file),"r") as inpipe:
        execute("DPMGen", log, inpipe)
    my_utils.move_file(my_utils.dpm_default_out, out_file)
    return out_file
    
def root_in(in_file,dpm_file):
    src_file = "%s/dpm_filter.C" % my_utils.macro_dir
    my_utils.copy_file_to_batch_dir(src_file)
    lines = [".L %s/dpm_filter.C++" % my_utils.macro_dir, "dpm_filter(\"%s\")" % dpm_file, ".q"]
    with open(in_file,"w") as inpipe:
        write_lines(inpipe,lines)
    return in_file

def filter_dpm(dpm_out):
    my_utils.dbg_msg("========= filtering for uniq_id = %d " % my_utils.uniq_id)
    (out_file,log_file,in_file) = my_utils.file_names("filt")
    with open(log_file,'w') as log, open(root_in(in_file, dpm_out),"r") as inpipe:
        execute("root -b", log, inpipe)
    my_utils.move_file(my_utils.filter_default_out,out_file)
    return out_file
    
def generate():
    if my_utils.uniq_id == None:
        my_utils.dbg_msg("my_utils.uniq_id should be set before calling this function")
    else:
        my_utils.dbg_msg("====== launching event generation for uniq_id = %d" % my_utils.uniq_id)

    (dummy,log_file,dummy) = my_utils.file_names("timer")
    
    dpm_out = run_dpm()
    filt_out = filter_dpm(dpm_out)

    # make delete_unfiltered_dpm golbal in my_utils at the next ocasion 
    if (my_utils.delete_unfiltered_dpm):
        if os.path.exists(dpm_out):
            os.remove(dpm_out)
    
    return filt_out
    
if __name__ == "__main__":
    my_utils.uniq_id = int(sys.argv[1])
    generate(False)

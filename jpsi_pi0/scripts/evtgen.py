#!/usr/bin/python

import os, subprocess
import my_utils

def test_exec(log,inpipe):
    my_utils.cd_to_batch_dir()
    first_input_line = inpipe.readline().strip()
    ref = ".L %s/dpm_filter.C++" % my_utils.macro_dir()
    my_utils.dbg_msg("first input line= %s" % first_input_line)
    my_utils.dbg_msg("reference= %s" % ref)
    if ( first_input_line == ref):
        my_utils.dbg_msg("its filter")
        out = my_utils.filter_default_out
        proc = subprocess.Popen(
            "ls>%s;echo %s>> %s; echo %f>>%s"%(out,inpipe.readlines(), out, my_utils.pbar_lab_mom ,out),
            shell=True, stdout=log, stderr=log);
    else:
        my_utils.dbg_msg("its dpm")
        out = my_utils.dpm_default_out
        proc = subprocess.Popen(
            "ls>%s;echo %s>> %s;echo plab %f>>%s"%(out,inpipe.readlines(), out, my_utils.pbar_lab_mom, out),
            shell=True, stdout=log, stderr=log);
    proc.wait()

def execute(cmd,log,inpipe):
    my_utils.cd_to_batch_dir()
    if not my_utils.test_run:
        proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log, stdin=inpipe);
        proc.wait()
    else:
        test_exec(log,inpipe)

def write_lines(file_handle, lines):
    for line in lines:
        file_handle.write(line + " \n")

def dpm_in(in_file):
    lines = [ str(my_utils.get_rnd_seed()), str(my_utils.pbar_lab_mom),
              str(my_utils.dpm_proc_selection), my_utils.dpm_part_status_mod, str(my_utils.dpm_nevt_per_file)]
    with open(in_file,"w") as inpipe:
        write_lines(inpipe,lines)
    return in_file

def run_dpm():
    my_utils.dbg_msg("========= running dpm for uniq_id = %d " % my_utils.uniq_id)
    (out_file,log_file,in_file) = my_utils.file_names("evt")
    if not os.path.exists(out_file):
        with open(log_file,'w') as log, open(dpm_in(in_file),"r") as inpipe:
            execute("DPMGen", log, inpipe)
    else:
        my_utils.dbg_msg("%s exists. will skip this step (evt)" % out_file)
    return out_file

def root_in(in_file,dpm_file):
    src_file = "%s/dpm_filter.C" % my_utils.macro_dir()
    my_utils.copy_file_to_batch_dir(src_file)
    lines = [".L %s/dpm_filter.C++" % my_utils.macro_dir(), "dpm_filter(\"%s\")" % dpm_file, ".q"]
    with open(in_file,"w") as inpipe:
        write_lines(inpipe,lines)
    return in_file

def filter_dpm():
    my_utils.dbg_msg("========= filtering for uniq_id = %d " % my_utils.uniq_id)
    (out_file,log_file,in_file) = my_utils.file_names("filt")
    if not os.path.exists(out_file):
        with open(log_file,'w') as log, open(root_in(in_file, my_utils.dpm_default_out),"r") as inpipe:
            execute("root -b", log, inpipe)
        my_utils.move_file(my_utils.filter_default_out,out_file)
    else:
        my_utils.dbg_msg("%s exists. will skip this step (filt)" % out_file)
    if not my_utils.delete_unfiltered_dpm:
        (evt_out_file,evt_log_file,evt_in_file) = my_utils.file_names("evt")
        my_utils.move_file(my_utils.dpm_default_out, evt_out_file)
    else:
        if os.path.exists(my_utils.dpm_default_out):
            os.remove(my_utils.dpm_default_out)
    return out_file

def generate():
    if my_utils.uniq_id == None:
        my_utils.dbg_msg("my_utils.uniq_id should be set before calling this function")
    else:
        my_utils.dbg_msg("====== launching event generation for uniq_id = %d" % my_utils.uniq_id)

    (dummy,log_file,dummy) = my_utils.file_names("timer")

    if my_utils.sim_type == my_utils.sim_bg:
        dpm_out = run_dpm()
        filt_out = filter_dpm()
    else:
        (out_file,log_file,in_file) = my_utils.file_names("filt")
        filt_out = out_file
        with open(out_file,'w') as out:
            out.write("this is just a place holder for signal simulation")

    return filt_out

if __name__ == "__main__":
    my_utils.uniq_id = int(sys.argv[1])
    generate(False)

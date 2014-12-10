#!/usr/bin/python

import os,sys,shutil,subprocess

plab = [5.513, 8., 12.]
nfile = [ [5, 5, 5], [10, 10, 10]]
ibrem = 1
ana_args  = [];

def run_ana(ii):
    args = ana_args[ii]
    logf_path = "log/ana/ana_tda_brem%d_%s_%3.1f_%d.log" % (ibrem, ("bg" if args[1]==0 else "jpsi"), plab[args[0]], args[2])
    print "logf_path= %s" % logf_path
    with open(logf_path, 'w') as logf:
        cmd = "root -b -q \"macros/ana.C(%d, %d, %d, %d, %d)\"" % (args[0], args[1], ibrem, args[2], args[3])
        print "cmd= %s" % cmd
        subprocess.Popen(cmd, shell=True, stdout=logf, stderr=logf);

if __name__ == "__main__":
    #arg1 = 0 for Bg 1 for sig
    #arg2 = plab idx
    itype = int(sys.argv[1])
    iplab = int(sys.argv[2])
    for ii in range(nfile[itype][iplab]):
        ana_args.append((iplab,itype,ii,10000))
    print ana_args;

    for ii in range(len(ana_args)):
        run_ana(ii)

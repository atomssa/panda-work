#!/usr/bin/python

import os,sys,shutil,subprocess

ana_args  = [];
ana_args.append((1,0,10000))
for i in range(10):
    ana_args.append((0,i,10000))
print ana_args;

def execute(logf, inpipe):
    proc = subprocess.Popen("root -b", shell=True, stdout=logf, stderr=logf, stdin=inpipe);

def run_ana(ii):
    args = ana_args[ii]
    logf_path = "log/ana/ana_tda_%s_%d.log" % (("bg" if args[0]==0 else "jpsi"), args[1])
    print "logf_path= %s" % logf_path
    with open(logf_path, 'w') as logf:
        cmd = "root -b \"macros/ana.C(%d,%d,%d)\"" % (args[0], args[1], args[2] )
        print "cmd= %s" % cmd
        subprocess.Popen(cmd, shell=True, stdout=logf, stderr=logf);

def count():
    sum = 0
    for i in range(1000):
        if i % 5 == 0 or i % 3 == 0 :
            sum += i
    print sum

if __name__ == "__main__":
    base = 4 * int(sys.argv[1])
    for i in range(4):
        if not (base == 8 and i == 3):
            run_ana(base+i)

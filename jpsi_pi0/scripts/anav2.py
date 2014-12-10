#!/usr/bin/python

import os,sys,shutil,subprocess,glob

nfile_dict = [
    [
        dict( [ (0, 9685),  (1, 9869),  (2, 9942),  (3, 9949),  (4, 9824),  (5, 9993),  (7, 9723),  ]),
        dict( [ (0, 10219),  (1, 10196),  (2, 10298),  (3, 10110),  (4, 10278),  (5, 10277),  (7, 10208),  ]),
        dict( [ (0, 9646),  (1, 9366),  (2, 9458),  (3, 9309),  (4, 9601),  (5, 9559),  (6, 9485),  ]),
    ],
    [
        dict( [ (0, 10000),  (1, 10000),  (10, 10000),  (11, 10000),  (12, 10000),  (13, 10000),  (14, 10000),  (15, 10000),  (16, 10000),  (17, 10000),  (18, 10000),  (19, 10000),  (2, 10000),  (3, 10000),  (4, 10000),  (5, 10000),  (6, 10000),  (7, 10000),  (8, 10000),  (9, 10000),  ]),
        dict( [ (0, 10000),  (1, 10000),  (10, 10000),  (11, 10000),  (12, 10000),  (13, 10000),  (15, 10000),  (2, 10000),  (3, 10000),  (4, 10000),  (5, 10000),  (6, 10000),  (7, 10000),  (8, 10000),  (9, 10000),  ]),
        dict( [ (0, 10000),  (1, 10000),  (10, 10000),  (2, 10000),  (3, 10000),  (4, 10000),  (5, 10000),  (6, 10000),  (7, 10000),  (8, 10000),  (9, 10000),  ]),
    ],
]

def test_dict():
    nfile = [[7, 7, 7], [20, 15, 11]]
    errcnt = 0
    for it,tt in enumerate(nfile):
        for ip,n in enumerate(tt):
            if len(nfile_dict[it][ip]) != nfile[it][ip]:
                print "ERROR: nfile[%d][%d]= %d,  EXPECTED: %d" %(it, ip, len(nfile_dict[it][ip]), nfile[it][ip])
                errcnt = errcnt+1;

    if nfile_dict[0][1][0] != 10219:
        print "ERROR: found n(bg,7)=%d expected 10208" % nfile_dict[0][1][7]
        errcnt = errcnt+1

    if errcnt == 0:
        print "TEST (number of entries) OK!"

ibrem = 1
plab = [5.513, 8., 12.]
#nfile = [ [4, 1, 1], [3, 5, 5]]
nevt = [[36000, 9000, 4500], [28040, 43013, 44525]]

def nfile(itype,iplab):
    n = 0
    nevt_tot = 0
    for idx in sorted(nfile_dict[itype][iplab].keys()):
        if (nevt_tot + nfile_dict[itype][iplab][idx]) < nevt[itype][iplab]:
            n += 1
            nevt_tot = nevt_tot + nfile_dict[itype][iplab][idx]
        else:
            break
    return n + 1

def hadd():
    for itype in range(2):
        for iplab in range(3):
            logf_path = "log/ana/hadd_anav2_tda_%s_%3.1f.log" % (("bg" if itype==0 else "jpsi"), plab[iplab])
            base = "anav2_%s_%s_plab%3.1f" % (("pip_pim" if itype==0 else "jpsi"), ("raw" if ibrem==0 else "brem"), plab[iplab])
            cmd = "[ -e test/%s.root] && rm test/%s.root; hadd test/%s.root test/%s_*.root" % (base, base, base, base)
            mrg_list = glob.glob("test/%s_*.root"%base)
            print mrg_list
            if len(mrg_list) != nfile(itype,iplab):
                print "Number of files to merge not correct %d expected" % nfile[itype][iplab]
            else:
                print cmd
                with open(logf_path, 'w') as logf:
                    subprocess.Popen(cmd, shell=True, stdout=logf);

def run_ana(args):
    logf_path = "log/ana/anav2_tda_%s_%3.1f_%d.log" % (("bg" if args[1]==0 else "jpsi"), plab[args[0]], args[3])
    print "logf_path= %s" % logf_path
    with open(logf_path, 'w') as logf:
        cmd = "root -b -q \"macros/anav2.C(%d, %d, %d, %d, %d)\"" % (args[0], args[1], args[2], args[3], args[4])
        print "cmd= %s" % cmd
        subprocess.Popen(cmd, shell=True, stdout=logf, stderr=logf);

def run_parallel():
    for itype in range(2):
        for iplab in range(3):
            #print "itype=%d iplab=%d nevt=%d" % (itype, iplab, nevt[itype][iplab])
            ana_args  = [];
            nevt_tot = 0;
            #print nfile_dict[itype][iplab]
            #print sorted(nfile_dict[itype][iplab].keys())
            for idx in sorted(nfile_dict[itype][iplab].keys()):
                #print "idx=%d nevt_tot=%d nevt_this=%d" % (idx,nevt_tot,nfile_dict[itype][iplab][idx])
                if (nevt_tot + nfile_dict[itype][iplab][idx]) < nevt[itype][iplab]:
                    ana_args.append((iplab,itype,ibrem,idx,nfile_dict[itype][iplab][idx]))
                    nevt_tot = nevt_tot + nfile_dict[itype][iplab][idx]
                else:
                    ana_args.append((iplab,itype,ibrem,idx,nevt[itype][iplab]-nevt_tot))
                    break

            print ana_args
            nevt_tot_check = 0
            for arg in ana_args:
                nevt_tot_check += arg[4]

            if (nevt_tot_check != nevt[itype][iplab]):
                print "ERROR: nevt_tot_check = %d != nevt_tot_expected= %d" % ( nevt_tot_check, nevt[itype][iplab])
            else:
                print "OK: nev tot checks out, launching ana on %d files!" % nfile(itype, iplab)
            print "========================"
            for args in ana_args:
                run_ana(args)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "hadd":
            print "calling hadd"
            hadd()
        else:
            print "Argument %s not supported" % sys.argv[1]
    else:
        run_parallel()

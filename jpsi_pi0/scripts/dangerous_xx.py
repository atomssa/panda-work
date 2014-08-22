#!/usr/bin/python

import subprocess,glob
from time import sleep

def du(idx):
    output = "*/*_%d.*" % i
    #if os.path.exists(output):
    count = glob.glob(output)
    #print count
    #if len(count) == 0:
    #    print "%s doesn't exist" % output
    #el
    if len(count) == 5:
        subprocess.Popen("echo -n \"%s \"; wc -l %s; echo" % (output,output), shell=True)
        subprocess.Popen("rm -vf %s" % output, shell=True)
    else:
        print "too many of %s -> %d " % (output,len(count))
    sleep(0.1)
        
s1 = range(1000, 1100, 3)
s2 = range(1001, 1100, 3)
#s3 = range(1002, 1100, 3)

print "s1 = %s" % s1
print "s2 = %s" % s2
#print "s3 = %s" % s3

for i in s1:
    du(i)

for i in s2:
    du(i)    

#for i in s3:
#    du(i)

#print "\n"

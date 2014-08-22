#!/usr/bin/python

import sys,os,glob

idx = int(sys.argv[1])

pattern = "*/*_%i.*"%idx

print "glob patern: %s"%pattern

files = glob.glob(pattern)

print ( "you requested to remove these files:")

for cand in files:
    print (cand)

if (len(files)==0):
    print ("nothing that matches %s was found"%pattern)
    sys.exit(1)
    
ans = raw_input("do you really want to remove these files?")
if ans == "yes_i_do":
    print ("you replied yes removing...")
    for cand in files:
        print os.remove(cand)
else:
   print ("type yes_i_do to confirm")

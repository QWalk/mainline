#!/usr/bin/python
""" Help help """

import sys
import getopt
import random

try:
    opts, args = getopt.getopt(sys.argv[1:], "hn:o:", ["help"])
except getopt.error, msg:
    print msg
    print "for help use --help"
    sys.exit(2)

outbase=""
np=0
for o, a in opts:
    if o in ("-h", "--help"):
	print __doc__
	sys.exit(0)
    if o in ("-o"):
	outbase=a
    if o in ("-n"):
	np=int(a)

if(outbase=="" or np==0):
    print __doc__
    sys.exit(0)

config_string="SAMPLE_POINT"
all_configs=[]
for file in args:
    f=open(file,"r")
    fstr=f.read()
    f.close()
    configs=fstr.split(config_string)
    #leave out the RANDNUM stuff
    all_configs.extend(configs[1:])

print "total number of configs ",len(all_configs)
print outbase,np

#we should generate new walkers rather than erroring out, by doing the 
#following in the different methods:
#VMC: just copy a few walkers at random; they'll decorrelate quickly
# this could cause problems for an optimize run right after, should
# perhaps warn about that..
#RMC: copy enough walkers, but warn that a block should be done to decorrelate
#DMC: copy walkers and either warn or just divide the weights by two
# of the walkers that we did copy, which is exact.
if(len(all_configs)%np!=0):
    print "Can't repartion for the moment, since the number of processors \
doesn't divide the number of configurations"
    sys.exit(1)



#now we can assume that len(all_configs)%np==0
nconfig=len(all_configs)/np
counter=0
for i in range(0,np):
    name=outbase+"_"+str(i)
    f=open(name,"w")
    is1=random.randint(1,10000000)
    is2=random.randint(1,10000000)
    f.write("RANDNUM  " + str(is1) + "  " + str(is2) + "  \n")
    for j in range(0,nconfig):
	f.write(config_string)
	f.write(all_configs[counter])
	counter+=1
    

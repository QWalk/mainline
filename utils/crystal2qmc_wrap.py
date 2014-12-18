#!/usr/bin/python
## A PYTHON script to convert CRYSTAL output to Qwalk input, 
##    based on Qwalk crystal2qmc 
## 
## Author:
##      Huihuo Zheng, 
##      Department of Physics, UIUC
##      hzheng8@illinois.edu / zhenghh04@gmail.com
##      Date: 02-02-2013
#############################################################################
from numpy import *
import os, sys
upf=1
if len(sys.argv)<2:
    print "Usage: crystal2qmc_wrap.py CRYSTALOUTPUT"
    exit()
f = open(sys.argv[1], "r")
print "read crystal output from file "+sys.argv[1]
def readstr(string, f):
    if string in f.read():
        f.seek(0)
        line=f.readline()
        while line.find(string)==-1:
            line=f.readline()
        return line
    else:
        return "NOT FOUND"
print "Finding number of k-points ..."
line = readstr("NUMBER OF K POINTS IN THE IBZ", f)
if (line !="NOT FOUND"):
    item = line.split()
    nkpts = int(item[len(item)-1])
    print "NUMBER OF K POINTS: ", nkpts
else:
    print "I don't know the number of k points"
    exit()
nline=(nkpts+3)/4
line = readstr("K POINTS COORDINATES", f)
kr=[]
ki=[]
m=0
n=0
if (line !="NOT FOUND"):
    for i in range(nline):
        items = f.readline().split()
        m = 0
        while ( n < nkpts and m < len(items)/4):
            if items[4*m].find("R")!=-1:
                kr.append(4*i+m)
            else:
                ki.append(4*i+m)
            n = n+1
            m = m+1
line = readstr("WEIGHT OF K POINTS - MONKHORST", f)
if (line != "NOT FOUND"):
    n = 0
    kweights = []
    kw = open("kweights.dat", "w")
    while ( n < nkpts):
        items = f.readline().split()
        for d in items:
            kweights.append(float(d))
            n +=1
            kw.write("%s\n" %float(d))
    kw.close()
else:
    print "I don't find the weights for k"
if len(sys.argv)>2 and sys.argv[2]=="weights":
    print "Generate k-points weights only, to kweights.dat\nnow exiting..."
    exit()
else:
    print "Generating qwalk input..."

print "real k-point index", kr
print "complex k-point index", ki

#line = readstr("SUMMED SPIN DENSITY", f)
#if (line !="NOT FOUND"):
#    items = line.split()
#    spin = float(items[len(items)-1])
#    print "TOTAL SPIN DENSITY: ", spin
'''
    if ( spin < 0.1 ):
        print "Treated as singlet states"
        upf = 1
    else:
        print "Treated as doublet..."
        upf = 2
'''
#    upf = int(spin) + 1
j=0
for i in kr:
    fout=open("crystal2qmc.inp", "w")
#    fout.write("%s" %upf+"\n"+"%s" %j)
    fout.write("%s" %j)
    fout.close()
    os.system(". ~/.bashrc; runcrystal2qmc -o qwalk_%s %s < crystal2qmc.inp" %(i, sys.argv[1]))
    os.system("cp qwalk_%s.jast2 qwalk.jast2; rm qwalk_%s.jast2" %(i, i))
    j +=1

j=0    
for i in ki:
    fout=open("crystal2qmc.inp", "w")
#    fout.write("%s" %upf+"\n"+"%s" %j)
    fout.close()
    os.system(". ~/.bashrc; runcrystal2qmc -c -o qwalk_%s %s < crystal2qmc.inp" %(i, sys.argv[1]))
    j +=1
    os.system("cp qwalk_%s.jast2 qwalk.jast2; rm qwalk_%s.jast2" %(i, i))
exit()

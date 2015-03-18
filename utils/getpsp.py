#!/usr/bin/python
from numpy import * 
from utils import * 
import sys, os
fin = open(sys.argv[1], 'r')
alst=sys.argv[2].split()
for alabel in alst:
    read_to_str(fin, "INPUT")
    ind=fin.readline().split()

    Z = float(ind[0])
    
    lst = []
    i=1
    while (ind[i]!='0' and i<len(ind)):
        lst.append(int(ind[i]))
        i=i+1
    psp=[]
    for n in lst:
        tmp=[]
        for i in range(n):
            tmp.append(fin.readline().split())
        psp.append(tmp)
    naip = 6
    if len(psp) > 2:
        naip = 12
    else:
        naip = 6
    print "PSEUDO {\n  %s"%alabel + "\n  AIP %s"%naip +"\n  BASIS { %s " %alabel+"\n RGAUSSIAN\n  OLDQMC {"
    print "0.0 %s" %len(psp)
    str=""

    for i in range(1, len(psp)):
        str = str + "%s " %lst[i]
    print str+ "%s" %lst[0]

    for i in range(1, len(psp)):
        a=psp[i]
        for j in range(len(a)):
            if a[j][2]=='0':
                print "   2 %s %s" %(a[j][0], a[j][1])
            elif a[j][2] == "1":
                print "   3 %s %s" %(a[j][0], a[j][1])
            else:
                print "   1 %s %s" %(a[j][0], a[j][1])
    a=psp[0]
    for j in range(len(a)):
        if a[j][2]=='0':
            print "   2 %s %s" %(a[j][0], a[j][1])
        elif a[j][2] == "1":
            print "   3 %s %s" %(a[j][0], a[j][1])
        else:
            print "   1 %s %s" %(a[j][0], a[j][1])
    print "   }\n  }\n}"

    print "\n"

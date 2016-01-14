#!/usr/bin/python
## This script is for computing finite size correction for solid in the potential energy part, based on PRL 97, 076404 (2006).  
import numpy as np
import matplotlib.pylab as plt
from fitlab import *
import os
import json
import copy
pi = np.pi
fjson="record.json" ## json file where we can read sk and supercell information
nresult=0 ## the run
fcif="si.cif" ## provide this if you want me to get out the volume 
dim=3 ## dimension of the system
maxiter=2 ## number of gaussion for fitting sk, please check skfit.pdf whether it is a good fit
compute_volume=True ## input volume v by hand in next line
volume=1.0 ## if compute_volume is True, this line is ignored
### To the correction, please make sure first that we have a good continue fit for the computed sk (please check skfit.pdf)
### We also calculated the zeta(0)/Omega in Eq.(8) of PRL 97, 076404 (2006)
def sort_with_x(x, y):
    return zip(*sorted(zip(x, y)))

def sk(fjson,nresult=0):
    f=open(fjson).read()
    data=json.loads(f)
    SK = np.array(data['qmc']['dmc'][u'results'][nresult][u'sk'])
    if SK==[]:
        print("Error, SK data is not founded in json file. Please make sure that you did the run.")
        return None
    SK = np.loadtxt("sk.out")
    sk =  SK[np.argsort(SK[:, 0])]
    return sk

def getvolume(fjson, fcif):
    f=open(fjson).read()
    data=json.loads(f)
    f = os.system("grep _cell_volume %s | sed 's/_cell_volume//g' > v.out" %fcif)
    return np.loadtxt("v.out")*np.linalg.det(data['supercell'])
    
def calccontinue(SK, dk=0.001, dim=3, direction='z', maxiter=2):
    k=[]
    sk=[]
    ske=[]
    n=0
    d = len(SK)
    c = 1.0
    SK0=[]
    for m in range(1, d-1):
        if (SK[d-m-1][1] - c) < 0.1:
            SK0.append(list(SK[d-m-1]))
            c=SK[d-m-1][1]
    for a in SK0:
        if a[1] < 1.1:
            if dim==3:
                pre = 1.0/(pi**3*16)*4.0*pi
                k.append(a[0])
                sk.append(a[1])
                ske.append(a[2])
            if dim==2:
                pre = 1.0/(pi**2*8)*2.0*pi
                if ( direction=='x' and a[3]==0.0 ) or \
                       ( direction=='y' and a[4]==0.0 ) or \
                       ( direction=='z' and a[5]==0.0 ):
                    k.append(a[0])
                    sk.append(a[1])
                    ske.append(a[2])
    p=FitXY(k, sk, maxiter=maxiter)
    kmax=max(k)
    kk = np.arange(0.0, kmax, dk)
    yy = GaussFun(p, kk)
    plt.xlabel(r"k [A$^{-1}$]")
    plt.ylabel(r"S(k)")
    plt.ylim(-0.049, 1.19)
    plt.savefig("skfit.pdf")
    plt.close()
    if sorted(k)[0] == 0.0:
        kmin = sorted(k)[1]
    else:
        kmin = min(k)
#    kkk = np.arange(0.0, kmin, dk)
#    yyy = GaussFun(p, kkk)
    if dim==2:
        vd = 2.0*pi
    elif dim==3:
        vd = 4.0*pi
    return sum(yy - 1.0)*dk*pre*vd #, p, sum(yyy - 1.0)*dk*pre*vd

def polyfitsk(SK, deg = 4):
    ''' assuming the k is already sorted '''
    d = len(SK)
    c = 1.0
    SK0=[]
    for m in range(1, d-1):
        if (SK[d-m-1][1] - c) < 0.1:
            SK0.append(list(SK[d-m-1]))
            c=SK[d-m-1][1]
    SK0 = np.array(SK0)
    k = [0.0, ] + list(SK0.transpose()[0][::-1])
    sk = [0.0,] + list(SK0.transpose()[1][::-1])
    i = 1
    while ( sk[i] < 0.9 ) and i < len(k):
        i += 1
    fit = np.polyfit(k[:i], sk[:i], deg = deg)
    kk = np.arange(0.0, k[i], 0.01)
    fsk = np.poly1d(fit)
    skk = [ fsk(a) for a in kk ]
    plt.plot(k, sk, 'o-')
    plt.plot(kk, skk, '-')
    plt.savefig("polyfit.pdf")
    return fit


def calcdiscret(SK, V = 1.0, dim=3, direction='z'):
    s=0.0
    pre = 1.0/2.0/V
    for a in SK:
        if a[0]>0.0:
            if dim==3:
                s+=(a[1] - 1.0)/a[0]**2*pre*4.0*pi
            if dim==2:
                if ( direction=='x' and a[3]==0.0 ) or \
                        ( direction=='y' and a[4]==0.0 ) or \
                        ( direction=='z' and a[5]==0.0 ):
                    s+=(a[1] - 1.0)/a[0]*pre*2.0*pi
    return s

def fitdiscret(SK, p, V = 1.0, dim=3, direction='z'):
    s=0.0
    pre = 1.0/2.0/V
    sk = GaussFun(p, SK.transpose()[0])
    n = 0
    plt.plot(SK.transpose()[0], sk, 'o', color='k')
    plt.plot(SK.transpose()[0], SK.transpose()[1], 'v', color='r')
    for a in SK:
        if a[0]>0.0:
            if dim==3:
                s+=(sk[n] - 1.0)/a[0]**2*pre*4.0*pi
            if dim==2:
                if ( direction=='x' and a[3]==0.0 ) or \
                        ( direction=='y' and a[4]==0.0 ) or \
                        ( direction=='z' and a[5]==0.0 ):
                    s+=(sk[n] - 1.0)/a[0]*pre*2.0*pi
        n+=1
    return s


SK = np.array(sk("record.json"))
if compute_volume==True:
    v = getvolume("record.json", "si.cif")
vd = calcdiscret(SK, V = v, dim=dim)
vc = calccontinue(SK, dim=dim, maxiter=maxiter)
p = polyfitsk(SK, deg=4)
print("**Zeta(0)/Omega: %s [Ha]" %(p[2]/2.0/v))
#vdf = fitdiscret(SK, p, V = v, dim=dim)
print("**Integral: %s [Ha] "  %vc)
print("**Discrete Sum: %s [Ha] "  %vd)
print("**Finite size correction: %s [Ha] " %(vc - vd ))
#print("**The contribution from ignoring k=0 part: %s [Ha] " %)
plt.savefig("sk_dat.pdf")

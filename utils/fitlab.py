import matplotlib.pylab as plt
import numpy as np
from scipy.optimize import leastsq, curve_fit
ms = 3.0
def GaussFun(p, x):
    yy = np.zeros(len(x))
    for m in range(len(x)):
        a = x[m]
        c = 0.0
        for i in range(int(len(p)/2)):
            c += np.exp(-p[2*i+1]*a**2)*p[2*i]
        c = 1.0 - c
        yy[m] = c
    return yy
def Residue(p, x, y):
    z = x
    yy = (y - GaussFun(p, x))/(x+0.01)
    return yy

def Residue2(p, x, y):
    z = x
    yy = (y - GaussFun(p, x))/(x+0.01)**2
    return yy

def Gauss(x, p0, p1):
    return np.exp(-p0*x**2)*p1

def FitXY(x, y, yerr=[], maxiter=3, color='Red', fmt='o', resfun=Residue, addr=""):
    p0 = [1.0, 1.0]
    xmax = max(x)
    if x[0] !=0.0:
        x = np.array([0.0, ] + list(x))
        y = np.array([0.0, ] + list(y))
        yerr = np.array([0.0, ] + list(yerr))
    xx = np.arange(0.0, xmax, 0.01)
    for j in range(maxiter):
        p, n = leastsq(resfun, p0, args=(x, y), maxfev=10000)
        p0 = list(p)+list(p0)
    yy = GaussFun(p,  xx)
    if len(yerr)<2:
        plt.plot(x, y, fmt, color=color, markersize=ms)
    else:
        plt.errorbar(x, y, yerr=yerr, fmt=fmt, color=color, markersize=ms)
    plt.plot(xx, yy, color=color, markersize=ms)
    return p


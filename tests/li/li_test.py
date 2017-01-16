import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad,dblquad,tplquad
import scipy as sp

#Active space down check
#Gradient: correct
#Gradderiv: still have to check
#Check over time: still have to check

def phi1(r):
	a1=1.
	return np.sqrt(a1/np.pi)*np.exp(-a1*r**2)

def phi2(r):
	a2=2.
	return np.sqrt(a2/np.pi)*np.exp(-a2*r**2)

def phi3(r):
	a3=3.
	return np.sqrt(a3/np.pi)*np.exp(-a3*r**2)


'''#Determinant down
def ddown(p1, p2, r):
	vec=np.array([phi3(r), phi2(r), phi1(r)])
	return np.dot(Rdown(p1, p2),vec)[0]

#Down rotation matrix
def Rdown(p1, p2):
	theta=np.zeros((3,3))
	theta[0,1]=p1
	theta[0,2]=p2
	theta[1,0]=-p1
	theta[2,0]=-p2
	return expm(theta)

#Parameter derivative
def pderivdown(p1, p2, r, num):
	vec=np.array([phi3(r), phi2(r), phi1(r)])
	if(num==2):
		return np.dot(Rdown(p1, p2),vec)[1]
	elif(num==3):
		return np.dot(Rdown(p1, p2),vec)[2]

#Laplacian 
def lap(p1, p2, r, num):
	dr=0.001
	val=(2./r)*(pderivdown(p1, p2, r+dr, num)-pderivdown(p1, p2, r, num))/dr
	val+=(pderivdown(p1, p2, r+dr, num)+pderivdown(p1, p2, r-dr, num)-2.*pderivdown(p1, p2, r, num))/(dr**2.)
	return val

#Determinant up
def dup(p1, p2, r1, r2):
	vec1=np.array([phi2(r1), phi1(r1), phi3(r1)])
	vec2=np.array([phi2(r2), phi1(r2), phi3(r2)])
	one=np.dot(Rup(p1,p2),vec1)
	two=np.dot(Rup(p1,p2),vec2)

	return one[0]*two[1]-one[1]*two[0]
#Rotation matrix
def Rup(p1, p2):
	theta=np.zeros((3,3))
	theta[0,2]=p1
	theta[1,2]=p2
	theta[2,0]=-p1
	theta[2,1]=-p2

	return expm(theta)

#Laplacian 
def lapu(r1, r2, e):
	dr=0.001
	if(e==0):
		return (2./r1)*(dup(0,0,r1+dr,r2)-dup(0,0,r1,r2))/dr + (dup(0,0,r1+dr,r2)+dup(0,0,r1-dr,r2)-2.*dup(0,0,r1,r2))/dr**2
	elif(e==1):
		return (2./r2)*(dup(0,0,r1,r2+dr)-dup(0,0,r1,r2))/dr + (dup(0,0,r1,r2+dr)+dup(0,0,r1,r2-dr)-2.*dup(0,0,r1,r2))/dr**2

#Integration 
p1=0.8529006671
p2=7.76876134
e=1
for num in range(2,4):
	denom,err=quad(lambda x: ddown(p1,p2,x)**2*x**2,0,np.inf)
	numer,err=quad(lambda x: pderivdown(p1,p2,x,num)*ddown(p1, p2, x)*x**2,0,np.inf)
	#numer,err=quad(lambda x: lap(p1,p2,x,num)*ddown(p1, p2, x)*x**2,0,np.inf)
	numer2,err=dblquad(lambda x,y: lapu(x,y,e)*dup(0,0,x,y)*x**2*y**2,0,np.inf,lambda x:0, lambda x:np.inf)
	denom2,err=dblquad(lambda x,y: dup(0,0,x,y)**2*x**2*y**2,0,np.inf,lambda x:0, lambda x:np.inf)
	print (numer*numer2/(denom*denom2))'''
	#print numer/denom


#Active space up check
#Determinant up
def dup(p1, p2, r1, r2):
	vec1=np.array([phi1(r1), phi3(r1), phi2(r1)])
	vec2=np.array([phi1(r2), phi3(r2), phi2(r2)])
	one=np.dot(Rup(p1,p2),vec1)
	two=np.dot(Rup(p1,p2),vec2)

	return one[0]*two[1]-one[1]*two[0]

#Derivatives
def pderivup(p1, p2, r1, r2, num):
	vec1=np.array([phi1(r1), phi3(r1), phi2(r1)])
	vec2=np.array([phi1(r2), phi3(r2), phi2(r2)])
	if(num==0):
		one=np.dot(Rup(p1,p2),vec1)
		two=np.dot(Rup(p1,p2),vec2)
		return one[2]*two[1]-one[1]*two[2]

	elif(num==1):
		one=np.dot(Rup(p1,p2),vec1)
		two=np.dot(Rup(p1,p2),vec2)
		return one[0]*two[2]-one[2]*two[0]

#Rotation matrix
def Rup(p1, p2):
	theta=np.zeros((3,3))
	theta[0,2]=p1
	theta[1,2]=p2
	theta[2,0]=-p1
	theta[2,1]=-p2

	return expm(theta)

#Def Laplacian
def lap(p1, p2, r1, r2, num, e):
	dr=0.001
	if(e==1):
		val=(2./r1)*(pderivup(p1, p2, r1+dr,r2, num)-pderivup(p1, p2, r1,r2, num))/dr
		val+=(pderivup(p1, p2, r1+dr,r2, num)+pderivup(p1, p2, r1-dr,r2, num)-2.*pderivup(p1, p2, r1,r2, num))/(dr**2.)
		return val
	elif(e==2):
		val=(2./r2)*(pderivup(p1, p2, r1, r2+dr, num)-pderivup(p1, p2, r1, r2, num))/dr
		val+=(pderivup(p1, p2, r1, r2+dr, num)+pderivup(p1, p2, r1, r2-dr, num)-2.*pderivup(p1, p2, r1, r2, num))/(dr**2.)
		return val

def ddown(p1, p2, r):
	vec=np.array([phi2(r), phi1(r), phi3(r)])
	return np.dot(Rdown(p1, p2),vec)[0]

#Down rotation matrix
def Rdown(p1, p2):
	theta=np.zeros((3,3))
	theta[0,1]=p1
	theta[0,2]=p2
	theta[1,0]=-p1
	theta[2,0]=-p2
	return expm(theta)

#Lapd
def lapd(r):
	dr=0.001
	return (2./r)*(ddown(0,0,r+dr)-ddown(0,0,r))/dr + (ddown(0,0,r+dr)+ddown(0,0,r-dr)-2.*ddown(0,0,r))/dr**2


#Local energy times wave function
def el(p1, p2, x, y, z):
	dr=0.001
	one=-0.5*((2./x)*(dup(p1,p2,x+dr,y)-dup(p1,p2,x,y))/dr + (dup(p1,p2,x+dr,y)+dup(p1,p2,x-dr,y)-2.*dup(p1,p2,x,y))/dr**2)*ddown(0,0,z)
	two=-0.5*((2./y)*(dup(p1,p2,x,y+dr)-dup(p1,p2,x,y))/dr + (dup(p1,p2,x,y+dr)+dup(p1,p2,x,y-dr)-2.*dup(p1,p2,x,y))/dr**2)*ddown(0,0,z)
	three=-0.5*((2./z)*(ddown(0,0,z+dr)-ddown(0,0,z))/dr + (ddown(0,0,z+dr)+ddown(0,0,z-dr)-2.*ddown(0,0,z))/dr**2)*dup(p1,p2,x,y)

	if(x-y!=0 and y-z!=0 and x-z!=0):
		four=dup(p1,p2,x,y)*ddown(0,0,z)*(1./abs(x-y)+1./abs(y-z)+1./abs(x-z)-1./x-1./y-1./z)
	else:
		return 0

	return one+two+three+four


#Integration
p1=1.13
p2=7.95

'''for e in range(2,3):
	for num in range(0,2):
		print("num=",num)
		print("e=",e)
		#numer,err=dblquad(lambda y,x:lap(p1,p2,x,y,num,e)*dup(p1,p2,x,y)*x**2*y**2,0,np.inf,lambda y:0, lambda y:np.inf)
		numer,err=dblquad(lambda y,x:pderivup(p1,p2,x,y,num)*dup(p1,p2,x,y)*x**2*y**2,0,np.inf,lambda y:0, lambda y:np.inf)
		denom,err=dblquad(lambda y,x:dup(p1,p2,x,y)**2*x**2*y**2,0,np.inf,lambda y:0, lambda y:np.inf)
		
		numer2,err=quad(lambda x: lapd(x)*ddown(0,0,x)*x**2,0,np.inf)
		denom2,err=quad(lambda x: ddown(0,0,x)**2*x**2,0,np.inf)

		print (numer*numer2)/(denom*denom2)
		#print numer/denom'''
'''
print("yea")
numer,err=dblquad(lambda y,x:pderivup(p1,p2,x,y,1)*pderivup(p1,p2,x,y,1)*x**2*y**2,0,np.inf,lambda y:0, lambda y:np.inf)
denom,err=dblquad(lambda y,x:dup(p1,p2,x,y)**2*x**2*y**2,0,np.inf,lambda y:0, lambda y:np.inf)
print(numer/denom)

'''
#lambda z,y,x:el(p1,p2,x,y,z)*ddown(p1,p2,x,y)*dup(0,0,z)*x**2*y**2*z**2
numer,err=tplquad(lambda z,y,x:el(p1,p2,x,y,z)*ddown(0,0,z)*dup(p1,p2,x,y)*x**2*y**2*z**2,1e-5,np.inf,lambda y:1e-5, lambda y:np.inf,lambda z,t:1e-5, lambda z,t:np.inf)
denom,err=tplquad(lambda z,y,x:dup(p1,p2,x,y)**2*ddown(0,0,z)**2*x**2*y**2*z**2,1e-5,np.inf,lambda y:1e-5, lambda y:np.inf,lambda z,t:1e-5, lambda z,t:np.inf)
print numer/denom






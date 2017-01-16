import numpy as np 
from scipy.linalg import eig

nparms=2
f=open('test.out','r')
data=np.zeros(2*(nparms+1)**2)
i=0
for line in f:
	data[i]=float(line)
	i+=1

H=np.zeros((nparms+1,nparms+1))
S=np.zeros((nparms+1,nparms+1))

for i in range(nparms+1):
	for j in range(nparms+1):
		H[i][j]=data[(nparms+1)*i+j]
		S[i][j]=data[(nparms+1)**2+(nparms+1)*i+j]

w,vr=eig(H,S)
print(min(w))
i=np.argmin(w)
print vr[:,i]
import numpy as np 
import matplotlib.pyplot as plt 

E=[]
Err=[]
data=np.zeros(10000)
f=open('en.out','r')
isErr=0
k=0
for line in f:
	data[k]=float(line)
	k+=1

isErr=0
for i in range(len(data)):
	en=data[i]
	if(isErr==1):
		if(en in Err):
			isErr=0
		else:
			Err.append(en)
			isErr=0
	if(en<0):
		if(en in E):
			pass
		else:
			E.append(en)
			isErr=1

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="New, exact")
plt.title("Li Test")
plt.ylabel("Energy, Ha")
plt.xlabel("Iteration")
plt.legend(loc=1)
plt.show()

import numpy as np 
import matplotlib.pyplot as plt 
E=[]
Err=[]
data=np.zeros(10000)
f=open('en2.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="New solver, exact derivatives w/ Jastrow")

E=[]
Err=[]
data=np.zeros(10000)
f=open('en3.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="New solver, exact derivatives w/ Jastrow")

E=[]
Err=[]
data=np.zeros(10000)
f=open('en1.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="New solver, exact derivatives w/out Jastrow")

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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test")

'''
E=[]
Err=[]
data=np.zeros(10000)
f=open('en4.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 4")
'''
plt.title(r"$H_2$ Test")
#plt.errorbar(np.arange(len(E)),E,yerr=np.ones(len(E))*0.01,label="Optimization, 2000 VMC steps, S-J 0 parameters Frozen")
plt.ylabel("Energy (Ha)")
plt.xlabel("Iteration")
plt.legend(loc=1)
plt.show()

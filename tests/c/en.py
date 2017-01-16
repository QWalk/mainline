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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 2")

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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 3")


E=[]
Err=[]
data=np.zeros(10000)
f=open('en8.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 8")

E=[]
Err=[]
data=np.zeros(10000)
f=open('en9.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 9")
'''
E=[]
Err=[]
data=np.zeros(10000)
f=open('en5.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 5")

E=[]
Err=[]
data=np.zeros(10000)
f=open('en6.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 6")

E=[]
Err=[]
data=np.zeros(10000)
f=open('en7.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 7")

E=[]
Err=[]
data=np.zeros(10000)
f=open('en8.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 8")
plt.plot(np.arange(len(E)), -5.37439033*np.ones(len(E)),'--',label="Test 7,8 VMC")

E=[]
Err=[]
data=np.zeros(10000)
f=open('en9.out','r')
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

plt.errorbar(np.arange(len(E)),E,yerr=Err,label="Test 9")
plt.plot(np.arange(len(E)),  -4.142107278*np.ones(len(E)),'--',label="Test 9 VMC")
'''
plt.title("2-Determinant C Test")
plt.ylabel("Energy (Ha)")
plt.xlabel("Iteration")
plt.legend(loc=1)
plt.show()

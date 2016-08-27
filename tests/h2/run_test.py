#!/usr/bin/env python
from __future__ import print_function
import subprocess
import json
import sys
import math
import os

sys.path.append("../")
from qwtest import *

reports=[]

print("""###########################################
Checking the linear method for the H2 molecule. 
This tests the LINEAR method, the Jastrow factor, analytic parameter derivatives..
The reference is a previous QWalk run. 
################################################""")

try:
  subprocess.check_output(['rm','qw.linear.o','qw.linear.config'])
except:
  pass
subprocess.check_output([QW,'qw.linear'])

report={'system':'h2', 'method':'linear'}
success={}

ref_data=[-1.068309651,0.0009113787988]

run_data=[0.0,0.0]
with open("qw.linear.o") as f:
  for line in f:
    if 'current energy' in line:
      spl=line.split()
      run_data[0]=float(spl[4])
      run_data[1]=float(spl[6])

report['description']='Checking the linear method for the H2 molecule. This tests the LINEAR method, the Jastrow factor, analytic parameter derivatives. The reference is a previous QWalk run.'
report['quantity']='total_energy'
report['result']=run_data[0]
report['error']=run_data[1]
report['reference']=ref_data[0]
report['err_ref']=ref_data[1]

#I added a factor of 2 fudge factor to account for differences in optimization.
success['total_energy']=check_errorbars(ref_data[0],run_data[0],2.*math.sqrt(ref_data[1]**2+run_data[1]**2))
report['passed']=success['total_energy']

reports.append(report)

print_results(reports)
save_results(reports)

if success['total_energy']:
  sys.exit(0)
sys.exit(1)

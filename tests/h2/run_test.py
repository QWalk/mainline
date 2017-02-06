#!/usr/bin/env python
from __future__ import print_function
import subprocess
import json
import sys
import os

sys.path.append("../")
from qwtest import *

reports=[]
allsuc=[]

print(""" Testing derivatives """)


test_str=str(subprocess.check_output([QW,'qw.testder']))

if test_str.count("FAILED") > 4:
  print("Failed")
else:
  print("Succeeded")



print("""###########################################
Checking the linear method for the H2 molecule. 
This tests the LINEAR method, the Jastrow factor, analytic parameter derivatives..
The reference is a previous QWalk run. 
################################################""")

ref_data={'total_energy':[-1.068309651,0.0009113787988]}
systems={'total_energy':'h2'}
methods={'total_energy':'linear'}
descriptions={'total_energy':'Checking the linear method for the H2 molecule. This tests the LINEAR method, the Jastrow factor, analytic parameter derivatives. The reference is a previous QWalk run.'}
#I added a factor of 2 fudge factor to account for differences in optimization.
sigmas={'total_energy':6.0}

try:
  subprocess.check_output(['rm','qw.linear.o','qw.linear.config'])
except:
  pass

subprocess.check_output([QW,'qw.linear'])

run_data=[0.0,0.0]
with open("qw.linear.o") as f:
  for line in f:
    if 'current energy' in line:
      spl=line.split()
      run_data[0]=float(spl[4])
      run_data[1]=float(spl[6])
dat_properties={'total_energy':{'value':[run_data[0]], 'error':[run_data[1]]}}

success=compare_result_ref(ref_data,dat_properties,sigmas)
for k,v in success.items():
  allsuc.append(v)

reports.extend(summarize_results(ref_data,dat_properties,success,systems,methods,descriptions))

print_results(reports)
save_results(reports)

if False in allsuc:
  sys.exit(1)
sys.exit(0)

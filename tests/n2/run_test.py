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

print("""###########################################
Checking VMC and Slater determinant for the stretched N2 dimer.
This tests the orbital evaluation routine, the VMC method, ECPs, and multiple determinants.
The reference is the CI result from GAMESS.
################################################""")

ref_data={'total_energy': [-19.3080264611,0.0],
          'kinetic':[12.2497549662,0.0],
          }
systems={}
methods={}
descriptions={}
sigmas={}
for k in ref_data.keys():
  systems[k]='n2'
  methods[k]='vmcci'
  descriptions[k]='Checking VMC and Slater determinant for the stretched N2 dimer.  This tests the orbital evaluation routine, the VMC method, ECPs, and multiple determinants.  The reference is the CI result from GAMESS.'
  sigmas[k]=3.0

try:
  subprocess.check_output(['rm','qw.ci.log','qw.ci.config'])
except:
  pass

subprocess.check_output([QW,'qw.ci'])

json_str=subprocess.check_output([GOS,'-json','qw.ci.log'])
dat=json.loads(json_str)
dat_properties=dat['properties']

success=compare_result_ref(ref_data,dat_properties,sigmas)
for k,v in success.items():
  allsuc.append(v)

reports.extend(summarize_results(ref_data,dat_properties,success,systems,methods,descriptions))

print_results(reports)
save_results(reports)

if False in allsuc:
  sys.exit(1)
sys.exit(0)

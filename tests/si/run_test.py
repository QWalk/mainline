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
Silicon Slater determinant
################################################""")

ref_data={}
ref_data['111']={'total_energy': [-7.538,0.02],
          'kinetic':[3.38,0.02],
          }
ref_data['222']={'total_energy': [-7.71,0.02],
          'kinetic':[3.12,0.02],
          }

for bc in ['111','222']:
  for orb in ['blas','cutoff']:
    fname='qw_'+bc+'.'+orb+'.hf'
    systems={}
    methods={}
    descriptions={}
    sigmas={}
    for k in ref_data[bc].keys():
      systems[k]='silicon'
      methods[k]='vmc'
      descriptions[k]='PBCs and k-points'
      sigmas[k]=3.0

    try:
      subprocess.check_output(['rm',fname+'.log'])
    except:
      pass

    subprocess.check_output([QW,fname])
    json_str=subprocess.check_output([GOS,'-json',fname+'.log'])
    dat=json.loads(json_str)
    dat_properties=dat['properties']

    success=compare_result_ref(ref_data[bc],dat_properties,sigmas)
    for k,v in success.items():
      allsuc.append(v)

    reports.extend(summarize_results(ref_data[bc],dat_properties,success,systems,methods,descriptions))

print_results(reports)
save_results(reports)

if False in allsuc:
  sys.exit(1)
sys.exit(0)

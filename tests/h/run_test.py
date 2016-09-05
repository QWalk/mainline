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
Checking VMC and Slater determinant for the H atom. 
This tests the orbital evaluation routine and the VMC method.
The reference is the Hartree-Fock result from GAMESS.
################################################""")

ref_data={'total_energy':[-0.4998098113,0.0],
          'kinetic':[0.4997883996,0.0],
          'potential':[-0.9995982109,0.0]
          }
systems={}
methods={}
descriptions={}
sigmas={}
for k in ref_data.keys():
  systems[k]='h'
  methods[k]='vmc'
  descriptions[k]='Checking VMC and Slater determinant for the H atom. This tests the orbital evaluation routine and the VMC method. The reference is the Hartree-Fock result from GAMESS.'
  sigmas[k]=3.0

try:
  subprocess.check_output(['rm','qw.hf.log','qw.hf.config'])
except:
  pass

subprocess.check_output([QW,'qw.hf'])

json_str=subprocess.check_output([GOS,'-json','qw.hf.log'])
dat=json.loads(json_str)
dat_properties=dat['properties']

success=compare_result_ref(ref_data,dat_properties,sigmas)
for k,v in success.items():
  allsuc.append(v)

reports.extend(summarize_results(ref_data,dat_properties,success,systems,methods,descriptions))

print("""###########################################
Checking DMC for the H atom. 
This tests the basic DMC algorithm.
The reference is the exact result.
################################################""")

ref_data={'total_energy':[-0.5,0.0] }
systems={'total_energy':'h'}
methods={'total_energy':'dmc'}
descriptions={'Checking DMC for the H atom.  This tests the basic DMC algorithm.  The reference is the exact result.'}
sigmas={'total_energy':3.0}


try:
  subprocess.check_output(['rm','qw.dmc.log','qw.dmc.config'])
except:
  pass

subprocess.check_output([QW,'qw.dmc'])

json_str=subprocess.check_output([GOS,'-json','qw.dmc.log'])
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

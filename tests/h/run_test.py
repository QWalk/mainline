#!/usr/bin/env python
from __future__ import print_function
import subprocess
import json
import sys
import csv
import os

sys.path.append("../")
from qwtest import *
QW="../../bin/qwalk"
GOS="../../bin/gosling"

def check_errorbars(a,b,berr,sigma=3):
  if abs(a-b)/berr > sigma:
    return False
  return True

def check_sane(a,b,berr,sigma=3):
  if berr > 0.1:
    return False
  return True


reports=[]

print("""###########################################
Checking VMC and Slater determinant for the H atom. 
This tests the orbital evaluation routine and the VMC method.
The reference is the Hartree-Fock result from GAMESS.
################################################""")

try:
  subprocess.check_output(['rm','qw.hf.log','qw.hf.config'])
except:
  pass
subprocess.check_output([QW,'qw.hf'])
json_str=subprocess.check_output([GOS,'-json','qw.hf.log'])

dat=json.loads(json_str)

ref_data={'total_energy':-0.4998098113,
          'kinetic':0.4997883996,
          'potential':-0.9995982109
          }


success={}
for k in ref_data.keys():
  success[k]=check_errorbars(ref_data[k],
                  dat['properties'][k]['value'][0],
                  dat['properties'][k]['error'][0])
  success[k+'sane']=check_sane(ref_data[k],
                  dat['properties'][k]['value'][0],
                  dat['properties'][k]['error'][0])
  report={'system':'h', 'method':'vmc'}
  report['quantity']=k
  report['description']='Checking VMC and Slater determinant for the H atom. This tests the orbital evaluation routine and the VMC method. The reference is the Hartree-Fock result from GAMESS.'
  report['result']=dat['properties'][k]['value'][0]
  report['error']=dat['properties'][k]['error'][0]
  report['reference']=ref_data[k]
  report['err_ref']=0.0
  if success[k] and success[k+'sane']:
    report['passed']=True
  else:
    report['passed']=False
  reports.append(report)

allsuc=[]
for k,v in success.items():
  allsuc.append(v)
  print(k,v)

print("""###########################################
Checking DMC for the H atom. 
This tests the basic DMC algorithm.
The reference is the exact result.
################################################""")
try:
  subprocess.check_output(['rm','qw.dmc.log','qw.dmc.config'])
except:
  pass
subprocess.check_output([QW,'qw.dmc'])

success={}

json_str=subprocess.check_output([GOS,'-json','qw.dmc.log'])
dat=json.loads(json_str)

report={'system':'h', 'method':'dmc'}
report['description']='Checking DMC for the H atom. This tests the basic DMC algorithm. The reference is the exact result.'
report['quantity']='total_energy'
report['result']=dat['properties']['total_energy']['value'][0]
report['error']=dat['properties']['total_energy']['error'][0]
report['reference']=-0.5
report['err_ref']=0.0

success['total_energy']=check_errorbars(-0.5,dat['properties']['total_energy']['value'][0],dat['properties']['total_energy']['error'][0])
report['passed']=success['total_energy']

reports.append(report)

print("Exact: -0.5")
print("Data: ",dat['properties']['total_energy']['value'][0],"+/-",dat['properties']['total_energy']['error'][0])

fieldnames = ['method','quantity','system','description','passed','result','error','reference','err_ref']
with open('report.csv','w') as csvfile:
  writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
  writer.writerow(dict(zip(fieldnames,fieldnames)))
  writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
  writer.writerows(reports)

if False in allsuc:
  exit(1)

if not success['total_energy']:
  print("FAILED")
  exit(1)
print("PASSED")
exit(0)

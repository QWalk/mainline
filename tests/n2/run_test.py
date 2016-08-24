#!/usr/bin/env python
from __future__ import print_function
import subprocess
import json
import sys
import csv
import os

sys.path.append("../")
from qwtest import *

reports=[]

print("""###########################################
Checking VMC and Slater determinant for the stretched N2 dimer.
This tests the orbital evaluation routine, the VMC method, ECPs, and multiple determinants.
The reference is the CI result from GAMESS.
################################################""")

try:
  subprocess.check_output(['rm','qw.ci.log','qw.ci.config'])
except:
  pass
subprocess.check_output([QW,'qw.ci'])
json_str=subprocess.check_output([GOS,'-json','qw.ci.log'])

dat=json.loads(json_str)

ref_data={'total_energy': -19.3080264611,
          'kinetic':12.2497549662,
          }


success={}
for k in ref_data.keys():
  success[k]=check_errorbars(ref_data[k],
                  dat['properties'][k]['value'][0],
                  dat['properties'][k]['error'][0])
  success[k+'sane']=check_sane(ref_data[k],
                  dat['properties'][k]['value'][0],
                  dat['properties'][k]['error'][0])
  report={'system':'n2', 'method':'vmcci'}
  report['quantity']=k
  report['description']='Checking VMC and Slater determinant for the stretched N2 dimer.  This tests the orbital evaluation routine, the VMC method, ECPs, and multiple determinants.  The reference is the CI result from GAMESS.'
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

fieldnames = ['method','quantity','system','description','passed','result','error','reference','err_ref']
if not os.path.exists('../report.csv'):
  with open('../report.csv','w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writerow(dict(zip(fieldnames,fieldnames)))
with open('../report.csv','a') as csvfile:
  writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
  writer.writerows(reports)

if False in allsuc:
  exit(1)


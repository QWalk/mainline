from __future__ import print_function
import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control as jc
import os
import json

import taub

# Now we define the job sequence.
element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:05:00",queue="test")))
element_list.append(runcrystal.RunProperties(
  submitter=taub.LocalTaubPropertiesSubmitter(
    nn=1,np=1,time="0:05:00",queue="test")))
#element_list.append(runcrystal.RunProperties())
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(
  submitter=taub.LocalTaubQwalkSubmitter(
    nn=1,time="0:05:00",queue="test")))
element_list.append(runqwalk.QWalkRunDMC(
  submitter=taub.LocalTaubBundleQwalkSubmitter(
    nn=2,time="0:30:00",queue="secondary")))

qchecker = taub.LocalTaubSubmitter()

# Specific defaults for this material.
default_job=jc.default_job_record("si.cif")
default_job['dft']['kmesh'] = [2,2,2]
default_job['dft']['functional']['hybrid'] = 0
default_job['dft']['tolinteg'] = [6,6,6,6,12]
default_job['dft']['maxcycle'] = 100
default_job['dft']['fmixing'] = 80
default_job['dft']['edifftol'] = 6
default_job['dft']['broyden'] = [0.1,60,20]
default_job['qmc']['dmc']['save_trace'] = True
default_job['total_spin'] = 0
idbase = "si_ag_"

# A demonstration of varying basis parameters.
count=1
results = []

# Full runs.
for a in [0.2,0.3,0.4]:
  name = idbase+str(count)
  job_record = copy.deepcopy(default_job)
  job_record['dft']['basis']=[a,3,3]
  job_record['control']['id']=name
  results.append(jc.execute(job_record,element_list))
  count+=1

# Reduce allowed run time.
element_list[1] = runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:00:10",queue="test"))

# Demonstrate killed-job error correction.
for a in [0.2,0.3,0.4]:
  name = idbase+str(count)
  job_record = copy.deepcopy(default_job)
  job_record['dft']['basis']=[a,3,3]
  job_record['control']['id']=name
  job_record['dft']['maxcycle'] = 100
  check_continue = jc.check_continue(name,qchecker)
  if check_continue == "continue":
    restart_file = "/home/brian/projects/si/autogen/"+jc.continue_job(name)
  results.append(jc.execute(job_record,element_list))
  count+=1

# Save the data, either in JSON:
json.dump(results,open("data.json",'w'))
# or in YAML (which is easier to read for some).
#yaml.dump(results,open("data.yaml",'w'),default_flow_style=False)

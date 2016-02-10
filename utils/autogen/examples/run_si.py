from __future__ import print_function
import sys
sys.path.append("../")
sys.path.append("/projects/wagner/apps/qwalk_src/utils/autogen")
import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control as jc
import os
import json
import yaml

import taub

# Now we define the job sequence.
element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:20:00",queue="secondary")))
element_list.append(runcrystal.RunProperties(
  submitter=taub.LocalTaubPropertiesSubmitter(
    nn=1,np=1,time="0:20:00",queue="secondary")))
#element_list.append(runcrystal.RunProperties())
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(
  submitter=taub.LocalTaubQwalkSubmitter(
    nn=1,time="0:20:00",queue="secondary")))
element_list.append(runqwalk.QWalkRunDMC(
  submitter=taub.LocalTaubBundleQwalkSubmitter(
    nn=2,time="1:00:00",queue="secondary")))

qchecker = taub.LocalTaubSubmitter()

# Specific defaults for this material.
default_job=jc.default_job_record("si.cif")
idbase = "si_ag_"

# A demonstration of varying basis parameters.
count=1
results = []
for a in [0.2,0.3,0.4]:
  name = idbase+str(count)
  job_record = copy.deepcopy(default_job)
  job_record['dft']['basis']=[a,3,3]
  job_record['control']['id']=name
  results.append(jc.execute(job_record,element_list))
  count+=1


# Save the data, either in JSON:
json.dump(results,open("data.json",'w'))
# or in YAML (which is easier to read for some).
yaml.dump(results,open("data.yaml",'w'),default_flow_style=False)

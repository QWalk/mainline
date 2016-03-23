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


cifs = os.listdir(os.getcwd()+'/Cifs/')
results = []
for cif in cifs:
  if cif.split('.')[1] == 'cif':
    job_record=jc.default_job_record(os.getcwd()+'/Cifs/'+cif)
    job_record['control']['id'] = cif.split('.')[0]
    results.append(jc.execute(copy.deepcopy(job_record),element_list))


# Save the data, either in JSON:
json.dump(results,open("data.json",'w'))
# or in YAML (which is easier to read for some).
yaml.dump(results,open("data.yaml",'w'),default_flow_style=False)

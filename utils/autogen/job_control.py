from __future__ import print_function
import os
import json
import shutil
def default_job_record(ciffile):
  job_record={}
  job_record['dft']={}
  job_record['qmc']={}
  job_record['control']={}

  #Set up Hamiltonian
  with open (ciffile, "r") as f:
    job_record['cif']=f.read()
  job_record['supercell']=[[1,0,0],[0,1,0],[0,0,1]]
  job_record['pseudopotential']='BFD'
  job_record['charge']=0
  job_record['total_spin']=0

  #DFT-specific options
  job_record['dft']['functional']={'exchange':'PBE','correlation':'PBE','hybrid':25}
  job_record['dft']['basis']=[0.2,3,3]
  job_record['dft']['kmesh']=[8,8,8]
  job_record['dft']['spin_polarized']=True
  job_record['dft']['initial_spin']=[]
  job_record['dft']['initial_charges']={} #For example, 'O':-2,'Mg':2 
  job_record['dft']['fmixing']=99
  job_record['dft']['broyden']=[0.1,60,8]


  #QMC-specific options
  job_record['qmc']['timestep']=0.02
  job_record['qmc']['jastrow']='twobody'
  job_record['qmc']['optimize']='variance'
  job_record['qmc']['localization']='tmoves'
  job_record['qmc']['target_error']=0.01

  #Control options
  job_record['control']['id']=1
  job_record['control']['elements']=[]
  job_record['control']['pretty_formula']=''
  job_record['control']['queue_id']=[]
  return job_record

def execute(job_list, element_list):
  updated_job_list=[]
  for record in job_list:
    currwd=os.getcwd()
    d=str(record['control']['id'])
    try:
      os.mkdir(d)
    except:
      pass
    os.chdir(d)
    jsonfile="record.json"
    if os.path.isfile(jsonfile):
      f=open('record.json','r')
      record=json.load(f)
      f.close()


    print("#######################ID",record['control']['id'])

    for element in element_list:
      status=element.check_status(record)
      print(element._name_,"status",status)
      if status=='not_started':
        status=element.run(record)
        print(element._name_,"status",status)
      if status=='not_finished':
        status=element.retry(record)
        print(element._name_,"status",status)

      if status!='ok':
        break
      record=element.output(record)

    f=open(jsonfile,'w')
    json.dump(record,f)
    os.chdir(currwd)
    updated_job_list.append(record)
  return updated_job_list


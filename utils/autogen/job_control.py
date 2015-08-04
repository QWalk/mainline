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
  job_record['dft']['tolinteg']=[12,12,12,12,20]
  job_record['dft']['spin_polarized']=True
  job_record['dft']['initial_spin']=[]
  job_record['dft']['initial_charges']={} #For example, 'O':-2,'Mg':2 
  job_record['dft']['edifftol']=10
  job_record['dft']['fmixing']=99
  job_record['dft']['broyden']=[0.01,60,8]
  job_record['dft']['maxcycle']=200
  # None = fresh run, else copy this path to fort.20;
  # e.g. job_record['dft']['restart_from'] = ../successful_run/fort.9
  job_record['dft']['restart_from']=None

  #QMC-specific options
  job_record['qmc']['dmc']={}
  job_record['qmc']['dmc']['timestep']=[0.02]
  job_record['qmc']['dmc']['jastrow']='twobody'
  job_record['qmc']['dmc']['nblock']=16
  job_record['qmc']['dmc']['optimizer']='variance' #or energy
  job_record['qmc']['dmc']['localization']=['tmoves']
  job_record['qmc']['dmc']['target_error']=0.01
  job_record['qmc']['dmc']['kpoints']='real' # or a list of numbers
  job_record['qmc']['dmc']['excitations']='no'#VBM-CBM or other..


  job_record['qmc']['variance_optimize']={}
  job_record['qmc']['variance_optimize']['niterations']=10
  job_record['qmc']['variance_optimize']['nruns']=2


  job_record['qmc']['energy_optimize']={}
  job_record['qmc']['energy_optimize']['threshold']=0.001
  job_record['qmc']['energy_optimize']['vmc_nstep']=1000


  #Control options
  job_record['control']['id']=1
  job_record['control']['elements']=[]
  job_record['control']['pretty_formula']=''
  job_record['control']['queue_id']=[]
  job_record['control']['incomplete'] = False
  # Currently, force_retry clashes with Cif2Crystal.check_status()
  # Possible fix: execute checks if element is RunCrystal, and edits
  # job_record['dft']['restart_from'] accordingly.
  job_record['control']['force_retry'] = False
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
      #We could do checking to make sure the 
      #definition hasn't changed.
      f=open('record.json','r')
      record_read=json.load(f)
      record['control']=record_read['control']
      # Read DFT input file and check available keys.
      # TODO Read QMC input file and check available keys.
      f.close()
    print("#######################ID",record['control']['id'])

    for element in element_list:
      status=element.check_status(record)
      print(element._name_,"status",status)
      if status=='not_started':
        status=element.run(record)
        print(element._name_,"status",status)
      if status=='not_finished':
        if record['control']['force_retry']:
          status=element.retry(record)
        else:
          record['control']['incomplete'] = True
        print(element._name_,"status",status)
      else:
        record['control']['incomplete'] = False
      if status != 'ok':
        break
      record=element.output(record) 

    f=open(jsonfile,'w')
    json.dump(record,f)
    os.chdir(currwd)
    updated_job_list.append(record)
  return updated_job_list

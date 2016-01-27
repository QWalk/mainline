from __future__ import print_function
import os
import json
import shutil
def default_job_record(filename):
  job_record={}
  job_record['dft']={}
  job_record['qmc']={}
  job_record['control']={}

  #Set up Hamiltonian

  with open (filename, "r") as f:
    suffix=filename.split('.')[-1]
    if suffix=='cif':
      job_record['cif']=f.read()
    elif suffix=='xyz':
      job_record['xyz']=f.read()
    else:
      print("ERROR: didn't understand file suffix",suffix)
      quit()
  
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
  job_record['qmc']['dmc']['save_trace']=False


  job_record['qmc']['variance_optimize']={}
  job_record['qmc']['variance_optimize']['niterations']=10
  job_record['qmc']['variance_optimize']['nruns']=3


  job_record['qmc']['energy_optimize']={}
  job_record['qmc']['energy_optimize']['threshold']=0.001
  job_record['qmc']['energy_optimize']['vmc_nstep']=1000


  #Control options
  job_record['control']['id']=1
  job_record['control']['elements']=[]
  job_record['control']['pretty_formula']=''
  job_record['control']['queue_id']=[None,None]
  return job_record

def execute(record, element_list):
  """
  Run element_list tasks on this job record
  """
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
    if status != 'ok':
      break
    record=element.output(record)

  with open(jsonfile,'w') as f:
    json.dump(record,f)
  os.chdir(currwd)
  return record

def restart_job(jobname):
  do_it = raw_input("Restart %s? (y/n)"%jobname)
  if do_it == 'y':
    try:
      os.remove(jobname + "/autogen.d12")
    except OSError:
      pass
    try:
      os.remove(jobname + "/autogen.d12.o")
    except OSError:
      pass
  else:
    print("Didn't do it")

# Currently only defined for CRYSTAL runs.
def check_continue(jobname,qchecker,reasonable_lastSCF=50.0):
  jobrecord = json.load(open(jobname+"/record.json",'r'))
  qstatus = qchecker.status(jobrecord)
  if qstatus == 'running': 
    print("Shouldn't continue %s because still running"%jobname)
    return "running"
  try:
    outf = open(jobname+"/autogen.d12.o",'r')
  except IOError:
    print("Can't continue %s because no output"%jobname)
    return False
  outlines = outf.read().split('\n')
  reslines = [line for line in outlines if "ENDED" in line]
  if len(reslines) > 0:
    print("Shouldn't continue %s because finished."%jobname)
    return "finished"
  detots = [float(line.split()[5]) for line in outlines if "DETOT" in line]
  if len(detots) == 0:
    print("Shouldn't continue %s because no SCF last time."%jobname)
    return False
  detots_net = sum(detots[1:])
  if detots_net > reasonable_lastSCF:
    print("Shouldn't continue %s because unreasonable last try (%.2f>%.2f)."%\
        (jobname,detots_net,reasonable_lastSCF))
    return "unreasonable"
  etots = [float(line.split()[3]) for line in outlines if "DETOT" in line]
  if etots[-1] > 0:
    # This case probably won't happen if this works as expected.
    print("Shouldn't continue %s because divergence (%.2f)."%\
        (jobname,etots[-1]))
    return "divergence"
  print("Should continue %s."%jobname)
  return "continue"

# Currently only defined for CRYSTAL runs. Returns name of restart file.
def continue_job(jobname):
  jobrecord = json.load(open(jobname+"/record.json",'r'))
  trynum = 0
  while os.path.isfile(jobname+"/"+str(trynum)+".autogen.d12.o"):
    trynum += 1
  for filename in ["autogen.d12","autogen.d12.o","fort.79"]:
    shutil.move(jobname+"/"+filename,jobname+"/"+str(trynum)+"."+filename)
  return jobname+"/"+str(trynum)+".fort.79"

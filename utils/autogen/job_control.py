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
  job_record['dft']['smear']=None
  # Values:
  # 'stubborn'     : resume if job is killed or ran out of SCF steps.
  # 'optimistic'   : resume if job is killed.
  # 'conservative' : never resume job.
  job_record['dft']['resume_mode']='conservative'

  #QMC-specific options
  job_record['qmc']['dmc']={}
  job_record['qmc']['dmc']['timestep']=[0.02]
  job_record['qmc']['dmc']['jastrow']=['twobody'] #or 'threebody'
  job_record['qmc']['dmc']['nblock']=16
  job_record['qmc']['dmc']['optimizer']='variance' #or energy
  job_record['qmc']['dmc']['localization']=['tmoves']
  job_record['qmc']['dmc']['target_error']=0.01
  job_record['qmc']['dmc']['kpoints']='real' # or a list of numbers
  job_record['qmc']['dmc']['excitations']='no'#VBM-CBM or other..
  job_record['qmc']['dmc']['save_trace']=True


  job_record['qmc']['variance_optimize']={}
  job_record['qmc']['variance_optimize']['niterations']=10
  job_record['qmc']['variance_optimize']['nruns']=3
  job_record['qmc']['variance_optimize']['reltol']=0.1
  job_record['qmc']['variance_optimize']['abstol']=1e3 # TODO better default.
  job_record['qmc']['variance_optimize']['jastrow']=['twobody']


  job_record['qmc']['energy_optimize']={}
  job_record['qmc']['energy_optimize']['threshold']=0.001
  job_record['qmc']['energy_optimize']['vmc_nstep']=1000
  job_record['qmc']['energy_optimize']['jastrow']=['twobody']



  job_record['qmc']['maximize'] = {}
  job_record['qmc']['maximize']['nconfig'] = [100]
  job_record['qmc']['maximize']['jastrow']=['twobody']

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
      status=element.resume(record)
      print(element._name_,"status",status)
    if status != 'ok':
      break
    record=element.output(record)

  with open(jsonfile,'w') as f:
    json.dump(record,f)
  os.chdir(currwd)
  return record

# This will be moved to RunCrystal after our bigger merge.
def restart_job(jobname):
  """
  Restart a crystal job from scratch. This means deleting all progress. Use it
  for redoing a job that may have been corrupted or you'd like to change
  something important.
  """
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
# This will be moved to RunCrystal after our bigger merge.
def check_continue(jobname,qchecker,reasonable_lastSCF=50.0):
  """
  Look at CRYSTAL output, and report results.

  Current return values:
  no_record, running, no_output, success, too_many_cycles, finished (fall-back),
  scf_fail, not_enough_decrease, divergence, continue

  "continue" suggests the calculation should call continue_job(), and is only
  returned when no other condition is found.
  """
  try:
    jobrecord = json.load(open(jobname+"/record.json",'r'))
  except IOError:
    print("JOB CONTROL: Shouldn't continue %s has no record."%jobname)
    return "no_record"
  qstatus = qchecker.status(jobrecord)
  if qstatus == 'running': 
    print("JOB CONTROL: Shouldn't continue %s because still running"%jobname)
    return "running"
  try:
    outf = open(jobname+"/autogen.d12.o",'r')
  except IOError:
    print("JOB CONTROL: Can't continue %s because no output"%jobname)
    return "no_output"
  outlines = outf.read().split('\n')
  reslines = [line for line in outlines if "ENDED" in line]
  if len(reslines) > 0:
    if "CONVERGENCE" in reslines[0]:
      print("JOB CONTROL: Shouldn't continue %s because successful."%jobname)
      return "success"
    elif "TOO MANY CYCLES" in reslines[0]:
      print("JOB CONTROL: check_continue found %s has 'too many cycles'."%jobname)
      return "too_many_cycles"
    else: # What else can happen?
      return "finished"
  detots = [float(line.split()[5]) for line in outlines if "DETOT" in line]
  if len(detots) == 0:
    print("JOB CONTROL: Shouldn't continue %s because no SCF last time."%jobname)
    return "scf_fail"
  detots_net = sum(detots[1:])
  if detots_net > reasonable_lastSCF:
    print("JOB CONTROL: Shouldn't continue %s because not enough decrease (%.2f>%.2f)."%\
        (jobname,detots_net,reasonable_lastSCF))
    return "not_enough_decrease"
  etots = [float(line.split()[3]) for line in outlines if "DETOT" in line]
  if etots[-1] > 0:
    # This case probably won't happen if this works as expected.
    print("JOB CONTROL: Shouldn't continue %s because divergence (%.2f)."%\
        (jobname,etots[-1]))
    return "divergence"
  print("JOB CONTROL: Should continue %s."%jobname)
  return "continue"

# Currently only defined for CRYSTAL runs. Returns name of restart file.
# This will be moved to RunCrystal after our bigger merge.
# TODO: Add max_continues option.
# TODO: The stdout of this is off because the line between jobs is drawn between
#       this output and the execute() output.
def continue_job(jobname):
  """ Continue a job that ran out of time."""
  jobrecord = json.load(open(jobname+"/record.json",'r'))
  trynum = 0
  while os.path.isfile(jobname+"/"+str(trynum)+".autogen.d12.o"):
    trynum += 1
  for filename in ["autogen.d12","autogen.d12.o"]:
    shutil.move(jobname+"/"+filename,jobname+"/"+str(trynum)+"."+filename)
  for filename in ["fort.79"]:
    shutil.copy(jobname+"/"+filename,jobname+"/"+str(trynum)+"."+filename)
  return jobname+"/"+str(trynum)+".fort.79"

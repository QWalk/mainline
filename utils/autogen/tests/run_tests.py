import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control as jc
import os
import json

import taub

element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:05:00",queue="test")))
element_list.append(runcrystal.RunProperties(
  submitter=taub.LocalTaubPropertiesSubmitter(
    nn=1,np=1,time="0:05:00",queue="test")))
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(
  submitter=taub.LocalTaubQwalkSubmitter(
    nn=1,time="0:05:00",queue="test")))
#element_list.append(runqwalk.QWalkEnergyOptimize(
#  submitter=taub.LocalTaubQwalkSubmitter(
#    nn=1time="0:10:00",queue="test")))
element_list.append(runqwalk.QWalkRunDMC(
  submitter=taub.LocalTaubBundleQwalkSubmitter(
    nn=2,time="0:10:00",queue="secondary")))

default_job=jc.default_job_record("si.cif")
default_job['dft']['kmesh'] = [2,2,2]
default_job['dft']['functional']['hybrid'] = 0
default_job['dft']['tolinteg'] = [6,6,6,6,12]
default_job['dft']['basis']=[0.2,2,2]
default_job['dft']['maxcycle'] = 100
default_job['dft']['fmixing'] = 80
default_job['dft']['edifftol'] = 6
default_job['dft']['broyden'] = [0.1,60,20]
default_job['qmc']['variance_optimize']['reltol']=0.1
default_job['qmc']['variance_optimize']['abstol']=10
default_job['qmc']['dmc']['save_trace'] = False
default_job['qmc']['dmc']['nblock']=5
default_job['qmc']['dmc']['target_error']=0.1
default_job['total_spin'] = 0
idbase = "test_si_"

count=1

# Simple run.
name = idbase+"simple"
job_record = copy.deepcopy(default_job)
job_record['control']['id']=name
jc.execute(job_record,element_list)
count+=1

# Restart DFT and change something. Also check redefinition of element_list
element_list[1] = runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:20:00",queue="secondary"
  ))
if os.path.getsize("%s/fort.9"%(idbase+"simple")) > 0:
  name = idbase+"edit"
  job_record = copy.deepcopy(default_job)
  job_record['dft']['functional']['hybrid'] = 25
  job_record['dft']['restart_from'] = "../%s/fort.9"%(idbase+"simple")
  job_record['control']['id']=name
  jc.execute(job_record,element_list)
  count+=1
element_list[1] = runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:05:00",queue="test"))

# Too-many cycles case.
name = idbase+"toomany"
job_record = copy.deepcopy(default_job)
job_record['control']['id']   = name
job_record['dft']['maxcycle'] = 10
job_record['dft']['edifftol'] = 10
jc.execute(job_record,element_list)
count+=1

name = idbase+"toomany_resumed"
job_record = copy.deepcopy(default_job)
job_record['control']['id']   = name
job_record['dft']['maxcycle'] = 10
job_record['dft']['edifftol'] = 10
job_record['dft']['resume_mode'] = 'stubborn'
jc.execute(job_record,element_list)
count+=1

# Reduce allowed run time.
element_list[1] = runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:00:30",queue="test"
  ))

# Demonstrate dft killed-job error correction.
name = idbase+"killed"
job_record = copy.deepcopy(default_job)
job_record['control']['id']   = name
job_record['dft']['edifftol'] = 14
jc.execute(job_record,element_list)
count+=1

name = idbase+"killed_revived"
job_record = copy.deepcopy(default_job)
job_record['control']['id']     = name
job_record['dft']['resume_mode'] = 'optimistic'
job_record['dft']['edifftol'] = 14
jc.execute(job_record,element_list)
count+=1

# Reset allowed run time.
element_list[1] = runcrystal.RunCrystal(
  submitter=taub.LocalTaubCrystalSubmitter(
    nn=1,time="0:05:00",queue="test"
  ))

# Variance optimize takes many tries.
name = idbase+"retryvar"
job_record = copy.deepcopy(default_job)
job_record['control']['id']=name
job_record['qmc']['variance_optimize']['reltol']=0.001
jc.execute(job_record,element_list)
count+=1

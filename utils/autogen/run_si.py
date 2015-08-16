from __future__ import print_function
import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control
#import data_processing as dp
import os
import json
import pxssh
import pexpect
import submission_tools

#hawk.py contains submitters for hawk.physics.illinois.edu
#from hawk import * 

#from veritas import * 

job_record=job_control.default_job_record("si.cif")
job_record['dft']['kmesh'] = [2,2,2]
job_record['dft']['functional']['hybrid'] = 0
job_record['dft']['tolinteg'] = [6,6,6,6,12]
job_record['dft']['maxcycle'] = 100
job_record['dft']['fmixing'] = 90
job_record['dft']['broyden'] = [0.1,60,20]
#job_record['control']['force_retry'] = True
job_record['total_spin'] = 0
idbase = "si_ag_"

#An example of varying basis parameters
job_list=[]
count=1
for alpha in [0.1]:
  job_record['dft']['basis']=[alpha,3,3]
  job_record['control']['id']=idbase+str(count)
  job_list.append(copy.deepcopy(job_record))
  count+=1

# for alpha in [0.1,0.2,0.3]:
#   job_record['dft']['basis']=[alpha,3,3]
#   job_record['control']['id']=idbase+str(count)
#   job_record['dft']['maxcycle'] = 10
#   #job_record['dft']['restart_from'] = "../si_ag_3/fort.9"
#   job_list.append(copy.deepcopy(job_record))
#   count+=1

#for job in job_list:
#  if job['control']['id'] in ['si_4','si_5','si_6']: 
  
#for alpha in [0.1,0.2,0.3]:
#  job_record['dft']['basis']=[alpha,2,2]
#  job_record['control']['id']=count
#  job_list.append(copy.deepcopy(job_record))
#  count+=1
#  
#for alpha in [0.1,0.2,0.3]:
#  job_record['dft']['basis']=[alpha,1,2]
#  job_record['control']['id']=count
#  job_list.append(copy.deepcopy(job_record))
#  count+=1
#
##after running the above basis convergence, we select the minimal basis that gives a low energy
##and calculate at a larger supercell
#job_record['dft']['basis']=[0.1,1,2]
#job_record['control']['id']=count
#job_record['supercell']=[[2,0,0],[0,2,0],[0,0,2]]
#job_list.append(copy.deepcopy(job_record))
#count+=1
#
##We decided to go back and check the mixing percentage for the primitive cell
#job_record['supercell']=[[1,0,0],[0,1,0],[0,0,1]]
#
#for mixing in [0,10,20,30,40]: 
#  job_record['control']['id']=count
#  job_record['dft']['functional']['hybrid']=mixing
#  job_list.append(copy.deepcopy(job_record))
#  count+=1


# PX Tools #
ssh_taub = pxssh.pxssh()
ssh_taub.login('taub.campuscluster.illinois.edu', 'jaschil2')
ftp_taub = pexpect.spawn('sftp jaschil2@taub.campuscluster.illinois.edu')

# Submitters #
sub_crystal = submission_tools.RemoteSubmitter(
  ssh_taub, 
  ftp_taub, 
  'runcrystal_taub', 
  '/projects/erg/jaschil2/QMCDB/')

sub_qmc_var = submission_tools.RemoteSubmitter(
  ssh_taub,
  ftp_taub,
  'run_v_opt_taub',
  '/projects/erg/jaschil2/QMCDB/')

sub_qmc_en = submission_tools.RemoteSubmitter(
  ssh_taub,
  ftp_taub,
  'run_en_opt_taub',
  '/projects/erg/jaschil2/QMCDB/')

sub_dmc = submission_tools.RemoteSubmitter(
  ssh_taub,
  ftp_taub,
  'run_dmc_taub',
  '/projects/erg/jaschil2/QMCDB/')

#now we define the job sequence
element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(sub_crystal))
element_list.append(runcrystal.RunProperties())
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(sub_qmc_var))
element_list.append(runqwalk.QWalkEnergyOptimize(sub_qmc_en))
element_list.append(runqwalk.QWalkRunDMC(sub_dmc))

job_list = job_control.execute(job_list,element_list)

#json.dump(dp.trace_analysis(
#  [job['control']['id']+'/autogen.d12.o' for job in job_list],
#  ids = [job['control']['id'] for job in job_list]),
#  open('broyden.json','w'))

#Print out a summary of the results
print("DFT energy")
for result in job_list:
  if 'total_energy' in result['dft'].keys():
    print(result['dft']['basis'],result['dft']['total_energy'])

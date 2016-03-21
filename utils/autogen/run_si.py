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

#An example to submit to the Taub campus cluster.
#Change these 
directory="this/is/the/remote/directory/" #the trailing slash is important!
username="johndoe3"
host = 'taub.campuscluster.illinois.edu'


job_record=job_control.default_job_record("si.cif")
job_record['dft']['kmesh'] = [4,4,4]
job_record['dft']['functional']['hybrid'] = 0
job_record['dft']['tolinteg'] = [6,6,6,6,12]
job_record['dft']['maxcycle'] = 100
job_record['dft']['fmixing'] = 90
job_record['dft']['broyden'] = [0.1,60,10]
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


# PX Tools #
#ssh_taub = pxssh.pxssh()
#ssh_taub.login(host, username)
#ftp_taub = pexpect.spawn('sftp '+ username + '@' + host)
#ftp_taub.delaybeforesend = 0.5

#ssh_taub.setecho(True)
#ftp_taub.setecho(True)

# Submitters #
#sub_crystal = submission_tools.RemoteSubmitter(
#  ssh_taub, 
#  ftp_taub, 
#  'subcrystal_taub', 
#  directory)

#sub_qmc_var = submission_tools.RemoteSubmitter(
#  ssh_taub,
#  ftp_taub,
#  'subqwalk_taub',
#  directory)

#sub_qmc_en = submission_tools.RemoteSubmitter(
#  ssh_taub,
#  ftp_taub,
#  'subqwalk_taub',
#  directory)

#sub_dmc = submission_tools.RemoteSubmitter(
#  ssh_taub,
#  ftp_taub,
#  'subqwalk_taub',
#  directory)

#now we define the job sequence
element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
#element_list.append(runcrystal.RunCrystal(sub_crystal))
#element_list.append(runcrystal.RunProperties())
#element_list.append(runqwalk.Crystal2QWalk())
#element_list.append(runqwalk.QWalkVarianceOptimize(sub_qmc_var))
#element_list.append(runqwalk.QWalkEnergyOptimize(sub_qmc_en))
#element_list.append(runqwalk.QWalkRunDMC(sub_dmc))

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

from __future__ import print_function
import sys
sys.path.append("../")
import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control
import os
import json
import yaml
import pexpect
from pexpect import pxssh
import submission_tools

#An example to submit to the Taub campus cluster.
#Change these 
directory="/this/is/the/remote/directory/" #the trailing slash is important!
username="johndoe3"
#directory="/home/lkwagner/autogen_tmp"
#username="lkwagner"
host = 'hawk.physics.illinois.edu'


# Here we set up our jobs. 
job_list=[]

#This is a 2-atom silicon run.
job_record=job_control.default_job_record("si.cif")
job_record['dft']['kmesh'] = [4,4,4]
job_record['dft']['functional']['hybrid'] = 0
job_record['dft']['tolinteg'] = [6,6,6,6,12]
job_record['dft']['maxcycle'] = 100
job_record['dft']['fmixing'] = 90
job_record['dft']['broyden'] = [0.1,60,10]
#job_record['control']['force_retry'] = True
job_record['total_spin'] = 0
job_record['dft']['basis']=[0.1,3,3]
job_record['control']['id']="silicon"
job_list.append(copy.deepcopy(job_record))


# PX Tools #
ssh_hawk = pxssh.pxssh()
ssh_hawk.login(host, username)
ftp_hawk = pexpect.spawn('sftp '+ username + '@' + host)
ftp_hawk.delaybeforesend = 0.5
ssh_hawk.setecho(True)
ftp_hawk.setecho(True)

# Submitters #
sub_crystal = submission_tools.RemoteSubmitter(
  ssh_hawk, ftp_hawk, 
  'subcrystal_hawk', 
  directory)

sub_qwalk = submission_tools.RemoteSubmitter(
  ssh_hawk,ftp_hawk,
  'subqwalk_hawk',
  directory)

#now we define the job sequence
element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(sub_crystal))
element_list.append(runcrystal.RunProperties())
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(sub_qwalk))
element_list.append(runqwalk.QWalkEnergyOptimize(sub_qwalk))
element_list.append(runqwalk.QWalkRunDMC(sub_qwalk))


#This executes the job
job_list = job_control.execute(job_list,element_list)


#Save the data, either in JSON:
json.dump(job_list,open("data.json",'w'))
#or in YAML (which is easier to read for some)
yaml.dump(job_list,open("data.yaml",'w'),default_flow_style=False)


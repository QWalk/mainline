from __future__ import print_function
import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control

#hawk.py contains submitters for hawk.physics.illinois.edu
#from hawk import * 

from remotetaub import * 

job_record=job_control.default_job_record("si.cif")

#An example of varying basis parameters
job_list=[]
count=1
for alpha in [0.1,0.2,0.3]:
  job_record['dft']['basis']=[alpha,3,3]
  job_record['control']['id']=count
  job_list.append(copy.deepcopy(job_record))
  count+=1
for alpha in [0.1,0.2,0.3]:
  job_record['dft']['basis']=[alpha,2,3]
  job_record['control']['id']=count
  job_list.append(copy.deepcopy(job_record))
  count+=1
  
for alpha in [0.1,0.2,0.3]:
  job_record['dft']['basis']=[alpha,2,2]
  job_record['control']['id']=count
  job_list.append(copy.deepcopy(job_record))
  count+=1
  
for alpha in [0.1,0.2,0.3]:
  job_record['dft']['basis']=[alpha,1,2]
  job_record['control']['id']=count
  job_list.append(copy.deepcopy(job_record))
  count+=1

#after running the above basis convergence, we select the minimal basis that gives a low energy
#and calculate at a larger supercell
job_record['dft']['basis']=[0.1,1,2]
job_record['control']['id']=count
job_record['supercell']=[[2,0,0],[0,2,0],[0,0,2]]
job_list.append(copy.deepcopy(job_record))
count+=1

#We decided to go back and check the mixing percentage for the primitive cell
job_record['supercell']=[[1,0,0],[0,1,0],[0,0,1]]

for mixing in [0,10,20,30,40]: 
  job_record['control']['id']=count
  job_record['dft']['functional']['hybrid']=mixing
  job_list.append(copy.deepcopy(job_record))
  count+=1


#now we define the job sequence
element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(submitter=MyTorqueCrystalSubmitter()))
element_list.append(runcrystal.RunProperties())
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(submitter=MyTorqueQWalkSubmitter()))
element_list.append(runqwalk.QWalkRunDMC(submitter=MyTorqueQWalkSubmitter()))

job_list=job_control.execute(job_list,element_list)


#Print out a summary of the results
print("DFT energy")
for result in job_list:
  if 'total_energy' in result['dft'].keys():
    print(result['dft']['basis'],result['dft']['total_energy'])
print("QMC energy")
for result in job_list:
  if 'total_energy' in result['qmc'].keys():
    print(result['dft']['functional']['hybrid'],result['dft']['basis'],
            result['qmc']['total_energy'],result['qmc']['total_energy_err'])


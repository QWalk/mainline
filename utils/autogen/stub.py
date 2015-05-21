from __future__ import print_function
import cif2crystal
import runcrystal
import runqwalk

job_record={}
job_record['dft']={}
job_record['qmc']={}
job_record['control']={}
##This is a testing stub that will set up a single calculation
with open ("si.cif", "r") as f:
    job_record['cif']=f.read()
job_record['supercell']=[[1,0,0],[0,1,0],[0,0,1]]
job_record['pseudopotential']='BFD'
job_record['charge']=0
job_record['total_spin']=0

#DFT-specific options
job_record['dft']['functional']={'exchange':'PBE','correlation':'PBE','hybrid':25}
job_record['dft']['basis']=[0.2,3,3]
job_record['dft']['kmesh']=[8,8,8]
job_record['dft']['spin_polarized']=False
job_record['dft']['initial_spin']=[]

#QMC-specific options
job_record['qmc']['timestep']=0.02
job_record['qmc']['jastrow']='twobody'
job_record['qmc']['optimize']='variance'
job_record['qmc']['localization']='tmoves'


#Control options
job_record['control']['id']=1
job_record['control']['elements']=['Si']
job_record['control']['pretty_formula']='Si'


#now we define the job sequence
element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal())
element_list.append(runcrystal.RunProperties())
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize())
element_list.append(runqwalk.QWalkRunDMC())


for element in element_list:
  status=element.check_status(job_record)
  print("status",status)
  if status=='not_started':
    element.run(job_record)
  elif status=='running':
    break
  elif status=='ok':
    job_record=element.output(job_record)
  else:
    print("Got unknown status:",status)
    quit()

print(job_record)


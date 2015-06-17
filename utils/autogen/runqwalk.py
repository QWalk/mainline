import os
import job_submission
####################################################

def crystal_patch_output(propname,outname,patchname):
  prop=open(propname,'r')
  shrink=[1,1,1]
  for line in prop:
    if "SHRINK FACTORS(MONK.)" in line:
      spl=line.split()
      shrink[0]=int(spl[2])
      shrink[1]=int(spl[3])
      shrink[2]=int(spl[4])


  patch=open(patchname,'w')

  out=open(outname,'r')
  for line in out:
    if "SHRINK. FACT.(MONKH.)" in line:
      patch.write("SHRINK. FACT.(MONKH.)    %i  %i  %i  NUMBER OF K POINTS IN THE IBZ    XXX\n"%(shrink[0],shrink[1],shrink[2]))
    else:
      patch.write(line)
    
    if "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT EDFT" in line:
      break
  out.close()

  prop=open(propname,'r')
  patch.write("NEWK EIGENVECTORS\n \n")
  found_hamil=False
  for line in prop:
    if "HAMILTONIAN EIGENVECTORS" in line:
      found_hamil=True
    if found_hamil:
      patch.write(line)

####################################################



class Crystal2QWalk:
  _name_="Crystal2QWalk"
  def run(self,job_record):
    crystal_patch_output("prop.in.o","autogen.d12.o","patched.o")
    os.system("crystal2qmc -o qw patched.o > crystal2qmc.stdout")
    return 'ok'
  def check_status(self,job_record):
    outfilename="qw_0.sys"
    if os.path.exists(outfilename):
      return 'ok'
    return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    return job_record

####################################################

class QWalkVarianceOptimize:
  _name_="QWalkVarianceOptimize"
  _submitter=job_submission.TorqueQWalkSubmitter()
  
  def __init__(self,submitter=job_submission.TorqueQWalkSubmitter()):
    self._submitter=submitter
  
  def run(self,job_record):
    f=open("qw_0.opt",'w')
    f.write("""method { optimize } 
include qw_0.sys
trialfunc { slater-jastrow
wf1 { include qw_0.slater } 
wf2 { include qw.jast2 } 
}
""")
    f.close()
    self._submitter.execute(job_record,
            ['qw_0.opt','qw_0.sys','qw_0.slater','qw_0.orb','qw.basis','qw.jast2'])
    
    return 'running'


  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if 'Wall' in line:
          return 'ok'
      return 'running'


  def check_status(self,job_record):
    outfilename="qw_0.opt.o"
      
    if self.check_outputfile(outfilename)=='ok':
      return 'ok'
    status=self._submitter.status(job_record,[outfilename,'qw_0.opt.wfout'])
    if status=='running':
      return status
    if self.check_outputfile(outfilename)=='ok':
      return 'ok'
      
  
    return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    outfilename="qw_0.opt.o"
    f=open(outfilename,'r')
    disp=0.00
    for line in f:
      if 'iteration' in line and 'dispersion' in line:
        spl=line.split()
        disp=float(spl[4])

    job_record['qmc']['optimized_variance']=disp
    return job_record

####################################################
class QWalkEnergyOptimize:
  _name_="QWalkEnergyOptimize"
  _submitter=job_submission.TorqueQWalkSubmitter()
  
  def __init__(self,submitter=job_submission.TorqueQWalkSubmitter()):
    self._submitter=submitter
  
  def run(self,job_record):
    if not os.path.isfile("qw_0.opt.wfout"):
      print("Could not find qw_0.opt.wfout")
      return "failed"
    os.system("sed s/OPTIMIZEBASIS//g qw_0.opt.wfout > qw_0.enopt.wfin")
    f=open("qw_0.enopt",'w')
    f.write("""method { LINEAR VMC_NSTEP 5000 } 
include qw_0.sys
trialfunc { include qw_0.enopt.wfin }
""")
    f.close()
    self._submitter.execute(job_record,
            ['qw_0.enopt','qw_0.enopt.wfin','qw_0.sys','qw_0.slater','qw_0.orb','qw.basis'])
    
    return 'running'


  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if 'Wall' in line:
          return 'ok'
      return 'running'


  def check_status(self,job_record):
    outfilename="qw_0.enopt.o"
      
    if self.check_outputfile(outfilename)=='ok':
      return 'ok'
    status=self._submitter.status(job_record,[outfilename,'qw_0.enopt.wfout'])
    if status=='running':
      return status
    if self.check_outputfile(outfilename)=='ok':
      return 'ok'
      
  
    return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    outfilename="qw_0.opt.o"
    f=open(outfilename,'r')
    disp=0.00
    for line in f:
      if 'iteration' in line and 'dispersion' in line:
        spl=line.split()
        disp=float(spl[4])

    job_record['qmc']['vmc_energy']=disp
    return job_record



####################################################

class QWalkRunDMC:
  _name_="QwalkRunDMC"
  _submitter=job_submission.TorqueQWalkSubmitter()
  def __init__(self,submitter=job_submission.TorqueQWalkSubmitter()):
    self._submitter=submitter
#-----------------------------------------------
  def run(self,job_record,restart=False):
    qmc_options=job_record['qmc']
    f=open("qw_0.dmc",'w')
    f.write("""method { DMC timestep %g nblock 16 %s  } 
include qw_0.sys
trialfunc { slater-jastrow
wf1 { include qw_0.slater } 
wf2 { include opt.jast } 
}
"""%(qmc_options['timestep'],qmc_options['localization']))
    f.close()
    os.system("separate_jastrow qw_0.opt.wfout > opt.jast")
    infiles=["qw_0.dmc","opt.jast",'qw_0.sys','qw_0.slater','qw_0.orb','qw.basis']
    if restart:
      infiles.extend(['qw_0.dmc.config','qw_0.dmc.log'])
    self._submitter.execute(job_record,infiles)
    return 'running'
#-----------------------------------------------
  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if 'Wall' in line:
          return 'ok'
      return 'running'
#-----------------------------------------------
  def energy(self,job_record):
    os.system("gosling qw_0.dmc.log > qw_0.dmc.log.stdout")
    f=open("qw_0.dmc.log.stdout")
    energy=0
    err=1e8
    for line in f:
      if "total_energy0" in line:
        spl=line.split()
        energy=float(spl[1])
        err=float(spl[3])
        return (energy,err)
    return (energy,err)

#-----------------------------------------------

  def check_status(self,job_record):
    outfilename="qw_0.dmc.o"
    #if self.check_outputfile(outfilename)=='ok':
    #  return 'ok'
    energy,err=self.energy(job_record)
    print("check",energy,err)
    
    if err< job_record['qmc']['target_error']:
      return 'ok'
    status=self._submitter.status(job_record,[outfilename,'qw_0.dmc.log','qw_0.dmc.config'])
    if status=='running':
      return status
    if not os.path.isfile("qw_0.dmc.log"):
      return 'not_started'

    energy,err=self.energy(job_record)
    if err < job_record['qmc']['target_error']:
      return 'ok'
    else:
      return 'not_finished'
#-----------------------------------------------
    
      
  def retry(self,job_record):
    return self.run(job_record,restart=True)
#-----------------------------------------------

  def output(self,job_record):
    energy,err=self.energy(job_record)
    job_record['qmc']['total_energy']=energy
    job_record['qmc']['total_energy_err']=err
    return job_record
####################################################


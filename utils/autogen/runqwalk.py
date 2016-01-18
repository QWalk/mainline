from __future__ import print_function
import os
import glob
import re
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
  
  def __init__(self,submitter):
    self._submitter=submitter
  
  def run(self,job_record):
    f=open("qw_0.opt",'w')
    nit=job_record['qmc']['variance_optimize']['niterations']
    nruns=job_record['qmc']['variance_optimize']['nruns']
    for i in range(0,nruns):
        f.write("method { optimize iterations %i } "%nit)
    f.write("""
include qw_0.sys
trialfunc { slater-jastrow
wf1 { include qw_0.slater } 
wf2 { include qw.jast2 } 
}
""")
    f.close()
    job_record['control'][self._name_+'_jobid'] = [self._submitter.execute(
      job_record, 
      ['qw_0.opt','qw_0.sys','qw_0.slater','qw_0.orb','qw.basis','qw.jast2'], 
      'qw_0.opt',
      'qw_0.opt.stdout')]
    
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
    status=self._submitter.status(job_record)
    self._submitter.output(job_record, ['qw_0.opt.o', 'qw_0.opt.wfout'])  
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
  
  def __init__(self,submitter):
    self._submitter=submitter
  
  def run(self,job_record,restart=False):
    if not os.path.isfile("qw_0.opt.wfout"):
      print("Could not find qw_0.opt.wfout")
      return "failed"
    if restart:
      os.system("cp qw_0.enopt.wfout qw_0.enopt.wfin")
    else:
      os.system("sed s/OPTIMIZEBASIS//g qw_0.opt.wfout > qw_0.enopt.wfin")

    enopt_options=job_record['qmc']['energy_optimize']

    f=open("qw_0.enopt",'w')
    f.write("""method { LINEAR VMC_NSTEP %i } 
include qw_0.sys
trialfunc { include qw_0.enopt.wfin }
"""%enopt_options['vmc_nstep'])
    f.close()
    job_record['control'][self._name_+'_jobid'] = [self._submitter.execute(
      job_record, 
      ['qw_0.enopt','qw_0.enopt.wfin','qw_0.sys','qw_0.slater','qw_0.orb','qw.basis'],
      'qw_0.enopt',
      'qw_0.enopt.stdout')]
    
    return 'running'


  def check_outputfile(self,outfilename, threshold=0.001):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      last_change=1e8
      for line in f:
        if 'Wall' in line:
          return 'ok'
        if 'step' in line and 'current energy' in line:
          spl=line.split()
          if len(spl) > 9:
            last_change=float(spl[9])
      print("energy optimize: last change",last_change, 'threshold',threshold)
      if abs(last_change) > threshold:
        return 'not_finished'
      return 'ok'
    return 'not_started'


  def check_status(self,job_record):
    outfilename="qw_0.enopt.o"
    self._submitter.output(job_record, [outfilename, 'qw_0.enopt.wfout'])
      
    status=self._submitter.status(job_record)
    if status=='running':
      return status

    if self.check_outputfile(outfilename)=='ok':
      return 'ok'
    
    status=self._submitter.status(job_record)
    if status=='running':
      return status
    return self.check_outputfile(outfilename,
            threshold=job_record['qmc']['energy_optimize']['threshold'])
      
  def retry(self,job_record):
    return self.run(job_record,restart=True)

  def output(self,job_record):
    outfilename="qw_0.enopt.o"
    f=open(outfilename,'r')
    last_energy=1e8
    last_energy_err=1e8
    for line in f:
      if 'step' in line:
        spl=line.split()
        if len(spl) > 9:
          last_energy=spl[4]
          last_energy_err=spl[6]

    job_record['qmc']['energy_optimize']['vmc_energy']=last_energy
    job_record['qmc']['energy_optimize']['vmc_energy_err']=last_energy_err
    return job_record



####################################################

class QWalkRunDMC:
  _name_="QwalkRunDMC"
  
  def __init__(self,submitter):
    self._submitter=submitter
#-----------------------------------------------
  def run(self,job_record,restart=False):
    qmc_options=job_record['qmc']
    kpts=self.get_kpts(job_record)
    if self._name_+'_jobid' not in job_record['control'].keys():
      job_record['control'][self._name_+'_jobid'] = []
      job_record['control']['queue_id'] = []
    
    #choose which wave function to use
    if not restart:
      if qmc_options['dmc']['optimizer']=='variance':
        os.system("separate_jastrow qw_0.opt.wfout > opt.jast")
      elif qmc_options['dmc']['optimizer']=='energy':
        os.system("separate_jastrow qw_0.enopt.wfout > opt.jast")

    #make and submit the runs.
    #this could be extended to bundle jobs
    for k in kpts:
      for t in qmc_options['dmc']['timestep']:
        for loc in qmc_options['dmc']['localization']:
          kname="qw_%i"%k
          basename="qwt%g%s_%i"%(t,loc,k)
          f=open(basename+".dmc",'w')
          f.write(self.dmcinput(t,loc,k,qmc_options['dmc']['nblock']))
          f.close()
          infiles=[basename+".dmc","opt.jast",kname+'.sys',kname+'.slater',
              kname+'.orb','qw.basis']
          if restart:
            infiles.extend([basename+'.dmc.config',basename+'.dmc.log'])
          job_id = self._submitter.execute(
            job_record,
            infiles,
            basename+".dmc",
            basename+".dmc.stdout")
          job_record['control'][self._name_+'_jobid'].append(job_id)

    job_record['control']['queue_id'] = job_record['control'][self._name_+'_jobid']
    return 'running'


#-----------------------------------------------
  def get_kpts(self, job_record):
    kpts=glob.glob("qw*.sys")
    kpt_num=[]
    for kp in kpts:
      kpt_num.append(int(re.findall(r'\d+',kp)[0]))
    return kpt_num


#-----------------------------------------------
  def dmcinput(self,timestep,localization,kpt_num=0,nblock=16):
    return """method { DMC timestep %g nblock %i %s 
    average { SK } 
} 
include qw_%i.sys
trialfunc { slater-jastrow
wf1 { include qw_%i.slater } 
wf2 { include opt.jast } 
}
"""%(timestep,nblock,localization,kpt_num,kpt_num)

#-----------------------------------------------
  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if 'Wall' in line:
          return 'ok'
      return 'running'
#-----------------------------------------------
  def energy(self,logfilename):
    os.system("gosling %s > %s.stdout"%(logfilename,logfilename))
    f=open(logfilename+".stdout",'r')
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
  def sk(self,logfilename):
    os.system("gosling %s | grep 'sk_out' | sed 's/sk_out//g' > %s_sk.out"%(logfilename,logfilename))
    SK=[]
    f=open("%s_sk.out" %logfilename)
    line=f.readline()
    while len(line)>0:
      SK.append([float(d) for d in line.split()[1:]])
      line = f.readline()
    f.close()
    return SK
#  def finitesize_correction(SK):
    
#-----------------------------------------------

  def collect_runs(self,job_record):
    ret=[]

    kpts=self.get_kpts(job_record)
    for k in kpts:
      for t in job_record['qmc']['dmc']['timestep']:
        for loc in job_record['qmc']['dmc']['localization']:
          kname="qw_%i"%k
          basename="qwt%g%s_%i"%(t,loc,k)
          entry={}
          entry['knum']=k
          entry['energy']=self.energy(basename+".dmc.log")
          entry['timestep']=t
          entry['localization']=loc
          entry['sk']=self.sk(basename+".dmc.log")
          ret.append(entry)
    return ret
          

    

#-----------------------------------------------

  def check_status(self,job_record):
    
    results=self.collect_runs(job_record)
    
    status='ok'
    for e in results:
      print("energy check",e['knum'],e['energy'])
        
      if e['energy'][1] >  job_record['qmc']['dmc']['target_error']:
        status='not_finished'
    print("initial status",status)
    if status=='ok':
      return status


    outfiles=[]
    for k in self.get_kpts(job_record):
      for t in job_record['qmc']['dmc']['timestep']:
        for local in job_record['qmc']['dmc']['localization']:
          basename="qwt%g%s_%i"%(t,local,k)
          outfiles.extend([basename+'.dmc.log',
                         basename+'.dmc.config',
                         basename+'.dmc.o'])
    print(outfiles)
    self._submitter.output(job_record, outfiles)
    status=self._submitter.status(job_record)
    print("status",status)
    if status=='running':
      return status


    if not os.path.isfile(outfiles[0]):
      return 'not_started'

    results=self.collect_runs(job_record)
    status='ok'
    for e in results:
      if e['energy'][1] >  job_record['qmc']['dmc']['target_error']:
        status='not_finished'
    return status
  
#-----------------------------------------------
    
      
  def retry(self,job_record):
    return self.run(job_record,restart=True)
#-----------------------------------------------

  def output(self,job_record):
    job_record['qmc']['dmc']['results']=self.collect_runs(job_record)
    return job_record
####################################################


import os
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
  def run(self,job_record):
    d=str(job_record['control']['id'])+"/"
    crystal_patch_output(d+"prop.in.o",d+"autogen.d12.o",d+"patched.o")
    os.system("cd %s; crystal2qmc -o qw patched.o > crystal2qmc.stdout"%d)
    return 'ok'
  def check_status(self,job_record):
    outfilename=str(job_record['control']['id'])+"/qw_0.sys"
    if os.path.exists(outfilename):
      return 'ok'
    #this will just always run. Ideally should check to see if the files are there
    return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    return job_record

####################################################

class QWalkVarianceOptimize:
  def run(self,job_record):
    d=str(job_record['control']['id'])+"/"
    f=open(d+"qw_0.opt",'w')
    f.write("""method { optimize } 
include qw_0.sys
trialfunc { slater-jastrow
wf1 { include qw_0.slater } 
wf2 { include qw.jast2 } 
}
""")
    f.close()
    os.system("cd %s; qwalk qw_0.opt >& qw_0.opt.stdout "%d)
    os.system("cd %s; separate_jastrow qw_0.opt.wfout > opt.jast"%d)
    return 'ok'


  def check_status(self,job_record):
    outfilename=str(job_record['control']['id'])+"/qw_0.opt.o"
    if os.path.isfile(outfilename):
      if os.path.isfile(str(job_record['control']['id'])+"/opt.jast"):
        return 'ok'
      else:
        return 'running'
    return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    outfilename=str(job_record['control']['id'])+"/qw_0.opt.o"
    f=open(outfilename,'r')
    disp=0.00
    for line in f:
      if 'iteration' in line and 'dispersion' in line:
        spl=line.split()
        disp=float(spl[4])

    job_record['qmc']['optimized_variance']=disp
    return job_record


####################################################

class QWalkRunDMC:
  def run(self,job_record):
    d=str(job_record['control']['id'])+"/"
    f=open(d+"qw_0.dmc",'w')
    f.write("""method { DMC timestep .02 nblock 16 } 
include qw_0.sys
trialfunc { slater-jastrow
wf1 { include qw_0.slater } 
wf2 { include opt.jast } 
}
""")
    f.close()
    os.system("cd %s; qwalk qw_0.dmc >& qw_0.dmc.stdout "%d)
    return 'ok'

  def check_status(self,job_record):
    outfilename=str(job_record['control']['id'])+"/qw_0.dmc.o"
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if "Total wall" in line:
          return 'ok'
      return 'running'
    return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)

  def output(self,job_record):
    d=str(job_record['control']['id'])+"/"      
    os.system("cd %s; gosling qw_0.dmc.log > qw_0.dmc.log.stdout"%d)
    f=open(d+"qw_0.dmc.log.stdout")
    for line in f:
      if "total_energy0" in line:
        spl=line.split()
        job_record['qmc']['total_energy']=float(spl[1])
        job_record['qmc']['total_energy_err']=float(spl[3])
    return job_record
####################################################


from __future__ import print_function
import sys
sys.path.append('..') # Other QWalk utils.
import read_numberfluct as nf
import dm_tools as dm
import os
import glob
import re
import shutil
from crystal2qmc import convert_crystal
import numpy as np
import json
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

class NewCrystal2QWalk:
  _name_="NewCrystal2QWalk"
  def run(self,job_record):
    job_record['qmc']['kpoint_weights'] = \
        convert_crystal(
            base="qw",
            kfmt='int',
            kset=job_record['qmc']['kpoints']).tolist()
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

class Crystal2QWalk:
  _name_="Crystal2QWalk"
  def run(self,job_record):
    crystal_patch_output("prop.in.o","autogen.d12.o","patched.o")
    if job_record['qmc']['kpoints'] == "complex":
      os.system("crystal2qmc -c -o qw patched.o > crystal2qmc.stdout")
    else:
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
  _name_ = "QWalkVarianceOptimize"
  
  def __init__(self,submitter):
    self._submitter = submitter
#------------------------------------------
  def run(self,job_record):
    infiles=[]
    jastfiles=[]
    for jast in job_record['qmc']['variance_optimize']['jastrow']:
      jast_suf=""
      if jast=='twobody':
        jast_suf = 'jast2'
      elif jast=='threebody':
        jast_suf = 'jast3'
      else:
        print("Didn't understand Jastrow",jast)
        quit()
      fname="qw_0.%s.opt"%jast
      f = open(fname,'w')
      nit = job_record['qmc']['variance_optimize']['niterations']
      nruns = job_record['qmc']['variance_optimize']['nruns']
      for i in range(0,nruns):
        f.write("method { optimize iterations %i } "%nit)
      f.write("\n".join([
          "include qw_0.sys",
          "trialfunc { slater-jastrow",
          "  wf1 { include qw_0.slater } ",
          "  wf2 { include qw.%s } "%jast_suf,
          "}"
        ]))
      infiles.append(fname)
      jastfiles.append("qw.%s"%jast_suf)
      f.close()
    outfiles=[]
    for fname in infiles:
      outfiles.append(fname+".stdout")

    self._submitter.execute(
      job_record, 
      infiles+['qw_0.sys','qw_0.slater','qw_0.orb','qw.basis']+jastfiles, 
      infiles,
      outfiles[0],
      self._name_)
    
    return 'running'
#-------------------------------------------

  def check_outputfile(self,outfilename,nruns,reltol=0.1,abstol=1e3):
    status = 'unknown'
    if os.path.isfile(outfilename):
      outf = open(outfilename,'r')
      outlines = outf.read().split('\n')
      finlines = [l for l in outlines if "Optimization finished" in l]
      if len(finlines) < nruns:
        return 'failed' # This function unstable if job was killed.
      displines = [l for l in outlines if "dispersion" in l]
      init_disps = [float(l.split()[4]) for l in displines if "iteration # 1 " in l]
      disps = [float(l.split()[4]) for l in displines]
      if len(disps) > 1:
        dispdiff = abs(disps[-1] - init_disps[-1])/init_disps[-1]
        if (dispdiff < reltol) and (disps[-1] < abstol):
          return 'ok'
        else:
          print("Variance optimization dispersion not converged:")
          print("rel_change(%.3f>%.3f) or abs(%.0f>%.0f)"\
              %(dispdiff,reltol,disps[-1],abstol))
          return 'not_finished'
      else:
        return 'failed'
    else:
      return 'not_started'
#--------------------------------------------------
  def check_status(self,job_record):
    # TODO check different output files (the ones that are requested)
    nruns=job_record['qmc']['variance_optimize']['nruns']
    reltol = job_record['qmc']['variance_optimize']['reltol']
    abstol = job_record['qmc']['variance_optimize']['abstol']
    fnames=[]
    for jast in job_record['qmc']['variance_optimize']['jastrow']:
      fnames.append("qw_0.%s.opt"%jast)

    outfnames=[]
    wfoutnames=[]
    for f in fnames:
      outfnames.append(f+".o")
      wfoutnames.append(f+".wfout")

    #First check if all the runs are ok.
    #If so, we return ok
    all_ok=True
    for outfilename in outfnames:
      if self.check_outputfile(outfilename,nruns,reltol,abstol)!='ok':
        all_ok=False
    if all_ok:
      return 'ok'

    #Check on the submitter. If still running report that.
    status=self._submitter.status(job_record,self._name_)
    if 'running' in status:
      return 'running'
    
    #If not running, try to transfer files.
    print(fnames,outfnames,wfoutnames)
    self._submitter.transfer_output(job_record, fnames+outfnames+wfoutnames)

    #Now check on the output files again
    statuses=[]
    for outfilename in outfnames:
      statuses.append(self.check_outputfile(outfilename,nruns,reltol,abstol))
    
    #Finally, decide what to do
    if len(set(statuses))==1:
      print("all statuses the same")
      return statuses[0]
    if 'not_finished' in statuses:
      return 'not_finished'
    
    #We may have some failed and some not..
    print("Not sure what to do right now..")
    print(statuses)
    quit()
#-------------------------------------------------      
  def resume(self,job_record,maxresume=5):
    #print("resume currently broken")
    #quit()
    infiles=[]
    jastfiles=[]
    for jast in job_record['qmc']['variance_optimize']['jastrow']:
      jast_suf=""
      if jast=='twobody':
        jast_suf = 'jast2'
      elif jast=='threebody':
        jast_suf = 'jast3'
      else:
        print("Didn't understand Jastrow",jast)
        quit()
      fname="qw_0.%s.opt"%jast
      f = open(fname,'w')
      nit = job_record['qmc']['variance_optimize']['niterations']
      nruns = job_record['qmc']['variance_optimize']['nruns']
      for i in range(0,nruns):
        f.write("method { optimize iterations %i } "%nit)
      f.write("\n".join([
          "include qw_0.sys",
          "trialfunc { slater-jastrow",
          "  wf1 { include qw_0.slater } ",
          "  wf2 { include qw.%s } "%jast_suf,
          "}"
        ]))
      infiles.append(fname)
      jastfiles.append("qw.%s"%jast_suf)
      f.close()
    outfiles=[]
    for fname in infiles:
      outfiles.append(fname+".stdout")

    self._submitter.execute(
      job_record, 
      infiles+['qw_0.sys','qw_0.slater','qw_0.orb','qw.basis']+jastfiles, 
      infiles,
      outfiles[0],
      self._name_)
    
    return 'running'
    ########OLD########
    fname="qw_0.%s.opt"%jast
    if not os.path.isfile(fname):
      return self.run(job_record)
    else: # Save previous output.
      trynum=0
      while os.path.isfile("%d.%s.o"%(trynum,fname)):
        trynum += 1
        if trynum > maxresume:
          print("Not resuming because resume limit reached ({}>{}).".format(
            trynum,maxresume))
          return 'failed'
      shutil.move("%s.o"%fname,"%d.%s.o"%(trynum,fname))

    nit=job_record['qmc']['variance_optimize']['niterations']
    nruns=job_record['qmc']['variance_optimize']['nruns']
    inplines = ["method { optimize iterations %i } "%nit for i in range(nruns)]
    inplines += [
        "include qw_0.sys",
        "trialfunc { include %s.wfout }"%fname
      ]
    with open(fname,'w') as inpf:
      inpf.write('\n'.join(inplines))

    job_record['control']['queue_id'] = self._submitter.execute(
      job_record, 
      ['qw_0.opt','qw_0.sys','qw_0.slater','qw_0.orb','qw.basis','qw.%s'%jast], 
      'qw_0.opt',
      'qw_0.opt.stdout',self._name_)
    
    return 'running'
#------------------------------------------------------
  def output(self,job_record):
    for jast in job_record['qmc']['variance_optimize']['jastrow']:
      job_record['qmc']['variance_optimize'][jast]={}
      
      outfilename="qw_0.%s.opt.o"%jast
      f=open(outfilename,'r')
      disp=[]
      for line in f:
        if 'iteration' in line and 'dispersion' in line:
          spl=line.split()
          disp.append(float(spl[4]))
      job_record['qmc']['variance_optimize'][jast]['sigma']=disp
    return job_record

####################################################
class QWalkEnergyOptimize:
  _name_="QWalkEnergyOptimize"
  
  def __init__(self,submitter):
    self._submitter=submitter
  
#-------------------------------------------------      
  def run(self,job_record,restart=False):
    infiles=[]
    jastfiles=[]
    for jast in job_record['qmc']['energy_optimize']['jastrow']:
      jast_suf=""
      if jast=='twobody':
        jast_suf = 'jast2'
      elif jast=='threebody':
        jast_suf = 'jast3'
      else:
        print("Didn't understand Jastrow",jast)
        quit()
    
      fname="qw_0.%s.enopt"%jast
      # TODO make restart work with 2 and 3-body jastrow
      if restart:
        if not os.path.isfile("qw_0.enopt.wfout"):
          print("Could not find qw_0.enopt.wfout")
          return "failed"

        os.system("cp qw_0.enopt.wfout qw_0.enopt.wfin")
      else:
        os.system("sed s/OPTIMIZEBASIS//g qw_0.%s.opt.wfout > %s.wfin"%(jast,fname))

      enopt_options=job_record['qmc']['energy_optimize']
      f=open(fname,'w')
      f.write("""method { LINEAR VMC_NSTEP %i } 
include qw_0.sys
trialfunc { include %s.wfin }
"""%(enopt_options['vmc_nstep'],fname))
      infiles.append(fname)
      jastfiles.append("qw.%s"%jast_suf)
      f.close()
    outfiles=[]
    for fname in infiles:
      outfiles.append(fname+".stdout")
    self._submitter.execute(
      job_record, 
      infiles+['%s.wfin'%f for f in infiles]+['qw_0.sys','qw_0.slater','qw_0.orb','qw.basis'],
      infiles,
      outfiles[0],
      self._name_)
    
    return 'running'


#-------------------------------------------------      
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


#-------------------------------------------------      
  def check_status(self,job_record):
    thresh=job_record['qmc']['energy_optimize']['threshold']
    fnames=[]
    for jast in job_record['qmc']['variance_optimize']['jastrow']:
      fnames.append("qw_0.%s.enopt"%jast)

    outfnames=[]
    wfoutnames=[]
    for f in fnames:
      outfnames.append(f+".o")
      wfoutnames.append(f+".wfout")

    #Check on the submitter. If still running report that.
    status=self._submitter.status(job_record,self._name_)
    if 'running' in status:
      return 'running'
    
    #If not running, try to transfer files.
    #print(fnames,outfnames,wfoutnames)
    self._submitter.transfer_output(job_record, fnames+outfnames+wfoutnames)

    #Now check on the output files again
    statuses=[]
    for outfilename in outfnames:
      statuses.append(self.check_outputfile(outfilename,thresh))
    
    #Finally, decide what to do
    if len(set(statuses))==1:
      return statuses[0]
    if 'not_finished' in statuses:
      return 'not_finished'
    #We may have some failed and some not..
    print("Not sure what to do right now..")
    print(statuses)
    quit()

      
#-------------------------------------------------      
  def resume(self,job_record):
    return self.run(job_record,restart=True)

#-------------------------------------------------      
  def output(self,job_record):
    
    for jast in job_record['qmc']['energy_optimize']['jastrow']:
      job_record['qmc']['energy_optimize'][jast]={}
      outfilename="qw_0.%s.enopt.o"%jast
      f=open(outfilename,'r')
      energy=[]
      energy_err=[]
      for line in f:
        if 'current energy' in line:
          spl=line.split()
          if len(spl) > 9:
            energy.append(float(spl[4]))
            energy_err.append(float(spl[6]))
      job_record['qmc']['energy_optimize'][jast]['energy']=energy
      job_record['qmc']['energy_optimize'][jast]['energy_err']=energy_err
    return job_record

####################################################

class QWalkRunVMC:
  _name_="QwalkRunVMC"
  
  #-----------------------------------------------
  def __init__(self,submitter):
    self._submitter=submitter
  #-----------------------------------------------
  def run(self,job_record,restart=False):
    qmc_options=job_record['qmc']
    kpts=self.get_kpts(job_record)
    #if self._name_+'_jobid' not in job_record['control'].keys():
    #  job_record['control'][self._name_+'_jobid'] = []
    #  job_record['control']['queue_id'] = []
    
    #choose which wave function to use
    if not restart:
      if qmc_options['vmc']['optimizer']=='variance':
        os.system("separate_jastrow qw_0.opt.wfout > opt.jast")
      elif qmc_options['vmc']['optimizer']=='energy':
        os.system("separate_jastrow qw_0.enopt.wfout > opt.jast")
      elif qmc_options['vmc']['optimizer']==None:
        shutil.copy("qw.jast2","opt.jast")
      else:
        print("Warning: didn't understand optimizer option!")

    # Make and submit the runs: bundle all jobs.
    infiles = []
    inpfns = []
    for k in kpts:
      kname="qw_%i"%k
      basename="qw_%i"%k
      f=open(basename+".vmc",'w')
      f.write(self.vmcinput(k,qmc_options['vmc']['nblock']))
      f.close()
      infiles.extend([basename+".vmc","opt.jast",kname+'.sys',kname+'.slater',
          kname+'.orb','qw.basis'])
      if restart:
        infiles.extend([basename+'.vmc.config',basename+'.vmc.log'])
      inpfns.append(basename+".vmc")

    qid = self._submitter.execute(
      job_record,
      infiles, # Dependencies. This naming needs to be less redundant.
      inpfns,  # Actual DMC inputs.
      "qw.vmc.stdout")
    return 'running'
  #-----------------------------------------------
  def get_kpts(self, job_record):
    kpts=glob.glob("qw*.sys")
    kpt_num=[]
    for kp in kpts:
      kpt_num.append(int(re.findall(r'\d+',kp)[0]))
    return kpt_num
  #-----------------------------------------------
  def vmcinput(self,kpt_num=0,nblock=16):
    outlist = [
        "method { VMC ",
        "nblock %i"%nblock
      ]
    outlist += [
        "}",
        "include qw_%i.sys"%kpt_num,
        "trialfunc { slater-jastrow",
        "wf1 { include qw_%i.slater } "%kpt_num,
        "wf2 { include opt.jast }",
        "}"
      ]
    outstr = '\n'.join(outlist)
    return outstr
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
  def collect_runs(self,job_record):
    ret=[]

    kpts=self.get_kpts(job_record)
    for k in kpts:
      kname="qw_%i"%k
      basename="qw_%i"%k
      entry={}
      entry['knum']=k
      entry['energy']=self.energy(basename+".vmc.log")
      ret.append(entry)
    return ret
  #-----------------------------------------------
  def check_status(self,job_record):
    results=self.collect_runs(job_record)
    
    status='ok'
    for e in results:
      print("energy check",e['knum'],e['energy'])
        
      if e['energy'][1] >  job_record['qmc']['vmc']['target_error']:
        status='not_finished'
    print("initial status",status)
    if status=='ok':
      return status

    outfiles=[]
    for k in self.get_kpts(job_record):
      basename="qw_%i"%k
      outfiles.extend([basename+'.vmc.log',
                     basename+'.vmc.config',
                     basename+'.vmc.o'])
    #print(outfiles)
    self._submitter.transfer_output(job_record, outfiles)
    status=self._submitter.status(job_record)
    print("status",status)
    if status == 'running':
      return status

    if not os.path.isfile(outfiles[0]):
      return 'not_started'

    results=self.collect_runs(job_record)
    status='ok'
    for e in results:
      if e['energy'][1] >  job_record['qmc']['vmc']['target_error']:
        status='not_finished'
    return status
  #-----------------------------------------------
  def resume(self,job_record):
    return self.run(job_record,restart=True)
  #-----------------------------------------------
  def output(self,job_record):
    job_record['qmc']['vmc']['results']=self.collect_runs(job_record)
    return job_record

####################################################

class QWalkRunDMC:
  _name_="QwalkRunDMC"
  
  def __init__(self,submitter):
    self._submitter=submitter
#-----------------------------------------------
  def run(self,job_record,restart=False):
    options=job_record['qmc']['dmc']
    kpts=self.get_kpts(job_record)
    
    calc_sk=False
    if 'cif' in job_record.keys():
      calc_sk=True
    # Make and submit the runs: bundle all jobs.
    depfns = []# Dependencies. 
    inpfns = [] #DMC inputs
    for k in kpts:
      for t in options['timestep']:
        for loc in options['localization']:
          for jast in options['jastrow']:
            for opt in options['optimizer']:
              kname="qw_%i"%k
              basename=self.gen_basename(k,t,loc,jast,opt)
              f=open(basename+".dmc",'w')
              f.write(self.dmcinput(k,t,loc,jast,opt,
                options['nblock'],
                options['save_trace'],calc_sk))
              f.close()

#Warning: remote may not be working with this..
              depfns.extend([basename+".dmc",
                             "opt.jast",
                             kname+'.sys',
                             kname+'.slater',
                             kname+'.orb',
                             'qw.basis'])
              if restart:
                depfns.extend([basename+'.dmc.config',basename+'.dmc.log'])
              inpfns.append(basename+".dmc")

    self._submitter.execute(
      job_record,
      depfns,
      inpfns,  # Actual DMC inputs.
      "qw.dmc.stdout",
      self._name_)
    return 'running'

#-----------------------------------------------
  def gen_basename(self,k,t,loc,jast,opt):
    return "qw_%i_%s_%g_%s_%s"%(k,jast,t,opt,loc)

#-----------------------------------------------
  def get_kpts(self, job_record):
    kpts=glob.glob("qw*.sys")
    kpt_num=[]
    for kp in kpts:
      kpt_num.append(int(re.findall(r'\d+',kp)[0]))
    return kpt_num

#-----------------------------------------------
  def dmcinput(self,k,t,loc,jast,opt,nblock=16,save_trace=False,sk=False):
    basename=self.gen_basename(k,t,loc,jast,opt)
    outlist = [
        "method { DMC ",
        "timestep %g"%t,
        "nblock %i"%nblock,
        loc
      ]
    if sk:
      outlist.append("average { SK } ")
    if save_trace:
      outlist += ["save_trace %s.trace"%basename]
    opt_trans={"energy":"enopt","variance":"opt"}
    outlist += [
        "}",
        "include qw_%i.sys"%k,
        "trialfunc { include qw_0.%s.%s.wfout"%(jast,opt_trans[opt]),
        "}"
      ]
    outstr = '\n'.join(outlist)
    return outstr

#-----------------------------------------------
  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if 'Wall' in line:
          return 'ok'
      return 'running'

#-----------------------------------------------
  def collect_runs(self,job_record):
    ret=[]
    options=job_record['qmc']['dmc']
    kpts=self.get_kpts(job_record)
    for k in kpts:
      for t in options['timestep']:
        for loc in options['localization']:
          for jast in options['jastrow']:
            for opt in options['optimizer']:
              basename=self.gen_basename(k,t,loc,jast,opt)
              if os.path.isfile("%s.dmc.log"%basename):
                entry={}
                entry['knum']=k
                entry['timestep']=t
                entry['localization']=loc
                entry['jastrow']=jast
                entry['optimizer']=opt
                os.system("gosling -json %s.dmc.log > %s.json"%(basename,basename))
                entry['results']=json.load(open("%s.json"%basename))
                ret.append(entry)
    return ret

#-----------------------------------------------
  def check_status(self,job_record):
    
    options=job_record['qmc']['dmc']
    kpts=self.get_kpts(job_record)
    infns = [] #DMC inputs
    for k in kpts:
      for t in options['timestep']:
        for loc in options['localization']:
          for jast in options['jastrow']:
            for opt in options['optimizer']:
              infns.append(self.gen_basename(k,t,loc,jast,opt))
    
    #Check on the submitter. If still running report that.
    status=self._submitter.status(job_record,self._name_)
    if 'running' in status:
      return 'running'
    
    #If not running, try to transfer files.
    self._submitter.transfer_output(job_record, infns)

    #Now check on the runs
    ret=self.collect_runs(job_record)
    if len(ret)==0:
      return "not_started"
    if len(ret) != len(infns):
      print("There are no jobs running and not enough .log files. Not sure what's going on.")
      quit()
    
    statuses=[]
    thresh=options['target_error']
    for r in ret:
      if r['results']['properties']['total_energy']['error'][0] < thresh:
        statuses.append("ok")
      else:
        status.append("not_finished")
    #Finally, decide what to do
    if len(set(statuses))==1:
      return statuses[0]
    if 'not_finished' in statuses:
      return 'not_finished'
    #We may have some failed and some not..
    print("Not sure what to do right now..")
    print(statuses)
    quit()
    
#-----------------------------------------------
    
      
  def resume(self,job_record):
    return self.run(job_record,restart=True)
#-----------------------------------------------

  def output(self,job_record):
    job_record['qmc']['dmc']['results']=self.collect_runs(job_record)
    return job_record

####################################################

class QWalkRunMaximize:
  _name_ = "QWalkRunMaximize"

  def __init__(self,submitter):
    self._submitter = submitter

  def make_basename(self, k, n, w, s):
    return "qw_%i.%s.n%i"%(k,w,n)

  def make_kname(self, k, s):
    return "qw_%i"%k

#-----------------------------------------------
  def run(self, job_record, restart=False):
    qmc_options = job_record['qmc']
    nconfiglist = qmc_options['maximize']['nconfig']
    w = qmc_options['maximize']['trialwf']
    s = qmc_options['maximize']['system']

    #choose which wave function to use
    if not restart:
      if qmc_options['dmc']['optimizer']=='variance':
        os.system("separate_jastrow qw_0.opt.wfout > opt.jast")
      elif qmc_options['dmc']['optimizer']=='energy':
        os.system("separate_jastrow qw_0.enopt.wfout > opt.jast")

    depfns = []
    inpfns = []
    kpts=self.get_kpts(job_record)

    for k in kpts:
      for n in nconfiglist:
          kname = self.make_kname(k, s)
          basename = self.make_basename(k, n, w, s)
          f = open(basename+".max",'w')
          f.write(self.maximizeinput(k, n, w, s))
          f.close()
          depfns.extend([basename+".max","opt.jast",kname+'.sys',kname+'.slater',kname+'.orb','qw.basis'])
          if restart:
            depfns.extend([basename+".max.config",basename+".max.log"])
          inpfns.append(basename+".max")

    qid = self._submitter.execute(
      job_record,
      depfns, # Dependencies. This naming needs to be less redundant
      inpfns, # Actual MAXIMIZE inputs
      "qw.max.sdout")

    return 'running'

#-----------------------------------------------
  def get_kpts(self, job_record):
    kpts=glob.glob("qw*.sys")
    kpt_num=[]
    for kp in kpts:
      kpt_num.append(int(re.findall(r'\d+',kp)[0]))
    return kpt_num

#-----------------------------------------------
  def maximizeinput(self, k, nconfig, wf, sys):
    outlist = [
      "method { MAXIMIZE "+
      "NCONFIG %i "%nconfig+
      "}"]
    outlist.append("include %s.sys"%self.make_kname(k,sys))

    if wf == "hf" or wf == "slater":
      outlist.append("trialfunc { include %s.slater }"%self.make_kname(k,sys))
    else:
      #outlist.append("trialfunc { include qw_0.opt.wfout }")
      outlist += [
        "trialfunc { slater-jastrow",
        "wf1 { include %s.slater } "%self.make_kname(k,sys),
        "wf2 { include opt.jast }",
        "}"
      ]

    outstr = '\n'.join(outlist)
    return outstr

#-----------------------------------------------
  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if 'Wall' in line:
          return 'ok'
      return 'running'

#-----------------------------------------------
  def extract_data(self,logfilename):
    data = np.loadtxt(logfilename,skiprows=1)
    amp = data[:,-1]
    locE = data[:,-2]
    conf = data[:,:-2]
    return amp, locE, conf

#-----------------------------------------------
  def collect_runs(self,job_record):
    ret=[]
    w = job_record['qmc']['maximize']['trialwf']
    kpts=self.get_kpts(job_record)

    for k in kpts:
      for n in job_record['qmc']['maximize']['nconfig']:
          logfilename = self.make_basename(k, n, w, s)+".max.table"
          entry = {}
          entry['nconfig'] = n
          amp, locE, conf = self.extract_data(logfilename)
          entry['psi'] = amp.tolist()
          entry['energies'] = locE.tolist()
          entry['configs'] = conf.tolist()
          ret.append(entry)
    return ret

#-----------------------------------------------
  def check_status(self,job_record):
    # TODO does this function do what it's supposed to?
    #results=self.collect_runs(job_record)

    status='ok'

    # TODO do we need something here?
    #for e in results:
      #print("energy check",e['knum'],e['energy'])

      #if e['energy'][1] >  job_record['qmc']['dmc']['target_error']:
      #  status='not_finished'
    #print("initial status",status)
    #if status=='ok':
    #  return 'ok'

    w = job_record['qmc']['maximize']['trialwf']
    s = job_record['qmc']['maximize']['system']
    kpts=self.get_kpts(job_record)
    outfiles=[]
    for k in kpts:
      for n in job_record['qmc']['maximize']['nconfig']:
          basename = self.make_basename(k, n, w, s)
          outfiles.extend([basename+".max.table",
                          basename+".max.log",
                          basename+".max.config",
                          basename+".max.o"])

    self._submitter.transfer_output(job_record, outfiles)
    status=self._submitter.status(job_record)
    print("status",status)
    if status=='running':
      return status


    if not os.path.isfile(outfiles[0]):
      return 'not_started'

    results=self.collect_runs(job_record)
    status='ok'
    # TODO do we need something here?
    #for e in results:
      #if e['energy'][1] >  job_record['qmc']['dmc']['target_error']:
      #  status='not_finished'
    return status

#-----------------------------------------------
  def retry(self,job_record):
    return self.run(job_record,restart=True)


#-----------------------------------------------
  def output(self,job_record):
    job_record['qmc']['maximize']['results']=self.collect_runs(job_record)
    return job_record

####################################################

class QWalkRunPostProcess:
  _name_="QwalkRunPostProcess"
  
  #-----------------------------------------------
  def __init__(self,submitter):
    self._submitter=submitter

  #-----------------------------------------------
  def run(self,job_record,restart=False):
    options=job_record['qmc']['dmc']
    kpts=self.get_kpts(job_record)
    
    calc_sk=False
    if 'cif' in job_record.keys():
      calc_sk=True
    # Make and submit the runs: bundle all jobs.
    depfns = [] # Dependencies. 
    inpfns = [] # DMC inputs
    for k in kpts:
      for t in options['timestep']:
        for loc in options['localization']:
          for jast in options['jastrow']:
            for opt in options['optimizer']:
              kname="qw_%i"%k
              basename=self.gen_basename(k,t,loc,jast,opt)
              f=open(basename+".post",'w')
              f.write(self.postprocessinput(k,t,loc,jast,opt,
                job_record['qmc']['postprocess']))
              f.close()

#Warning: remote may not be working with this..
              depfns.extend([basename+".post",
                             "opt.jast",
                             kname+'.sys',
                             kname+'.slater',
                             kname+'.orb',
                             'qw.basis'])
              inpfns.append(basename+".post")

    self._submitter.execute(
      job_record,
      depfns, inpfns,
      "qw.post.stdout",
      self._name_)
    return 'running'

  #-----------------------------------------------
  def gen_basename(self,k,t,loc,jast,opt):
    return "qw_%i_%s_%g_%s_%s"%(k,jast,t,opt,loc)

  #-----------------------------------------------
  def get_kpts(self, job_record):
    kpts=glob.glob("qw*.sys")
    kpt_num=[]
    for kp in kpts:
      kpt_num.append(int(re.findall(r'\d+',kp)[0]))
    return kpt_num

  #-----------------------------------------------
  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if 'Wall' in line:
          return 'ok'
      return 'running'

  #-----------------------------------------------
  def get_warmup(self,logfilename):
    os.system("gosling %s > %s.stdout"%(logfilename,logfilename))
    f=open(logfilename+".stdout",'r')
    nwarm=0
    for line in f:
      if "Threw out" in line:
        spl=line.split()
        nwarm=int(spl[4])
        return nwarm
    return nwarm

  #-----------------------------------------------
  def postprocessinput(self,k,t,loc,jast,opt,ppr_options):
    basename=self.gen_basename(k,t,loc,jast,opt)
    nwarmup = self.get_warmup("%s.dmc.log"%basename)
    print("nwarmup:",nwarmup)
    outlines = [
        "method { postprocess ",
        "readconfig %s.trace"%basename,
        "noenergy",
        "average { region_fluctuation }"
      ]
    nskip = nwarmup*2048 # TODO generalize 2048.
    if ppr_options['region_fluctuation'] == True:
      outlines += ["average { region_fluctuation }"]
    if ppr_options['density'] == True:
      outlines += [
          "density { density up   outputfile %s.dmc.up.cube }"%basename,
          "density { density down outputfile %s.dmc.dn.cube }"%basename
        ]
    if ppr_options['obdm'] == True:
      if (ppr_options['basis'] != None) and (ppr_options['orb'] != None):
        shutil.copy(ppr_options['basis'],ppr_options['basis'].replace("../",""))
        shutil.copy(ppr_options['orb'],ppr_options['orb'].replace("../",""))
        outlines += [
            "average { tbdm_basis ",
            "mode obdm",
            "include atomic.basis", #TODO generalize naming.
            "}"
          ]
      else:
        print("Since min basis wasn't defined, postprocess can't do OBDM.")
    opt_trans={"energy":"enopt","variance":"opt"}
    outlist += [
        "}",
        "include qw_%i.sys"%k,
        "trialfunc { include qw_0.%s.%s.wfout"%(jast,opt_trans[opt]),
        "}"
      ]
    outstr = '\n'.join(outlines)
    return outstr

  #-----------------------------------------------
  def collect_runs(self,job_record):
    ret=[]
    options=job_record['qmc']['dmc']
    kpts=self.get_kpts(job_record)
    for k in kpts:
      for t in options['timestep']:
        for loc in options['localization']:
          for jast in options['jastrow']:
            for opt in options['optimizer']:
              basename=self.gen_basename(k,t,loc,jast,opt)
              print("Debug basename",basename)
              if os.path.isfile("%s.post.o"%basename):
                entry={}
                with open(basename+".post.o") as postin:
                  nfres = list(map(lambda x:x.tolist(),
                      nf.read_number_dens(postin)))
                entry['number_fluctuation']=nfres
                with open(basename+".post.o") as postin:
                  nmres = dm.read_dm(postin)
                for key in nmres.keys(): nmres[key] = nmres[key].tolist()
              else:
                entry['number_fluctuation']=(None,None)
                entry['ordm']=(None,None)
                entry['tbdm']=(None,None)
                entry['density_matrices'] = nmres
                entry['timestep']=t
                entry['localization']=loc
                ret.append(entry)
    return ret

  #-----------------------------------------------
  def check_status(self,job_record):
    options=job_record['qmc']['dmc']
    kpts=self.get_kpts(job_record)
    infns = [] #DMC inputs
    for k in kpts:
      for t in options['timestep']:
        for loc in options['localization']:
          for jast in options['jastrow']:
            for opt in options['optimizer']:
              infns.append(self.gen_basename(k,t,loc,jast,opt))
    
    status='ok'
    for e in results:
      print("Postprocess",e['knum'],e['number_fluctuation'] != (None,None))
      print(e['number_fluctuation'])
      if e['number_fluctuation'][0] is None:
        status='failed'
    print("initial status",status)
    if status=='ok':
      return status

    outfiles=[]
    for k in self.get_kpts(job_record):
      for t in job_record['qmc']['dmc']['timestep']:
        for local in job_record['qmc']['dmc']['localization']:
          basename="qwt%g%s_%i"%(t,local,k)
          outfiles.extend([basename+'.post.o'])
    self._submitter.transfer_output(job_record, outfiles)
    status=self._submitter.status(job_record)
    print("status",status)
    if status=='running':
      return status

    #Now check on the runs
    ret=self.collect_runs(job_record)
    if len(ret)==0:
      return "not_started"
    if len(ret) != len(infns):
      print("There are no jobs running and not enough .log files. Not sure what's going on.")
      quit()
    
    statuses=[]
    thresh=options['target_error']
    for r in ret:
      if 'number_fluctuation' in r.keys():
        statuses.append("ok")
      else:
        status.append("not_finished")
    #Finally, decide what to do
    if len(set(statuses))==1:
      return statuses[0]
    if 'not_finished' in statuses:
      return 'not_finished'
    #We may have some failed and some not..
    print("Not sure what to do right now..")
    print(statuses)
    quit()

  #-----------------------------------------------
  def resume(self,job_record):
    return self.run(job_record,restart=True)

  #-----------------------------------------------
  def output(self,job_record):
    job_record['qmc']['postprocess']['results']=self.collect_runs(job_record)
    return job_record

####################################################

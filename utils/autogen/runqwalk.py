from __future__ import print_function
import sys
import os
import glob
import re
import shutil
from crystal2qmc import convert_crystal
import subprocess as sub
import numpy as np
import json
import yaml
import time

# If you need the swap_endian option, you need to set this to the correct location.
swap_endian_exe = "/home/busemey2/bin/swap_endian"

####################################################

def extract_jastrow(f):
  tokens=f.readlines()
  jastrow=""
  in_jastrow=False
  nopen=0
  nclose=0
  for line in tokens:
    if line.find("JASTROW2") != -1:
      in_jastrow=True
    if in_jastrow:
      nopen+=line.count("{")
      nclose+=line.count("}")
    if in_jastrow and nopen >= nclose:
      jastrow+=line
  return jastrow
  
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
    if job_record['qmc']['kpoints']=='real':
      os.system("crystal2qmc -o qw patched.o > crystal2qmc.stdout")
    elif job_record['qmc']['kpoints']=='all':
      os.system("crystal2qmc -c -o qw patched.o > crystal2qmc.stdout")
    else:
      print('Error in kpoints input.')
      quit()
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
      displines = [l for l in outlines if "dispersion" in l]
      if len(displines) < 4:
        print("Only completed four optimization routines. May want to increase the queue time.")
        return "not_finished" 
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

    #Check on the submitter. If still running report that.
    status=self._submitter.status(job_record,self._name_)
    if 'running' in status:
      return 'running'
    
    #If not running, try to transfer files.
    print(fnames,outfnames,wfoutnames)
    self._submitter.transfer_output(job_record, fnames+outfnames+wfoutnames)

    #Now check on the output files 
    statuses=[]
    for outfilename in outfnames:
      statuses.append(self.check_outputfile(outfilename,nruns,reltol,abstol))
    
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
  def resume(self,job_record,maxresume=5):
    infiles=[]
    for jast in job_record['qmc']['variance_optimize']['jastrow']:
      infiles.append("qw_0.%s.opt"%jast)
    
    for inf in infiles:
      if os.path.isfile(inf+".wfout"): #save previous output
        trynum=0
        while os.path.isfile("%d.qw_0.opt.o"%trynum):
          trynum += 1
          if trynum > maxresume:
            print("Not resuming because resume limit reached ({}>{}).".format(
              trynum,maxresume))
            return 'failed'
        shutil.move(inf+".o","%d."%trynum + inf+".o")
      else: 
        print("VarianceOptimize: asked to resume a job which didn't run")
        quit()

      nit=job_record['qmc']['variance_optimize']['niterations']
      nruns=job_record['qmc']['variance_optimize']['nruns']
      inplines = ["method { optimize iterations %i } "%nit for i in range(nruns)]
      inplines += [
          "include qw_0.sys",
          "trialfunc { include %s.wfout }"%inf
        ]
      with open(inf,'w') as inpf:
        inpf.write('\n'.join(inplines))

    wffiles=[]
    for inf in infiles:
      wffiles.append(inf+".wfout")
    self._submitter.execute(job_record, 
        ['qw_0.sys','qw_0.slater','qw_0.orb','qw.basis']+infiles+wffiles, 
         infiles,infiles[0]+'.stdout',self._name_)
    
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
      if restart:
        if not os.path.isfile("%s.wfout"%fname):
          print("Could not find %s.wfout"%fname)
          return "failed"

        os.system("cp %s.wfout %s.wfin"%(fname,fname))
        #Remove the .config file because sometimes the number of
        #processors changes.
        try:
          os.remove("%s.config"%fname)
        except:
          pass
      else:
        os.system("sed s/OPTIMIZEBASIS//g qw_0.%s.opt.wfout > %s.wfin"%(jast,fname))

      enopt_options=job_record['qmc']['energy_optimize']
      f=open(fname,'w')
      f.write("""method { LINEAR TOTAL_NSTEP %i } 
include qw_0.sys
trialfunc { include %s.wfin }
"""%(enopt_options['total_nstep'],fname))
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
  _name_ = "QWalkRunVMC"

  def __init__(self,submitter):
    self._submitter = submitter

#-----------------------------------------------
  def gen_basename(self,k,jast,opt):
    return "qw_%i_%s_%s"%(k,jast,opt)

#-----------------------------------------------
  def run(self, job_record, restart=False):
    options=job_record['qmc']['vmc']
    kpts=self.get_kpts(job_record)
    
    depfns = [] # Dependencies. 
    inpfns = [] # DMC inputs
    for k in kpts:
      for jast in options['jastrow']:
        for opt in options['optimizer']:
          kname="qw_%i"%k
          basename=self.gen_basename(k,jast,opt)
          f=open(basename+".vmc",'w')
          f.write(self.vmcinput(k,jast,opt,options['nblock']))
          f.close()

          depfns.extend([basename+".dmc",
                         "opt.jast",
                         kname+'.sys',
                         kname+'.slater',
                         kname+'.orb',
                         'qw.basis'])
          if restart:
            depfns.extend([basename+'.vmc.config',basename+'.vmc.log'])
          inpfns.append(basename+".vmc")

    self._submitter.execute(
      job_record,
      depfns,
      inpfns,  # Actual VMC inputs.
      "qw.vmc.stdout",
      self._name_)
    return 'running'

#-----------------------------------------------
  def get_kpts(self, job_record):
    kpts=glob.glob("qw*.sys")
    kpt_num=[]
    for kp in kpts:
      kpt_num.append(int(re.findall(r'\d+',kp)[0]))
    return kpt_num

#-----------------------------------------------
  def vmcinput(self,k,jast,opt,nblock):
    #outlist = [
    #    "method { VMC ",
    #    "nblock %i"%nblock,
    #  ]
    #opt_trans={"energy":"enopt","variance":"opt"}
    #jast_inp=extract_jastrow(open("qw_0.%s.%s.wfout"%(jast,opt_trans[opt])))
    #outlist += [
    #    "}",
    #    "include qw_%i.sys"%k,
    #    "trialfunc { slater-jastrow ",
    #       "wf1 { include qw_%i.slater } "%(k),
    #       "wf2 { ",
    #      jast_inp,
    #    " } ",
    #    "}"
    #  ]
    #outstr = '\n'.join(outlist)
    #return outstr
    outlist = [
      "method { VMC nblock %i }"%nblock,
      "include qw_%i.sys"%k
      ]
    if jast=="none" or jast=='':
      outlist.append("trialfunc { include qw_%i.slater }"%k)
    else:
      opt_trans={"energy":"enopt","variance":"opt"}
      jast_inp=extract_jastrow(open("qw_0.%s.%s.wfout"%(jast,opt_trans[opt])))
      outlist += [
        "trialfunc { slater-jastrow ",
        "  wf1 { include qw_%i.slater }"%(k),
        "  wf2 { ", jast_inp, " } ",
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
    options=job_record['qmc']['vmc']
    kpts=self.get_kpts(job_record)
    for k in kpts:
      for jast in options['jastrow']:
        for opt in options['optimizer']:
          basename=self.gen_basename(k,jast,opt)
          if os.path.isfile("%s.vmc.log"%basename):
            entry={}
            entry['knum']=k
            entry['jastrow']=jast
            entry['optimizer']=opt
            os.system("gosling -json %s.vmc.log > %s.json"%(basename,basename))
            entry['results']=json.load(open("%s.json"%basename))
            ret.append(entry)
    return ret

#-----------------------------------------------
  def check_status(self,job_record):

    options=job_record['qmc']['vmc']
    kpts=self.get_kpts(job_record)
    infns = [] #VMC inputs
    for k in kpts:
      for jast in options['jastrow']:
        for opt in options['optimizer']:
          infns.append(self.gen_basename(k,jast,opt))
    
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
        statuses.append("not_finished")
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
  def retry(self,job_record):
    return self.run(job_record,restart=True)

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
    if restart:
      ret = self.collect_runs(job_record)
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
              dep=[basename+".dmc",
                             "opt.jast",
                             kname+'.sys',
                             kname+'.slater',
                             kname+'.orb',
                             'qw.basis']
              if restart:
                results = None
                for r in ret:
                  if r['knum'] == k and r['timestep']==t and r['localization']==loc\
                    and r['jastrow']==jast and r['optimizer']==opt:
                    results = r['results']
                thresh = job_record['qmc']['dmc']['target_error']
                if results == None or results['properties']['total_energy']['error'][0] >= thresh:
                  print('%s not finished '%(basename))
                  depfns.extend(dep)
                  depfns.extend([basename+'.dmc.config',basename+'.dmc.log'])
                  inpfns.append(basename+".dmc")
              else:
                depfns.extend(dep)
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
    jast_inp=extract_jastrow(open("qw_0.%s.%s.wfout"%(jast,opt_trans[opt])))
    outlist += [
        "}",
        "include qw_%i.sys"%k,
        "trialfunc { slater-jastrow ",
           "wf1 { include qw_%i.slater } "%(k),
           "wf2 { ",
          jast_inp,
          " } ",
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
              entry={}
              entry['knum']=k
              entry['timestep']=t
              entry['localization']=loc
              entry['jastrow']=jast
              entry['optimizer']=opt
              entry['results'] = None
              if os.path.isfile("%s.dmc.log"%basename):
                try:
                  os.system("gosling -json %s.dmc.log > %s.json"%(basename,basename))
                  entry['results']=json.load(open("%s.json"%basename))
                except:
                  print("trouble processing",basename)
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
      print("WARNING: There are missing log files. Expected",len(infns), "found", len(ret) )
      return "not_finished"
    
    statuses=[]
    thresh=options['target_error']
    for r in ret:
      if r['results'] == None:
        basename = self.gen_basename(r['knum'],r['timestep'],r['localization'],r['jastrow'],r['optimizer'])
        print('DMC Calculation not finished: %s'%basename)
        statuses.append('not_finished')
      elif r['results']['properties']['total_energy']['error'][0] < thresh:
        statuses.append("ok")
      else:
        print("Stochastic error too large: {0:.2e}>{1:.2e}".format(
            r['results']['properties']['total_energy']['error'][0],
            thresh
          ))
        statuses.append("not_finished")
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

  #def make_basename(self, k, n, w):
  #  return "qw_%i.%s.n%i"%(k,w,n)
  
  def gen_basename(self,k,n,jast,opt="variance"):
    return "qw_%i_%s_%i_%s"%(k,jast,n,opt)

  def make_kname(self, k):
    return "qw_%i"%k

#-----------------------------------------------
  def run(self, job_record, restart=False):
    options=job_record['qmc']['maximize']
    kpts=self.get_kpts(job_record)

    #choose which wave function to use
    #if not restart:
    #  if qmc_options['dmc']['optimizer']=='variance':
    #    os.system("separate_jastrow qw_0.opt.wfout > opt.jast")
    #  elif qmc_options['dmc']['optimizer']=='energy':
    #    os.system("separate_jastrow qw_0.enopt.wfout > opt.jast")

    infiles = [] # Dependencies.
    inpfns = [] # MAXIMIZE inputs

    for k in kpts:
      for n in options['nconfig']:
        for jast in options['jastrow']:
          for opt in options['optimizer']:
            kname="qw_%i"%k
            basename=self.gen_basename(k,n,jast,opt)
            f = open(basename+".max",'w')
            f.write(self.maximizeinput(k,n,jast,opt))
            f.close()
            infiles.extend([basename+".max",
                            "opt.jast",
                            kname+'.sys',
                            kname+'.slater',
                            kname+'.orb',
                            'qw.basis'])
            if restart:
              infiles.extend([basename+".max.config",basename+".max.log"])
            inpfns.append(basename+".max")

    self._submitter.execute(
      job_record,
      infiles, # Dependencies. 
      inpfns, # MAXIMIZE inputs
      "qw.max.sdout",
      self._name_)

    return 'running'

#-----------------------------------------------
  def get_kpts(self, job_record):
    kpts=glob.glob("qw*.sys")
    kpt_num=[]
    for kp in kpts:
      kpt_num.append(int(re.findall(r'\d+',kp)[0]))
    #return kpt_num
    return [0] # only run k=0 for now

#-----------------------------------------------
  def maximizeinput(self, k, n, jast, opt):
    outlist = [
      "method { MAXIMIZE "+
      "NCONFIG %i "%n+
      "}"
      ]
    outlist.append("include qw_%i.sys"%k)

    if jast=="none" or jast=='':
      outlist.append("trialfunc { include qw_%i.slater }"%k)
    else:
      opt_trans={"energy":"enopt","variance":"opt"}
      outlist.append("trialfunc { include qw_%i.%s.%s.wfout }"%(k,jast,opt_trans[opt] ))

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
    amp = []
    locE = []
    conf = []
    if os.path.isfile(logfilename):
      data = np.loadtxt(logfilename,skiprows=1)
      amp = data[:,-1].tolist()
      locE = data[:,-2].tolist()
      conf = data[:,:-2].tolist()
      
    return amp, locE, conf

#-----------------------------------------------
  def collect_runs(self,job_record):
    print_progress = False
    if print_progress: 
      t = time.time()
      print("Time",0,"Starting timer for maximize collect_runs method")

    ret=[]
    options=job_record['qmc']['maximize']
    kpts=self.get_kpts(job_record)
    
    for k in kpts:
      for n in options['nconfig']:
        for jast in options['jastrow']:
          for opt in options['optimizer']:
            logfilename=self.gen_basename(k,n,jast,opt)+".max.json"
            if os.path.isfile(logfilename):
              if print_progress: print("Time",time.time()-t,'k',k,'nconfig',n,'jastrow',jast,'optimizer',opt)
              entry = {}
              entry['knum']=k
              entry['nconfig']=n
              entry['jastrow']=jast
              entry['optimizer']=opt
              #amp, locE, conf = self.extract_data(logfilename)
              jsonfile = open(logfilename,'r')
              if print_progress: print("Time",time.time()-t,"JSON file opened, loading dict...")
              data = json.load(jsonfile)
              if print_progress: print("Time",time.time()-t,"JSON data loaded")

              entry['psi'] = [r['psi'] for r in data]
              if print_progress: print("Time",time.time()-t,"Extracted psi list")
              entry['error'] = [r['error'] for r in data]
              if print_progress: print("Time",time.time()-t,"Extracted error list")
              entry['energies'] = [r['energy'] for r in data]
              if print_progress: print("Time",time.time()-t,"Extracted energy list")
              entry['configs'] = [r['config'] for r in data]
              if print_progress: print("Time",time.time()-t,"Extracted config list")
              #entry['hessian'] = [r['hessian'] for r in data]
              #if print_progress: print("Time",time.time()-t,"Extracted hessian list")
              
              #l = [[r['psi'],r['energy'],r['config'],r['hessian']] for r in data]
              #if print_progress: print("Time",time.time()-t,"Extracted list of all data, need to 'transpose'")
              #[entry['psi'],entry['energies'],entry['configs'],entry['hessian']] = list(map(list,zip(*l)))
              #if print_progress: print("Time",time.time()-t,"Reshaped into lists for each quantity")
               
              ret.append(entry)
    return ret

#-----------------------------------------------
  def check_status(self,job_record):
    
    options=job_record['qmc']['maximize']
    kpts=self.get_kpts(job_record)
    infns = [] #MAXIMIZE inputs
    outfiles=[]
    for k in kpts:
      for n in options['nconfig']:
        for jast in options['jastrow']:
          for opt in options['optimizer']:
            infns.append(self.gen_basename(k,n,jast,opt))
            #basename = self.make_basename(k, n, w)
            #outfiles.extend([basename+".max.table",
            #                basename+".max.log",
            #                basename+".max.config",
            #                basename+".max.o"])

    #Check on the submitter. If still running report that.
    status=self._submitter.status(job_record,self._name_)
    print("status",status)
    if 'running' in status:
      return 'running'

    #If not running, try to transfer files.
    self._submitter.transfer_output(job_record, infns)

    #Now check on the runs
    ret=self.collect_runs(job_record)
    if len(ret)==0:
      return "not_started"
    if len(ret) != len(infns):
      print("There are no jobs running and not enough .json files. Not sure what's going on.")
      quit()
    return 'ok'

#-----------------------------------------------
  def resume(self,job_record):
    return self.run(job_record,restart=True)

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
              if not os.path.exists(basename+".trace"):
                print("You need a trace file to run postprocess.")
                return "failed"
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
  def swap_endian(self,tracefn):
    return newtracefn

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
      # Sometimes output files are large, and this takes forever. 
      # Instead just search the end.
      if "Wall" in str(sub.check_output(["tail",outfilename])):
        return 'ok'
      else:
        return 'running'
    else:
      return 'not_started'

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
    tracefn = basename+".trace"
    if ppr_options['swap_endian']:
      newtracefn = tracefn.replace(".trace",".swap.trace")
      if not os.path.exists(newtracefn):
        print("Swapping."+sub.check_output([swap_endian_exe,tracefn,newtracefn]))
      tracefn = newtracefn
    outlines = [
        "method { postprocess ",
        "readconfig %s"%tracefn,
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
    jast_inp=extract_jastrow(open("qw_0.%s.%s.wfout"%(jast,opt_trans[opt])))
    outlines += [
        "}",
        "include qw_%i.sys"%k,
        "trialfunc { slater-jastrow ",
           "wf1 { include qw_%i.slater } "%(k),
           "wf2 { ",
          jast_inp,
          " } ",
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
              if os.path.isfile("%s.post.json"%basename):
                entry={}
                entry['knum']=k
                entry['timestep']=t
                entry['localization']=loc
                entry['jastrow']=jast
                entry['optimizer']=opt
                entry['results']=json.load(open("%s.post.json"%basename))
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
    elif len(ret) == len(infns):
      return "ok"
    else:
      print("There are no jobs running and not enough .log files. Not sure what's going on.")
      quit()

  def resume(self,job_record):
    return self.run(job_record,restart=True)

  #-----------------------------------------------
  def output(self,job_record):
    job_record['qmc']['postprocess']['results']=self.collect_runs(job_record)
    return job_record

####################################################

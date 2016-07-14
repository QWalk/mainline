from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
####################################################

class RunCrystal:
  _name_="RunCrystal"
  def __init__(self, submitter):
    self._submitter = submitter
#-------------------------------------------------      
  def run(self, job_record):
    self._submitter.execute(job_record, ['autogen.d12'],
          'autogen.d12', 'autogen.d12.o',self._name_)
    return 'running'
#-------------------------------------------------      
  def output(self,job_record):
    """ Collect results from output."""
    if os.path.isfile('autogen.d12.o'):
      f = open('autogen.d12.o', 'r')
      lines = f.readlines()
      for li,line in enumerate(lines):
        if 'SCF ENDED' in line:
          job_record['dft']['total_energy']=float(line.split()[8])    
        elif 'TOTAL ATOMIC SPINS' in line:
          moms = []
          shift = 1
          while "TTT" not in lines[li+shift]:
            moms += map(float,lines[li+shift].split())
            shift += 1
          job_record['dft']['mag_moments']=moms
      
    return job_record
#-------------------------------------------------      
  # This can be made more efficient if it's a problem: searches whole file for
  # each query.
  def check_outputfile(self,outfilename,acceptable_scf=10.0):
    """ Check output file. 

    Current return values:
    no_record, not_started, ok, too_many_cycles, finished (fall-back),
    scf_fail, not_enough_decrease, divergence, not_finished
    """
    if os.path.isfile("autogen.d12.o"):
      outf = open("autogen.d12.o",'r')
    else:
      return "not_started"
    outlines = outf.read().split('\n')
    reslines = [line for line in outlines if "ENDED" in line]
    if len(reslines) > 0:
      if "CONVERGENCE" in reslines[0]:
        return "ok"
      elif "TOO MANY CYCLES" in reslines[0]:
        print("RunCrystal output: Too many cycles.")
        return "too_many_cycles"
      else: # What else can happen?
        print("RunCrystal output: Finished, but unknown state.")
        return "finished"
    detots = [float(line.split()[5]) for line in outlines if "DETOT" in line]
    if len(detots) == 0:
      print("RunCrystal output: Last run completed no cycles.")
      return "scf_fail"
    detots_net = sum(detots[1:])
    if detots_net > acceptable_scf:
      print("RunCrystal output: Last run performed poorly.")
      return "not_enough_decrease"
    etots = [float(line.split()[3]) for line in outlines if "DETOT" in line]
    if etots[-1] > 0:
      print("RunCrystal output: Energy divergence.")
      return "divergence"
    print("RunCrystal output: Not finished.")
    return "not_finished"
#-------------------------------------------------      
  # Diagnose routines basically decide 'not_finished' or 'failed'
  def stubborn_diagnose(self,status):
    if status in ['too_many_cycles','not_finished']:
      return 'not_finished'
    else:
      return 'failed'
#-------------------------------------------------      
  def optimistic_diagnose(self,status):
    if status == 'not_finished':
      return 'not_finished'
    else:
      return 'failed'
#-------------------------------------------------      
  def conservative_diagnose(self,status):
    return 'failed'
#-------------------------------------------------      
  def check_status(self,job_record):
    """ Decide status of job (in queue or otherwise). """
    outfilename="autogen.d12.o"
    diagnose_options = {
        'stubborn':self.stubborn_diagnose,
        'optimistic':self.optimistic_diagnose,
        'conservative':self.conservative_diagnose
      }
    if job_record['dft']['resume_mode'] not in diagnose_options.keys():
      print("RunCrystal: Diagnose option not recognized: falling back on 'conservative'")
      diagnose = diagnose_options['conservative']
    else:
      diagnose = diagnose_options[job_record['dft']['resume_mode']]

    status=self._submitter.status(job_record,self._name_)
    if 'running' in status:
      return 'running'

    self._submitter.transfer_output(job_record, ['autogen.d12.o', 'fort.9'])
    status=self.check_outputfile(outfilename)
    print("status",status)
    if status == 'not_started':
      return 'not_started'
    elif status == 'ok':
      if 'initial_spin' in job_record['dft'].keys():
        if not self._consistent_spins(outfilename,job_record['dft']['initial_spin']):
          print("Error: Spin changed from initial configuration!")
          return self.conservative_diagnose("changed_spins")
      return 'ok'
    status=diagnose(status)
    if status in ['ok','not_finished','failed']:
      return status
    # This case shouldn't happen:
    if not os.path.isfile(outfilename):
      print("Warning: author of this code didn't expect this to occur!")
      return 'not_started'

    return 'failed'
#-------------------------------------------------      
  def _consistent_spins(self,outfn,init_spins,small_spin=0.1):
    f = open(outfn, 'r')
    lines = f.readlines()
    moms=[]
    for li,line in enumerate(lines):
      if 'TOTAL ATOMIC SPINS' in line:
        moms = []
        shift = 1
        while "TTT" not in lines[li+shift]:
          moms += map(float,lines[li+shift].split())
          shift += 1
    moms = np.array(moms)
    zs = abs(moms) < small_spin
    up = moms > 0.
    dn = moms < 0.
    moms.dtype = int
    moms[up] = 1
    moms[dn] = -1
    moms[zs] = 0
    moms=np.resize(moms,(len(init_spins)))
    if len(init_spins)==0:
      return True
    else:
      return np.array_equal(moms,np.array(init_spins))
#-------------------------------------------------      
  def add_guessp(self,inpfn):
    inplines = open(inpfn,'r').read().split('\n')
    if "GUESSP" not in inplines:
      cursor = -1
      while "END" not in inplines[cursor]:
        cursor -= 1
      inplines.insert(cursor,"GUESSP")
    with open(inpfn,'w') as outf:
      outf.write('\n'.join(inplines))
    return None
#-------------------------------------------------      
  def resume(self,job_record,maxresume=5):
    """ Continue a crystal run using GUESSP."""
    jobname = job_record['control']['id']
    trynum = 0
    while os.path.isfile("%d.autogen.d12.o"%trynum):
      trynum += 1
      if trynum > maxresume:
        print("Not resuming because resume limit reached ({}>{}).".format(
          trynum,maxresume))
        return 'failed'
    for filename in ["autogen.d12","autogen.d12.o","fort.79"]:
      shutil.copy(filename,"%d.%s"%(trynum,filename))
    shutil.copy("fort.79","fort.20")
    job_record['dft']['restart_from'] = os.getcwd()+"%d.fort.79"%trynum
    self.add_guessp("autogen.d12")
    return self.run(job_record)

####################################################

class RunProperties:
  _name_="RunProperties"
  def __init__(self,submitter=None):
    # 'None' implies to run in command line.
    self._submitter=submitter
  def run(self,job_record):
    f=open("prop.in",'w')
    if 'cif' in job_record.keys():
      kmax = max(job_record['dft']['kmesh'])
      out = '\n'.join([
        "NEWK",
        "%d %d"%(kmax,2*kmax),
        "1 1",
        "67 999",
        "END"
      ])
      f.write(out)
    else:
      out = '\n'.join([
        "NEWK",
        "1 1",
        "67 999",
        "END"
      ])
      f.write(out)
    f.close()

    if self._submitter==None:
      os.system("properties < prop.in > prop.in.o")
      return 'ok'
    else:
      self._submitter.execute(job_record,["prop.in"], 'prop.in', 'prop.in.o',self._name_)
      return 'running'

  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      # Sometimes prop.in.o files are GB large, and this takes forever. 
      # Instead just search the end.
      if "ENDPROP" in str(sub.check_output(["tail",outfilename])):
        return 'ok'
      else:
        return 'running'
      #for line in f:
      #  if "ENDPROP" in line:
      #    return 'ok'
      #return 'running'
    else:
      return 'not_started'

  def check_status(self,job_record):
    outfilename="prop.in.o"
    status=self.check_outputfile(outfilename)
    if status=='ok' or status=='failed':
      return status

    if self._submitter!=None:
      status=self._submitter.status(job_record,self._name_)
      if 'running' in status:
        return 'running'
      self._submitter.transfer_output(job_record, [outfilename, 'fort.9'])
      status=self.check_outputfile(outfilename)
      if status=='ok':
        self._submitter.cancel(job_record['control']['queue_id'])
        return status
      elif status=='not_finished' or status=='failed':
        return status
    
    if not os.path.isfile(outfilename):
      return 'not_started'

    return 'failed'
      
  def retry(self,job_record):
    return self.run(job_record)

  def output(self,job_record):
    return job_record

####################################################

class NewRunProperties:
  _name_="RunProperties"
  def __init__(self,submitter=None):
    # 'None' implies to run in command line.
    self._submitter=submitter
  def run(self,job_record):
    f=open("prop.in",'w')
    if 'cif' in job_record.keys():
      kmax = max(job_record['dft']['kmesh'])
      out = '\n'.join([
        "NEWK",
        "%d %d"%(kmax,2*kmax),
        "1 0",
        "CRYAPI_OUT",
        "END"
      ])
      f.write(out)
    else:
      out = '\n'.join([
        "NEWK",
        "1 0",
        "CRYAPI_OUT",
        "END"
      ])
      f.write(out)
    f.close()

    if self._submitter==None:
      os.system("properties < prop.in > prop.in.o")
      return 'ok'
    else:
      self._submitter.execute(job_record,["prop.in"], 'prop.in', 'prop.in.o',self._name_)
      return 'running'

  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if "ENDPROP" in line:
          return 'ok'
      return 'running'
    else:
      return 'not_started'

  def check_status(self,job_record):
    outfilename="prop.in.o"
    status=self.check_outputfile(outfilename)
    if status=='ok' or status=='failed':
      return status

    if self._submitter!=None:
      status=self._submitter.status(job_record,self._name_)
      if status=='running':
        return status
      self._submitter.transfer_output(job_record, [outfilename, 'fort.9'])
      status=self.check_outputfile(outfilename)
      if status=='not_finished' or status=='failed':
        return status
    
    if not os.path.isfile(outfilename):
      return 'not_started'

    return 'failed'
      
  def retry(self,job_record):
    return self.run(job_record)

  def output(self,job_record):
    return job_record

####################################################
  

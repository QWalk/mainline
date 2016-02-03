from __future__ import print_function
import os
import shutil
####################################################

class RunCrystal:
  _name_="RunCrystal"
  def __init__(self, submitter):
    self._submitter = submitter

  def run(self, job_record):
    job_record['control'][self._name_+'_jobid'] = \
        [self._submitter.execute(job_record, ['autogen.d12'],
          'autogen.d12', 'autogen.d12.o')]
    return 'running'

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

  # This can be made more efficient if it's a problem: searches whole file for
  # each query.
  def check_outputfile(self,outfilename,acceptable_scf=0.0):
    """ Check output file. 
    Current return values:
    no_record, not_started, ok, too_many_cycles, finished (fall-back),
    scf_fail, not_enough_decrease, divergence, not_finished
    """
    try:
      outf = open("autogen.d12.o",'r')
    except IOError:
      return "not_started"
    outlines = outf.read().split('\n')
    reslines = [line for line in outlines if "ENDED" in line]
    if len(reslines) > 0:
      if "CONVERGENCE" in reslines[0]:
        return "ok"
      elif "TOO MANY CYCLES" in reslines[0]:
        return "too_many_cycles"
      else: # What else can happen?
        return "finished"
    detots = [float(line.split()[5]) for line in outlines if "DETOT" in line]
    if len(detots) == 0:
      return "scf_fail"
    detots_net = sum(detots[1:])
    if detots_net > acceptable_scf:
      return "not_enough_decrease"
    etots = [float(line.split()[3]) for line in outlines if "DETOT" in line]
    if etots[-1] > 0:
      return "divergence"
    return "not_finished"

  # Diagnose routines basically decide 'not_finished' or 'failed'
  def stubborn_diagnose(self,status):
    if status in ['too_many_cycles','not_finished']:
      return 'not_finished'
    else:
      return 'failed'

  def optimistic_diagnose(self,status):
    if status == 'not_finished':
      return 'not_finished'
    else:
      return 'failed'

  def conservative_diagnose(self,status):
    return 'failed'

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

    status=self.check_outputfile(outfilename)
    if status in ['ok','failed']:
      return status

    self._submitter.transfer_output(job_record, ['autogen.d12.o', 'fort.9'])
    status=self._submitter.status(job_record)
    print("Submitter:",status)
    if status=='running':
      return status
    status=self.check_outputfile(outfilename)
    if status == 'not_started':
      return status
    status=diagnose(status)
    if status in ['ok','not_finished','failed']:
      return status
    # This case shouldn't happen:
    if not os.path.isfile(outfilename):
      print("Turns out this case happens!")
      return 'not_started'

    return 'failed'


  def resume(self,job_record,maxresume=5):
    """ Continue a crystal run using GUESSP."""
    jobname = job_record['control']['id']
    trynum = 0
    while os.path.isfile(str(trynum)+".autogen.d12.o"):
      trynum += 1
      if trynum > maxresume:
        print("Not resuming because resume limit reached ({}>{}).".format(
          trynum,maxresume))
        return 'failed'
    for filename in ["autogen.d12","autogen.d12.o"]:
      shutil.move(filename,str(trynum)+"."+filename)
    for filename in ["fort.79"]:
      shutil.copy(filename,str(trynum)+"."+filename)
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
      job_record['control'][self._name_+'_jobid'] = \
          [self._submitter.execute(job_record,["prop.in"], 'prop.in', 'prop.in.o')]
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
      status=self._submitter.status(job_record)
      if status=='running':
        return status
      self._submitter.transfer_output(job_record, [outfilename, 'fort.9'])
      status=self.check_outputfile(outfilename)
      if status=='ok':
        self._submitter.cancel(job_record['control'][self._name_+'_jobid'])
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
  

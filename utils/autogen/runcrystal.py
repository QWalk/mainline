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

  def check_outputfile(self,outfilename):
    """ Check if outputfile reports sucess. """
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if "SCF ENDED" in line:
          if "TOO MANY CYCLES" in line:
            print("Crystal failed: too many cycles.")
            return 'not_finished'
          energy = float(line.split()[8])
          if energy > 0.0:
            print("Crystal failed: energy divergence.") 
            return 'failed'

          return 'ok'
    else:
      return 'not_started'

  def check_status(self,job_record):
    """ Decide status of job (in queue or otherwise). """
    outfilename="autogen.d12.o"

    status=self.check_outputfile(outfilename)
    if status=='failed':
      return status
    elif status=='ok':
      return status

    self._submitter.transfer_output(job_record, ['autogen.d12.o', 'fort.9'])
    status=self._submitter.status(job_record)
    if status=='running':
      return status
    status=self.check_outputfile(outfilename)
    if status=='ok':
      return status
    elif status=='not_finished' or status=='failed':
      return status
  
    if not os.path.isfile(outfilename):
      return 'not_started'

    return 'failed'

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
  

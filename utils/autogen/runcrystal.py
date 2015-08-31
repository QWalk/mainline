from __future__ import print_function
import os
import job_submission
import shutil
####################################################

class RunCrystal:
  _name_="RunCrystal"
  _submitter=job_submission.TorqueCrystalSubmitter()

  def __init__(self, submitter=job_submission.TorqueCrystalSubmitter()):
    self._submitter = submitter

  def run(self, job_record):
    job_record['control'][self._name_+'_jobid'] = [self._submitter.execute(job_record, ['autogen.d12'], 'autogen.d12', 'autogen.d12.o')]
    return 'running'

  def output(self,job_record):
    if os.path.isfile('autogen.d12.o'):
      f = open('autogen.d12.o', 'r')
      lines = f.readlines()
      for l in lines:
        if 'SCF ENDED' in l:
          job_record['dft']['total_energy']=float(l.split()[8])    
      
    return job_record


  def check_outputfile(self,outfilename):
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
    outfilename="autogen.d12.o"

    status=self.check_outputfile(outfilename)
    if status=='failed':
      return status
    elif status=='ok':
      self._submitter.cancel(job_record['control'][self._name_+'_jobid'])
      return status

    self._submitter.output(job_record, ['autogen.d12.o', 'fort.9'])
    status=self._submitter.status(job_record)
    if status=='running':
      return status
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
    """Copy fort.9 to fort.20 and add GUESSP if it isn't already there."""
    #Removing this behavior for now as it doesn't seem to help too much.
    #shutil.copy('fort.9','fort.20')
    #with open('autogen.d12','r') as d12f:
    #  lines = d12f.read().split('\n')
    #if not any(["GUESSP" in line for line in lines]):
      # Currently autogen doesn't end the file with \n (e.g. "END\n"), 
      # this will fail if in the future it does.
    #  lines[-1] = "GUESSP\nEND"
    #  with open('autogen.d12','w') as d12f:
    #    d12f.write('\n'.join(lines))
    #job_record['dft']['nretries'] += 1
    return "did_not_retry" #self.run(job_record)

####################################################

class RunProperties:
  _name_="RunProperties"
  #_submitter=job_submission.TorquePropertiesSubmitter()
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
      job_record['control'][self._name_+'_jobid'] = self._submitter.execute(job_record,["prop.in"], 'prop.in', 'prop.in.o')
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
      self._submitter.output(job_record, [outfilename, 'fort.9'])
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
  

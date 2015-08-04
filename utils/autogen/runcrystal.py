from __future__ import print_function
import os
import job_submission
import shutil
####################################################

class RunCrystal:
  _name_="RunCrystal"
  _submitter=job_submission.TorqueCrystalSubmitter()
  def __init__(self,submitter=job_submission.TorqueCrystalSubmitter()):
    self._submitter=submitter

  def run(self,job_record):
    self._submitter.execute(job_record,['autogen.d12'])
    return 'running'

  def check_outputfile(self,outfilename):
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if "SCF ENDED" in line:
          if "TOO MANY CYCLES" in line:
            print("Crystal failed: too many cycles.")
            return 'not_finished'
          return 'ok'


  def check_status(self,job_record):
    d=str(job_record['control']['id'])
    outfilename="autogen.d12.o"
    status=self.check_outputfile(outfilename)
    if status=='ok' or status=='failed':
      return status

    status=self._submitter.status(job_record,[outfilename,'fort.9'])
    if status=='running':
      return status
    status=self.check_outputfile(outfilename)
    if status=='ok' or status=='not_finished' or status=='failed':
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

  def output(self,job_record):
    outfilename="autogen.d12.o"
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if "SCF ENDED" in line:
          job_record['dft']['total_energy']=float(line.split()[8])
    return job_record

####################################################

class RunProperties:
  _name_="RunProperties"
  def run(self,job_record):
    f=open("prop.in",'w')
    if 'cif' in job_record.keys():
      f.write("""NEWK 
4 4 
1 1
67 999
END 
""")
    else:
      f.write("""NEWK 
1 1
67 999
END 
""")

    f.close()
    os.system("properties < prop.in > prop.in.o")
    return 'ok'
  def check_status(self,job_record):
    outfilename="prop.in.o"
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if "ENDPROP" in line:
          return 'ok'
      return 'running'
    else:
      return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    return job_record

####################################################
  

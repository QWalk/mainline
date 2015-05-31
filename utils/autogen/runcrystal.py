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
          return 'ok'


  def check_status(self,job_record):
    d=str(job_record['control']['id'])
    outfilename="autogen.d12.o"
    if self.check_outputfile(outfilename)=='ok':
      return 'ok'

    status=self._submitter.status(job_record,[outfilename,'fort.9'])
    if status=='running':
      return status
    if self.check_outputfile(outfilename)=='ok':
      return 'ok'

    return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)

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
  

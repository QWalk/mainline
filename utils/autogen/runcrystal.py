import os
####################################################

class RunCrystal:
  def run(self,job_record):
    d=str(job_record['control']['id'])+"/"
    os.system("cd %s; crystal < autogen.d12 > autogen.d12.o"%d)
    return 'ok'
  def check_status(self,job_record):
    outfilename=str(job_record['control']['id'])+"/autogen.d12.o"
    if os.path.isfile(outfilename):
      f=open(outfilename,'r')
      for line in f:
        if "SCF ENDED" in line:
          return 'ok'
      return 'running'
    else:
      return 'not_started'
      
  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    return job_record
####################################################

class RunProperties:
  def run(self,job_record):
    d=str(job_record['control']['id'])+"/"
    f=open(d+"prop.in",'w')
    f.write("""NEWK 
4 4 
1 1
67 999
END 
""")
    f.close()
    os.system("cd %s; properties < prop.in > prop.in.o"%d)
    return 'ok'
  def check_status(self,job_record):
    outfilename=str(job_record['control']['id'])+"/prop.in.o"
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
  

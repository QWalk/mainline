import subprocess
import importlib
import os
import shutil
import sys
import time

#####################################################################################

class RemoteSubmitter:
  def __init__(self, px_ssh, px_ftp, module, remotePath):
    sys.path.append(os.getcwd() + '/Specific/')
    self.px_ssh = px_ssh
    self.px_ftp = px_ftp
    self.remotePath = remotePath
    self.module = importlib.import_module(module)
#---------------------------------------------------
  def execute(self, job_record, infiles, runfile, outfile):
    print("remotepath",self.remotePath)
    job_id = str(job_record['control']['id'])
    self.px_ssh.sendline('mkdir -p ' + self.remotePath + job_id)
    self.px_ssh.prompt()
    self.px_ssh.sendline('cd ' + self.remotePath + job_id)
    self.px_ssh.prompt()

    for dep in infiles:  
      print('Transferring file _to_ remote cluster:   ' + dep)     
      self.px_ftp.expect('sftp>', timeout=None)
      self.px_ftp.sendline('put ' + os.getcwd() + '/' + dep + ' ' + self.remotePath + job_id)
      

    print('Submitting to queue...')
    job_record['control']['queue_id'] = [self.module.execute(self.px_ssh, runfile, outfile, str(job_record['control']['id']))]
    return job_record['control']['queue_id'][0]
#---------------------------------------------------
  def status(self, job_record):
    job_id = str(job_record['control']['id'])
    self.px_ssh.sendline('cd ' + self.remotePath + job_id)
    self.px_ssh.prompt()
    statuses = []
    for q_id in job_record['control']['queue_id']:
      statuses.append(self.module.status(self.px_ssh, q_id))
    if 'running' in statuses:
      return 'running'
    else:
      return 'not_running'
#---------------------------------------------------
  def transfer_output(self, job_record, outfiles):
    job_id = str(job_record['control']['id'])
    self.px_ssh.sendline('cd ' + self.remotePath + job_id)
    self.px_ssh.prompt()
    self.px_ssh.sendline('ls')
    self.px_ssh.prompt()
    remote_files = self.px_ssh.before.decode('utf-8')
    
    for cop in outfiles:
      if cop in remote_files:
        print('Transferring file _from_ remote cluster:   ' + cop)
        self.px_ftp.expect('sftp>', timeout=None)
        self.px_ftp.sendline('get ' + self.remotePath + job_id + '/' + cop + ' ' + os.getcwd())
       
    self.px_ftp.expect('sftp>', timeout=None)
    self.px_ftp.sendline('pwd')
#---------------------------------------------------      
  def cancel(self, queue_id):
    output = []
    for q_id in queue_id:
      output.append(self.module.cancel(self.px_ssh, q_id))

#####################################################################################

class LocalSubmitter:
  """Abstract submission class. Child classes must define:
  __init__: 
     can have any parameters, but should set up any variables like queue time and number of
     processor cores
  _job_status(self,
              queue_id : a string that identifies the queue id of the job
              )
  _submit_job(self,
              inpfns : list of input filenames (list of strings)
              outfn  : where stdout should be sent (string)
              jobname : job name for the queue (string)
              loc :  directory where job should be run (string)
            )
            returns a list of queue ids (list of strings)
  """
#---------------------------------------------------  
  def execute(self,job_record,dependencies,inpfns,outfn,name):
    """Generate qsub file for this job, run it, and return qid from qsub
    transaction. 
    Adds a an element in job_record['control']['queue_id'] with the queue id
    qid: is a list
    """
    qid = self._submit_job(
        inpfns,
        outfn = outfn,
        jobname = job_record['control']['id'],
        loc = os.getcwd()
      )
    if not 'queue_id' in job_record['control'].keys():
      job_record['control']['queue_id']=[]
    for q in qid:
      job_record['control']['queue_id'].append([name,q])
    return qid
#---------------------------------------------------
  def status(self,job_record,name):
    """ Returns a list of job status elements """
    if not 'queue_id' in job_record['control']:
      return ['unknown']
    status=[]
    for q in job_record['control']['queue_id']:
      if q[0]==name:
        status.append(self._job_status(q[1]))
    return status
#---------------------------------------------------
  def transfer_output(self,job_record,outfiles):
    pass # Files should be already available locally.
#---------------------------------------------------
  def cancel(self, queue_id):
    output = []
    for q_id in queue_id:
      output.append(self._job_cancel(q_id))
#####################################################################################
      

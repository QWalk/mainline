import subprocess
import importlib
import os
import shutil
import sys
import time

class RemoteSubmitter:
  def __init__(self, px_ssh, px_ftp, module, remotePath):
    sys.path.append(os.getcwd() + '/Specific/')
    self.px_ssh = px_ssh
    self.px_ftp = px_ftp
    self.remotePath = remotePath
    self.module = importlib.import_module(module)

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
      
  def cancel(self, queue_id):
    output = []
    for q_id in queue_id:
      output.append(self.module.cancel(self.px_ssh, q_id))

class LocalSubmitter:
  """Abstract submission class. Child classes must define initialization and
  internal functions like _submit_job"""
  
  def execute(self,job_record,dependencies,inpfns,outfn):
    """Generate qsub file for this job, run it, and return qid from qsub
    transaction."""
    qid = job_record['control']['queue_id'] = self._submit_job(
        inpfns,
        outfn = outfn,
        jobname = job_record['control']['id'],
        loc = os.getcwd()
      )
    job_record['control']['queue_id'] = qid
    #job_record['control']['queue_id'] = [qid, job_record['control']['id']]
    return qid

  def status(self,job_record):
    try:
      status = self._job_status(job_record['control']['queue_id'])
    except KeyError:
      status = "unknown"
    return status

  def transfer_output(self,job_record,outfiles):
    pass # Files should be already available locally.

  def cancel(self, queue_id):
    output = []
    for q_id in queue_id:
      output.append(self._job_cancel(q_id))

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
      # print('Transferring file _to_ remote cluster:   ' + dep)     
      self.px_ftp.sendline('put ' + os.getcwd() + '/' + dep + ' ' + self.remotePath + job_id)
      self.px_ftp.expect('sftp>')

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

  def output(self, job_record, outfiles):
    job_id = str(job_record['control']['id'])
    self.px_ssh.sendline('cd ' + self.remotePath + job_id)
    self.px_ssh.prompt()
    self.px_ssh.sendline('ls')
    self.px_ssh.prompt()
    remote_files = self.px_ssh.before.decode('utf-8')
    
    for cop in outfiles:
      if cop in remote_files:
        # print('Transferring file _from_ remote cluster:   ' + cop)
        self.px_ftp.sendline('get ' + self.remotePath + job_id + '/' + cop + ' ' + os.getcwd())
        self.px_ftp.expect('sftp>')
      
  def cancel(self, queue_id):
    output = []
    for q_id in queue_id:
      output.append(self.module.cancel(self.px_ssh, q_id))

class LocalSubmitter:
  def __init__(self, module):
    sys.path.append(os.getcwd() + '/Specific/')
    self.module = importlib.import_module(module)

  def execute(self,job_record,infiles,runfile, outfile):
    job_record['control']['queue_id'] = self.module.execute(runfile, outfile)

  def status(self,job_record):
    statuses = []
    for q_id in job_record['control']['queue_id']:
      statuses.append(self.module.status(self.px_ssh, q_id))
    if 'running' in statuses:
      return 'running'
    else:
      return 'not_running'

  def output(self,job_record,outfiles):
    pass

  def cancel(self, queue_id):
    output = []
    for q_id in queue_id:
      output.append(self.module.cancel(self.px_ssh, q_id))


# RunCrystal for Mira
import subprocess as sub
from submission_tools import LocalSubmitter
import os

class LocalMiraSubmitter(LocalSubmitter):
  """Abstract submission class. Defines interaction with the queuing system."""
  def __init__(self,time='72:00:00',nn=1,np='allprocs',mode=16,queue='batch'):
    """ Initialize a submitter object. 

    This part should contain all parts of calculation that user should care
    about, and all parts that are common to calculations of a certain type (all
    crystal SCF calculations, for instance)."""
    self.time  = time
    self.nn    = nn
    self.np    = np
    self.queue = queue
    self.mode  = mode

  def _job_status(self,queue_id):
    status = "unknown"
    queue_id=str(queue_id)
    print("Queue id",queue_id)
    qstat=None
    if 'campuscluster' in queue_id:
      return 'finished'
    try:
      qstat_out = sub.check_output(
          "qstat %s"%queue_id, shell=True
        ).decode().split('\n')
      if(len(qstat_out) > 2):
        qstat=qstat_out[-2].split()[4]
      print(qstat)
    except sub.CalledProcessError:
      #print("Called process error")
      return "unknown"
    if qstat == "running" or qstat == "queued":
      print("Found running")
      return "running"
    if qstat == "exiting" or qstat == "killed":
      print("Found finished")
      return "finished"
    return status

  def _job_cancel(self,queue_id):
    print("Cancel was called, but not implemented")

class LocalMiraCrystalSubmitter(LocalMiraSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    print("Crystal should not be run on Mira! Skipping this!")
    return None

class LocalMiraPropertiesSubmitter(LocalMiraSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    print("Properties should not be run on Mira! Skipping this!")
    return None

class LocalMiraQwalkSubmitter(LocalMiraSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    print("Non-DMC should not be run on Mira! Skipping this!")
    return None

class LocalMiraDMCSubmitter(LocalMiraSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfns,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    time=self.time
    nproc=512*self.nn
    qsub = []
    qsub.append('qsub')
#    qsub.append('-q prod')
    qsub.append('-A '+self.queue)
    qsub.append('-t {time}'.format(time=time))
    qsub.append('-n {nproc}'.format(nproc=nproc))
    qsub.append('--mode c%d'%self.mode)
    qsub.append('-o {outfn}'.format(outfn=outfn))
    qsub.append('--jobname {jobname}'.format(jobname=jobname))
    qsub.append('~/bin/qwalk')
    qsub += inpfns
    qin = ' '.join(qsub)
    print(qin); qid = None
    qid = sub.check_output(qin,shell=True).decode().strip()
    print("Submitted as %s"%qid)
    return [qid]

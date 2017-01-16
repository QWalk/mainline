# RunCrystal for Taub
import subprocess as sub
from submission_tools import LocalSubmitter
import os

# Where are your executables stored?

##########################################################
class LocalTaubSubmitter(LocalSubmitter):
  """Abstract submission class. Defines interaction with the queuing system."""
#-------------------------------------------------------  
  def __init__(self,time='72:00:00',nn=1,np='allprocs', queue='batch',
               bin="~/bin/",
               pre_commands="module load openmpi/1.4-gcc+ifort"):
    """ Initialize a submitter object. 

    This part should contain all parts of calculation that user should care
    about, and all parts that are common to calculations of a certain type (all
    crystal SCF calculations, for instance)."""
    self.time  = time
    self.nn    = nn
    self.np    = np
    self.queue = queue
    self.bin   = bin
    self.pre_commands=pre_commands
#-------------------------------------------------------
  def _job_status(self,queue_id):
    status = "unknown"
    try:
      qstat = sub.check_output(
          "qstat %s"%queue_id, stderr=sub.STDOUT, shell=True
        ).decode().split('\n')[-2].split()[4]
    except sub.CalledProcessError:
      return "unknown"
    if qstat == "R" or qstat == "Q":
      return "running"
    if qstat == "C" or qstat == "E":
      return "finished"
    return status
#-------------------------------------------------------
  def _job_cancel(self,queue_id):
    print("Cancel was called, but not implemented")
#-------------------------------------------------------
  def _qsub(self,exe,prep_commands=[],final_commands=[],
      name="",stdout="",loc=""):
    """ Helper function for executable submitters. 
    Should work in most cases to simplify code."""

    if stdout=="": stdout="stdout"
    if loc=="": loc=os.getcwd()
    if name=="": name=stdout
    header = []
    header.append("#!/bin/bash")
    if self.np=="allprocs": 
      header.append("#PBS -l nodes=%d,flags=allprocs"%self.nn)
    else:
      header.append("#PBS -l nodes=%d:ppn=%d"%(self.nn,self.np))
    header += [
      "#PBS -q %s"%self.queue,
      "#PBS -l walltime=%s"%self.time,
      "#PBS -j oe",
      "#PBS -m n",
      "#PBS -N %s"%name,
      "#PBS -o {0}".format(loc+"/qsub.out"),
      self.pre_commands
    ]
    if self.np=="allprocs":
      exeline = "mpirun %s &> %s"%(exe, stdout)
    elif self.nn*self.np > 1:
      exeline = "mpirun -n %d %s &> %s"%(self.nn*self.np, exe, stdout)
    else:
      exeline = "%s &> %s"%(exe, stdout)
    commands = header + ["cd %s"%loc] + prep_commands + [exeline] + final_commands
    outstr = '\n'.join(commands)
    with open(loc+"/qsub.in",'w') as qsin:
      qsin.write(outstr)
    result = sub.check_output("qsub %s"%(loc+"/qsub.in"),shell=True)
    qid = result.decode().split()[0]
    print("Submitted as %s"%qid)
    return [qid]

###############################################################
class LocalTaubNullSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Do Nothing. Useful for getting autogen to gather output but submit
    no new jobs.
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    return [None]

###############################################################
class LocalTaubCrystalSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    exe = self.bin+"Pcrystal"
    prep_commands=["cp %s INPUT"%inpfn]
    final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    qid = self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)
    return qid

###############################################################
class LocalTaubPropertiesSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    exe = self.bin+"properties < %s"%inpfn
    prep_commands = []
    final_commands = []

    if self.nn != 1 or self.np != 1:
      print("Refusing to run properties in parallel!")
      self.nn = 1
      self.np = 1

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    qid = self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)
    return qid

#####################################################################
class LocalTaubQwalkSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfns,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    prep_commands=[]
    final_commands=[]
    qid=[]
    
    for f in inpfns:
      exe = self.bin+"qwalk %s"%f
      fulloutfn = f+outfn # Ugly but more functional.

      if jobname == "":
        jobname = outfn
      if loc == "":
        loc = os.getcwd()

      qid+=self._qsub(exe,prep_commands,final_commands,jobname,fulloutfn,loc)
    return qid

#####################################################################
class LocalTaubBundleQwalkSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfns,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    exe = " ".join([self.bin+"qwalk"]+inpfns)
    prep_commands=[]
    final_commands=[]

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    qid = self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)
    return qid
  

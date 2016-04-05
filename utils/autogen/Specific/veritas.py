<<<<<<< HEAD
# RunCrystal for Veritas
=======
# RunCrystal for Taub
>>>>>>> master
import subprocess as sub
from submission_tools import LocalSubmitter
import os

<<<<<<< HEAD
# Where are your executables stored?
BIN = "/home/brian/bin/"

if BIN[-1] != '/': BIN += '/'
##########################################################
class LocalVeritasSubmitter(LocalSubmitter):
  """Abstract submission class. Defines interaction with the queuing system."""
#-------------------------------------------------------  
=======
#####################################################################
class LocalVeritasSubmitter(LocalSubmitter):
  """Abstract submission class. Defines interaction with the queuing system."""
  # --------------------------------------------------------------
>>>>>>> master
  def __init__(self,time='72:00:00',nn=1,np='allprocs', queue='batch'):
    """ Initialize a submitter object. 

    This part should contain all parts of calculation that user should care
    about, and all parts that are common to calculations of a certain type (all
    crystal SCF calculations, for instance)."""
    self.time  = time,
    self.nn    = nn
    self.np    = np
    self.queue = queue
<<<<<<< HEAD
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
      status = "finished"
    return status
#-------------------------------------------------------
  def _job_cancel(self,queue_id):
    print("Cancel was called, but not implemented")
#-------------------------------------------------------
  def _qsub(self,exe,prep_commands=[],final_commands=[],
      name="",stdout="",loc=""):
    """ Helper function for executable submitters. 
=======

  # --------------------------------------------------------------
  # Finds if any jobs in queue_ids are still in the queue.
  def _job_status(self,queue_ids):
    if type(queue_ids) != type([]):
      # Try to keep a standard data type for the queue_ids.
      # It's a list of strings, all ids must be done before continuing.
      print("Warning: queue_ids had to be cast to list!")
      print("Debug:",queue_ids)
      qids = [queue_ids]
    else: qids = queue_ids
    status = "finished"
    for qid in qids:
      try:
        qstat = sub.check_output(
            "qstat %s"%qid, stderr=sub.STDOUT, shell=True
          ).decode().split('\n')[-2].split()[4]
      except sub.CalledProcessError:
        # Non-bundled jobs might finish and be removed before others.
        qstat = "C"
      if qstat == "R" or qstat == "Q":
        return "running"
    if qstat == "C" or qstat == "E":
      status = "finished"
    return status

  # --------------------------------------------------------------
  def _job_cancel(self,queue_id):
    print("Cancel was called, but not implemented")

  # --------------------------------------------------------------
  def _qsub(self,exe,prep_commands=[],final_commands=[],
      name="",stdout="",loc=""):
    """ Helper function for executible submitters. 
>>>>>>> master
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
<<<<<<< HEAD
      "#PBS -o {0}".format(loc+"/qsub.out"),
=======
      "#PBS -o {0}".format(loc+"/qsub.out")
>>>>>>> master
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
<<<<<<< HEAD
    return [qid]
###############################################################
=======
    return qid
>>>>>>> master

#####################################################################
class LocalVeritasCrystalSubmitter(LocalVeritasSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
<<<<<<< HEAD
    exe = BIN+"Pcrystal"
=======
    exe = "/home/brian/bin/Pcrystal"
>>>>>>> master
    prep_commands=["cp %s INPUT"%inpfn]
    final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    qid = self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)
    return qid
<<<<<<< HEAD
###############################################################
=======

#####################################################################
>>>>>>> master
class LocalVeritasPropertiesSubmitter(LocalVeritasSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
<<<<<<< HEAD
    exe = BIN+"properties < %s"%inpfn
=======
    exe = "/home/brian/bin/properties < %s"%inpfn
>>>>>>> master
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
<<<<<<< HEAD
###############################################################
=======
>>>>>>> master

#####################################################################
class LocalVeritasQwalkSubmitter(LocalVeritasSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfns,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
<<<<<<< HEAD
    prep_commands=[]
    final_commands=[]
    qid=[]
    
    for f in inpfns:
      exe = BIN+"qwalk %s"%f
=======
    qids = []
    for inpfn in inpfns:
      exe = "/home/brian/bin/qwalk %s"%inpfn
      prep_commands=[]
      final_commands=[]
>>>>>>> master

      if jobname == "":
        jobname = outfn
      if loc == "":
        loc = os.getcwd()

<<<<<<< HEAD
      qid+=self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)
    return qid
###############################################################
=======
      qids.append(self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc))
    return qids
>>>>>>> master

#####################################################################
class LocalVeritasBundleQwalkSubmitter(LocalVeritasSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfns,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
<<<<<<< HEAD
    exe = " ".join([BIN+"bin/qwalk"]+inpfns)
=======
    exe = " ".join(["/home/brian/bin/qwalk"]+inpfns)
>>>>>>> master
    prep_commands=[]
    final_commands=[]

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    qid = [self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)]
    return qid
<<<<<<< HEAD
###############################################################
  
=======
>>>>>>> master

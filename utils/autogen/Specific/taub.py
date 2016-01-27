# RunCrystal for Taub
import subprocess as sub
from submission_tools import LocalSubmitter
import os

class LocalTaubSubmitter(LocalSubmitter):
  """Abstract submission class. Defines interaction with the queuing system."""
  def __init__(self,time='72:00:00',nn=1,np='allprocs', queue='batch'):
    """ Initialize a submitter object. 

    This part should contain all parts of calculation that user should care
    about, and all parts that are common to calculations of a certain type (all
    crystal SCF calculations, for instance)."""
    self.time  = time,
    self.nn    = nn
    self.np    = np
    self.queue = queue

  def _job_status(self,queue_id):
    status = "unknown"
    try:
      qstat = sub.check_output(
          "qstat %s"%queue_id, stderr=sub.STDOUT, shell=True
        ).split('\n')[-2].split()[4]
    except sub.CalledProcessError:
      return "unknown"
    if qstat == "R" or qstat == "Q":
      return "running"
    if qstat == "C" or qstat == "E":
      return "finished"
    return status

  def _job_cancel(self,queue_id):
    print "Cancel was called, but not implemented"

class LocalTaubCrystalSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    exe = "/home/busemey2/bin/Pcrystal"
    prep_commands=[
        "module load openmpi/1.4-gcc+ifort",
        "cp %s INPUT"%inpfn]
    final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    header = []
    header.append('#!/bin/bash')
    if self.np=='allprocs': 
      header.append('#PBS -l nodes=%d,flags=allprocs'%self.nn)
    else:                  
      header.append('#PBS -l nodes=%d:ppn=%d'%(self.nn,self.np))
    header.append('#PBS -q %s'%self.queue)
    header.append('#PBS -l walltime=%s'%self.time)
    header.append('#PBS -j oe')
    header.append('#PBS -m n')
    header.append('#PBS -N %s'%jobname)
    header.append('#PBS -o {0}'.format(loc+'/qsub.out'))
    if self.np=='allprocs':
      exeline = 'mpirun %s &> %s'%(exe, outfn)
    elif self.nn*self.np > 1:
      exeline = 'mpirun -n %d %s &> %s'%(self.nn*self.np, exe, outfn)
    else:
      exeline = '%s &> %s'%(exe, outfn)
    commands = header + ['cd %s'%loc] + prep_commands + [exeline] + final_commands
    qsubstr = '\n'.join(commands)

    with open('qsub.in','w') as qsin:
      qsin.write(qsubstr)
    result = sub.check_output("qsub %s"%(loc+"/qsub.in"),shell=True)
    qid = result.split()[0]
    print "Submitted as %s"%qid
    return qid

class LocalTaubPropertiesSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    exe = "/home/busemey2/bin/properties < %s"%inpfn
    prep_commands = ["module load openmpi/1.4-gcc+ifort"]
    final_commands = []

    if self.nn != 1 or self.np != 1:
      print "Currently refusing to run properties in parallel!"
      self.nn = 1
      self.np = 1

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    header = []
    header.append('#!/bin/bash')
    if self.np=='allprocs': 
      header.append('#PBS -l nodes=%d,flags=allprocs'%self.nn)
    else:                  
      header.append('#PBS -l nodes=%d:ppn=%d'%(self.nn,self.np))
    header.append('#PBS -q %s'%self.queue)
    header.append('#PBS -l walltime=%s'%self.time)
    header.append('#PBS -j oe')
    header.append('#PBS -m n')
    header.append('#PBS -N %s'%jobname)
    header.append('#PBS -o {0}'.format(loc+'/qsub.out'))
    #if self.np=='allprocs':
    #  exeline = 'mpirun %s &> %s'%(exe, outfn)
    #elif self.nn*self.np > 1:
    #  exeline = 'mpirun -n %d %s &> %s'%(self.nn*self.np, exe, outfn)
    #else:
    exeline = '%s &> %s'%(exe, outfn)
    commands = header + ['cd %s'%loc] + prep_commands + [exeline] + final_commands
    qsubstr = '\n'.join(commands)

    with open('qsub.in','w') as qsin:
      qsin.write(qsubstr)
    result = sub.check_output("qsub %s"%(loc+"/qsub.in"),shell=True)
    qid = result.split()[0]
    print "Submitted as %s"%qid
    return qid

class LocalTaubQwalkSubmitter(LocalTaubSubmitter):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    exe = "/home/busemey2/bin/qwalk %s"%inpfn
    prep_commands = ["module load openmpi/1.4-gcc+ifort"]
    final_commands=[]

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    header = []
    header.append('#!/bin/bash')
    if self.np=='allprocs': 
      header.append('#PBS -l nodes=%d,flags=allprocs'%self.nn)
    else:                  
      header.append('#PBS -l nodes=%d:ppn=%d'%(self.nn,self.np))
    header.append('#PBS -q %s'%self.queue)
    header.append('#PBS -l walltime=%s'%self.time)
    header.append('#PBS -j oe')
    header.append('#PBS -m n')
    header.append('#PBS -N %s'%jobname)
    header.append('#PBS -o {0}'.format(loc+'/qsub.out'))
    if self.np=='allprocs':
      exeline = 'mpirun %s &> %s'%(exe, outfn)
    elif self.nn*self.np > 1:
      exeline = 'mpirun -n %d %s &> %s'%(self.nn*self.np, exe, outfn)
    else:
      exeline = '%s &> %s'%(exe, outfn)
    commands = header + ['cd %s'%loc] + prep_commands + [exeline] + final_commands
    qsubstr = '\n'.join(commands)

    with open('qsub.in','w') as qsin:
      qsin.write(qsubstr)
    result = sub.check_output("qsub %s"%(loc+"/qsub.in"),shell=True)
    qid = result.split()[0]
    print "Submitted as %s"%qid
    return qid

class LocalNullSubmitter(LocalTaubSubmitter):
  """NullSubmitter will not submit any jobs to the queue."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    return None

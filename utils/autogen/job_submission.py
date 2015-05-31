from __future__ import print_function
import subprocess
import os
import signal
import shutil

class LocalSubmitter:
  def execute(self, executable, argument,stdout,stdin):
    #os.system(executable+" "+argument+" > " + stdout + " &")
    stderr_f=open("STDERR",'w')
    stdout_f=open(stdout,'w')
    stdin_f=open(stdin,'r')
    pid=subprocess.Popen([executable,argument], 
            stdout=stdout_f,
            stderr=stderr_f,stdin=stdin_f).pid
    return pid
  def cancel(self, handle):
    print("handle",handle)
    pid=int(handle[0])
    os.kill(pid, signal.SIGQUIT) #or signal.SIGKILL 
    return "did_not_cancel"
  def status(self, handle):
    pid=int(handle[0])
    try:
        os.kill(pid, 0)
    except OSError, e:
        return "not_started"
    else:
        return "running"
  

class TorqueCrystalSubmitter:
  def execute(self, argument):
    shutil.copy(argument,"INPUT")
    f=open("QSUB",'w')
    f.write("""#PBS -l nodes=1:ppn=4
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -N autogen
#PBS -o QSUB.stdout
cd ${PBS_O_WORKDIR}
mpirun -np 4 /home/apps/bin/Pcrystal >& %s \n"""%(argument+".o"))
    f.close()
    output=subprocess.check_output(["qsub", "QSUB"])
    return output
  def cancel(self, handle):
    return "did_not_cancel"
  def status(self, handle):
    jobsign=handle[0]
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["qstat"])
    if jobnum in output:
      return "running"
    else:
      return "not_started"
  

  
class TorqueQWalkSubmitter:
  def execute(self, argument):
    f=open("QSUB",'w')
    f.write("""#PBS -l nodes=1:ppn=4
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -N autogen
#PBS -o QSUB.stdout
cd ${PBS_O_WORKDIR}
mpirun -np 4 ~/qwalk/bin/qwalk %s > %s \n"""%(argument,argument+".stdout"))
    f.close()
    output=subprocess.check_output(["qsub", "QSUB"])
    return output
  def cancel(self, handle):
    return "did_not_cancel"
  def status(self, handle):
    jobsign=handle[0]
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["qstat"])
    if jobnum in output:
      return "running"
    else:
      return "not_started"
  

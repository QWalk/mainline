from __future__ import print_function
import subprocess
import shutil
#------------Define how we run the calculation

class MyTorqueCrystalSubmitter:
  def execute(self, job_record,input_files):
    shutil.copy(input_files[0],"INPUT")
    f=open("QSUB",'w')
    f.write("""#PBS -l nodes=1:ppn=2
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -N autogen
#PBS -o QSUB.stdout
cd ${PBS_O_WORKDIR}
mpirun -np 2 /home/apps/bin/Pcrystal >& %s \n"""%(input_files[0]+".o"))
    f.close()
    print("submitting...")
    output=subprocess.check_output(["qsub", "QSUB"])
    job_record['control']['crystal_jobid']=output

  def cancel(self, handle):
    return "did_not_cancel"
  def status(self, job_record,output_files):
    if not 'crystal_jobid' in job_record['control']:
      return "not_started"
    jobsign=job_record['control']['crystal_jobid']
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["qstat"])
    if jobnum in output:
      return "running"
    else:
      return "not_started"


class MyTorqueQWalkSubmitter:
  def execute(self, job_record,input_files):
    argument=input_files[0]
    f=open("QSUB",'w')
    f.write("""#PBS -l nodes=1:ppn=6
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -N autogen
#PBS -o QSUB.stdout
cd ${PBS_O_WORKDIR}
mpirun -np 6 ~/qwalk/bin/qwalk %s > %s \n"""%(argument,argument+".stdout"))
    f.close()
    print("submitting...")
    output=subprocess.check_output(["qsub", "QSUB"])
    job_record['control']['qwalk_jobid']=output
    
    return output
  def cancel(self, handle):
    return "did_not_cancel"

  def status(self, job_record,output_files):
    if not 'qwalk_jobid' in job_record['control']:
      return "not_started"
    jobsign=job_record['control']['qwalk_jobid']
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["qstat"])
    if jobnum in output:
      return "running"
    else:
      return "not_started"



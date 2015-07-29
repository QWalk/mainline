from __future__ import print_function
import subprocess
import shutil
#------------Define how we run the calculation

class MyTorqueCrystalSubmitter:
  def __init__(self,nodes=1,ppn=16,time="48:00:00",queue="wagner"):
    self.nodes = nodes
    self.ppn   = ppn
    self.time  = time
    self.queue = queue
  def execute(self, job_record,input_files):
    shutil.copy(input_files[0],"INPUT")
    jobid=str(job_record['control']['id'])

    f=open("QSUB",'w')
    qf = '\n'.join([
      "#PBS -q %s"%self.queue,
      "#PBS -l nodes=%d:ppn=%d"%(self.nodes,self.ppn),
      "#PBS -l walltime=%s"%self.time,
      "#PBS -j oe",
      "#PBS -N autogen_%s"%jobid,
      "#PBS -o QSUB.stdout",
      "module load openmpi/1.4-gcc+ifort",
      "cd ${PBS_O_WORKDIR}",
      "mpirun -np %d ~/bin/Pcrystal >& %s "%(self.nodes*self.ppn,input_files[0]+".o"),
      "rm fort.*.pe*"
    ])
    f.write(qf)
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

class MyTorquePropertiesSubmitter:
  def __init__(self,nodes=1,ppn=1,time="4:00:00",queue="secondary"):
    self.nodes = 1
    self.ppn   = 1
    self.time  = time
    self.queue = queue
    if (nodes!=1) or (ppn!=1):
      print("Notice: RunProperties disobeys you, and will only be running on one processor.")
  def execute(self, job_record,input_files):
    jobid=str(job_record['control']['id'])

    f=open("QSUB",'w')
    qf = '\n'.join([
      "#PBS -q %s"%self.queue,
      "#PBS -l nodes=%d:ppn=%d"%(self.nodes,self.ppn),
      "#PBS -l walltime=%s"%self.time,
      "#PBS -j oe",
      "#PBS -N autogen_%s"%jobid,
      "#PBS -o QSUB.stdout",
      "module load openmpi/1.4-gcc+ifort",
      "cd ${PBS_O_WORKDIR}",
      "~/bin/properties < %s > %s"%(input_files[0],input_files[0]+".o")
    ])
    f.write(qf)
    f.close()
    print("submitting...")
    output=subprocess.check_output(["qsub", "QSUB"])
    job_record['control']['properties_jobid']=output

  def cancel(self, handle):
    return "did_not_cancel"
  def status(self, job_record,output_files):
    if not 'properties_jobid' in job_record['control']:
      return "not_started"
    jobsign=job_record['control']['properties_jobid']
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["qstat"])
    if jobnum in output:
      return "running"
    else:
      return "not_started"

class MyTorqueQWalkSubmitter:
  def __init__(self,nodes=1,ppn=16,time="48:00:00",queue="wagner"):
    self.nodes = nodes
    self.ppn   = ppn
    self.time  = time
    self.queue = queue
  def execute(self, job_record,input_files):
    argument=input_files[0]
    jobid=str(job_record['control']['id'])
    
    f=open("QSUB",'w')
    qf = '\n'.join([
      "#PBS -l nodes=%d:ppn=%d"%(self.nodes,self.ppn),
      "#PBS -q %s"%self.queue,
      "#PBS -l walltime=%s"%self.time,
      "#PBS -j oe",
      "#PBS -N autogen_%s"%jobid,
      "#PBS -o QSUB.stdout",
      "cd ${PBS_O_WORKDIR}",
      "module load openmpi/1.6.5-gcc-4.7.1 intel/14.0",
      "mpirun -np %d ~/qwalk/bin/qwalk %s > %s"%(self.nodes*self.ppn,argument,argument+".stdout")
    ])
    f.write(qf)
    f.close()
    print("submitting...")
    output=subprocess.check_output(["qsub", "QSUB"])
    if not 'qwalk_jobid' in job_record['control']:
      job_record['control']['qwalk_jobid']=[]

    job_record['control']['qwalk_jobid'].append(output)

  def cancel(self, handle):
    return "did_not_cancel"

  # Brian: This doesn't seem to be working, any reason this should be different
  # from MyTorqueCrystalSubmitter?
  #def status(self, job_record,output_files):
  #  if not 'qwalk_jobid' in job_record['control']:
  #    return "not_started"
  #  jobsign=job_record['control']['qwalk_jobid']
  #  output=subprocess.check_output(["qstat"])
  #  for line in output.split('\n'):
  #    for j in jobsign: 
  #      jobnum=j.split('.')[0]
  #      if jobnum in line:
  #        print("job",line)
  #        spl=line.split()
  #        if len(spl) > 9:
  #          s=line.split()[9 ]
  #          print("s:",s)
  #          if s!='C':
  #            print("Found running job")
  #            return "running"

  def status(self, job_record,output_files):
    if not 'qwalk_jobid' in job_record['control']:
      return "not_started"
    jobsign=job_record['control']['crystal_jobid']
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["qstat"])
    if jobnum in output:
      return "running"
    else:
      return "not_started"
    
    return "not_running"

import subprocess
import shutil
#------------Define how we run the calculation

class MyTorqueCrystalSubmitter:
  hostname="taub.campuscluster.illinois.edu"
  username="lkwagner"
  
  def execute(self,job_record,input_files):
    shutil.copy(input_files[0],"INPUT")
    directory="project-cse/autogen/"+str(job_record['control']['id'])
    
    f=open("QSUB",'w')
    f.write("""#PBS -q secondary
#PBS -l nodes=1:ppn=16
#PBS -l walltime=4:00:00
#PBS -j oe
#PBS -N autogen
#PBS -o QSUB.stdout
module load openmpi/1.4-gcc+ifort
cd ${PBS_O_WORKDIR}
mpirun -np 16 ~/bin/Pcrystal >& %s \n"""%(input_files[0]+".o"))
    f.close()
    input_files.append('QSUB')
    input_files.append('INPUT')
    #transfer input files to taub (making directory if needed)
    subprocess.check_output(['ssh',self.hostname,'mkdir -p %s'%directory])
    for fi in input_files:
      subprocess.check_output(['scp',fi,self.hostname+':'+directory])
    #submit job
    output=subprocess.check_output(["ssh", self.hostname,"cd "+directory+"; qsub QSUB"])
    job_record['control']['crystal_jobid']=output
    
  def cancel(self, handle):
    return "did_not_cancel"

  def status(self, job_record,output_files):
    directory="project-cse/autogen/"+str(job_record['control']['id'])
      
    if not 'crystal_jobid' in job_record['control']:
      return "not_started"
    jobsign=job_record['control']['crystal_jobid']
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["ssh",self.hostname,"qstat -u "+self.username])
    for line in output.split('\n'):
      if jobnum in line:
        print("job",line)
        spl=line.split()
        if len(spl) > 9:
          s=line.split()[9 ]
          if s!='C':
            return "running"
    
    #try to transfer output files back
    for fi in output_files:
      subprocess.check_output(['scp',self.hostname+':'+directory+"/"+fi,"."])
      
    return "not_running"
  


class MyTorqueQWalkSubmitter:
  hostname="taub.campuscluster.illinois.edu"
  username="lkwagner"
  
  def execute(self,job_record,input_files):
    shutil.copy(input_files[0],"INPUT")
    directory="project-cse/autogen/"+str(job_record['control']['id'])
    
    f=open("QSUB",'w')
    f.write("""#PBS -q secondary
#PBS -l nodes=1:ppn=16
#PBS -l walltime=4:00:00
#PBS -j oe
#PBS -N autogen
#PBS -o QSUB.stdout
module load openmpi/1.6.5-gcc-4.7.1 intel/14.0
cd ${PBS_O_WORKDIR}
mpirun -np 16 ~/qwalk/bin/qwalk %s >& %s \n"""%(input_files[0],input_files[0]+".stdout"))
    f.close()
    input_files.append('QSUB')
    #transfer input files to taub (making directory if needed)
    subprocess.check_output(['ssh',self.hostname,'mkdir -p %s'%directory])
    for fi in input_files:
      subprocess.check_output(['scp',fi,self.hostname+':'+directory])
    #submit job
    output=subprocess.check_output(["ssh", self.hostname,"cd "+directory+"; qsub QSUB"])
    job_record['control']['qwalk_jobid']=output
    
  def cancel(self, handle):
    return "did_not_cancel"

  def status(self, job_record,output_files):
    directory="project-cse/autogen/"+str(job_record['control']['id'])
      
    if not 'qwalk_jobid' in job_record['control']:
      return "not_started"
    jobsign=job_record['control']['qwalk_jobid']
    jobnum=jobsign.split('.')[0]
    output=subprocess.check_output(["ssh",self.hostname,"qstat -u "+self.username])
    for line in output.split('\n'):
      if jobnum in line:
        print("job",line)
        spl=line.split()
        if len(spl) > 9:
          s=line.split()[9 ]
          if s!='C':
            return "running"
    
    #try to transfer output files back
    try:
      for fi in output_files:
        subprocess.check_output(['scp',self.hostname+':'+directory+"/"+fi,"."])
    except:
      return "not_started"
      
    return "not_running"
  
    

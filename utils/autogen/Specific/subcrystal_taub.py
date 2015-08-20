# RunCrystal for Taub

def execute(px_ssh, infile, outfile, job_name):
  px_ssh.prompt()

  qscript = """#PBS -l nodes=1:ppn=12
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -q secondary
#PBS -o QSUB.stdout
#PBS -N %s

cd ${PBS_O_WORKDIR}
cp %s INPUT
module load openmpi/1.6.4-intel-13.1
mpirun -np 12 ~/bin/Pcrystal >& %s
rm fort.*.pe*
"""%(job_name, infile,outfile)
  
  px_ssh.sendline("echo '" + qscript + "' > batch_script")
  px_ssh.prompt()
  px_ssh.sendline('qsub batch_script')
  px_ssh.prompt()
  qsub_out=px_ssh.before.decode('utf-8').split('\n')[-2]
  print(qsub_out)
  return qsub_out.split('.')[0]

def status(px_ssh, queue_id):
  px_ssh.sendline("qstat | grep " + str(queue_id))
  px_ssh.prompt()

  try:
    stat = px_ssh.before.decode('utf-8').strip().split('\r\n')[1].split()[4]
  except:
    stat = ''

  if stat == 'R' or stat == 'Q':
    return 'running'
  else:
    return 'not_started'

def cancel(px_ssh, queue_id):
  px_ssh.sendline("qstat | grep " + str(queue_id))
  px_ssh.prompt()

  if queue_id in px_ssh.before.decode('utf-8'):
    px_ssh.sendline("qdel " + str(queue_id))
    px_ssh.prompt() 


# RunVarianceOptimization for Taub

def execute(px_ssh, infile, outfile, job_name):
	qscript = """#PBS -l nodes=1:ppn=16
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -q secondary
#PBS -o QSUB.stdout
#PBS -N %s

cd ${PBS_O_WORKDIR}
module load openmpi/1.6.5-gcc-4.7.1 intel/14.0
mpirun -np 16 ~/bin/qwalk %s >& %s"""%(job_name, infile, outfile)
	
	px_ssh.sendline("echo '" + qscript + "' > batch_script")
	px_ssh.prompt()
	px_ssh.sendline('qsub batch_script')
	px_ssh.prompt()
	return px_ssh.before.decode('utf-8').strip().split('\r\n')[1].split('.')[0]

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
		return 'not_running'

def cancel(px_ssh, queue_id):
	px_ssh.sendline("qstat")
	px_ssh.prompt()

	if queue_id in px_ssh.before.decode('utf-8'):
		px_ssh.sendline("qdel " + str(queue_id))
		px_ssh.prompt()	

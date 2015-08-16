# RunVarianceOptimization for Taub

def execute(px_ssh, infile, outfile):
	qscript = """#PBS -l nodes=1:ppn=12
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -q secondary
#PBS -o QSUB.stdout

cd ${PBS_O_WORKDIR}
module unload openmpi/1.6.4-intel-13.1
module load openmpi/1.4-intel
mpirun ~/bin/qwalk %s >& %s"""%(infile, outfile)
	
	px_ssh.sendline("echo '" + qscript + "' > batch_script")
	px_ssh.prompt()
	px_ssh.sendline('qsub batch_script')
	px_ssh.prompt()
	return px_ssh.before.strip().split('\r\n')[1].split('.')[0]

def status(px_ssh, queue_id):
	px_ssh.sendline("qstat | grep " + str(queue_id))
	px_ssh.prompt()

	try:
		stat = px_ssh.before.strip().split('\r\n')[1].split()[4]
	except:
		stat = ''

	if stat == 'R' or stat == 'Q':
		return 'running'
	else:
		return 'not_started'

def cancel(px_ssh, queue_id):
	px_ssh.sendline("qstat")
	px_ssh.prompt()

	if queue_id in px_ssh.before:
		px_ssh.sendline("qdel " + str(queue_id))
		px_ssh.prompt()		
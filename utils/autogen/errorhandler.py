


class ErrorHandler:
	def __init__(self,job_record):
		self.diagnose_options = {
		'stubborn':self.stubborn_diagnose,
		'optimistic':self.optimistic_diagnose,
		'conservative':self.conservative_diagnose
		}

		if job_record['dft']['resume_mode'] not in self.diagnose_options.keys():
			print("RunCrystal: Diagnose option not recognized: falling back on 'conservative'")
			self.diagnose = self.diagnose_options['conservative']
		else:
			self.diagnose = self.diagnose_options[job_record['dft']['resume_mode']]

	# Diagnose routines basically decide 'not_finished' or 'failed'
	def stubborn_diagnose(self,status):
		if status in ['too_many_cycles','not_finished']:
			return 'not_finished'
		else:
			return 'failed'
#-------------------------------------------------      
	def optimistic_diagnose(self,status):
		if status == 'not_finished':
			return 'not_finished'
		else:
			return 'failed'
#-------------------------------------------------      
	def conservative_diagnose(self,status):
		return 'failed'
#-------------------------------------------------   
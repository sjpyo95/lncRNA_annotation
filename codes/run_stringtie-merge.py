import sys, os
import module as mdl
import config as cf

stringtiePath = '/share/apps/programs/stringtie/stringtie-1.3.4a.Linux_x86_64/stringtie'

inputdir = sys.argv[1]
outputdir = sys.argv[2]
if not os.path.exists(outputdir): os.makedirs(outputdir)

gtffiles = filter(lambda x: '.gtf' in x, os.listdir(inputdir))
reps = dict()
for gtf in gtffiles:
	sample = gtf.split('_rep')[0]
	if not reps.has_key(sample):
		reps[sample] = []
	reps[sample].append(gtf)

for rep in reps.keys():
	i = 0
	gtfs = [inputdir + s for s in reps[rep]]
	sample = rep.split('_rep')[0]
	
	merge_CM = stringtiePath + ' --merge -p ' + cf.runThread_string + ' '.join(gtfs) + ' -o ' + outputdir + sample + '.M.gtf'
	queue = cf.queue('biglab', i)
	cf.qsub_time('PYO_job', jobN, 60)
	cf.qsub_execute('PYO_job', queue, merge_CM, outputdir+sample)
	i += 1

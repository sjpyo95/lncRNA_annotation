import sys, os
import config as cf
samtlPath = '/share/apps/programs/samtools/1.7/bin/samtools'

inputdir = sys.argv[1]
rpds = sys.argv[2]
server = sys.argv[3]
jobN = int(sys.argv[4])
samples = filter(lambda x: not '.t' in x and not '.out' in x and not '.log' in x, os.listdir(inputdir))
for i in xrange(len(samples)):
	sample = samples[i]
	if rpds == 'rpds':
		sinputdir = inputdir + '/' + sample + '/2.mapping/rpds_star2/'
	else:
		sinputdir = inputdir + '/'+sample+'/2.mapping/star/'
	bamfile = sinputdir + filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir))[0]
	indexCmm = samtlPath + ' index -b ' + bamfile + ' ' + bamfile + '.bai'
	queue = cf.queue(server, i)
	cf.qsub_time('index', jobN, 50)
	cf.qsub_execute('index', queue, indexCmm, sinputdir)


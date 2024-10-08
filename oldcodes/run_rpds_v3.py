import sys, os
import config as cf
rpdsPath = '/home/sjpyo/new/codes/cafe/rpds/rpds.py'
inputdir = sys.argv[1]+'/'
specie = sys.argv[2]	#uman or mouse
stranded = sys.argv[3]
outputdir = sys.argv[4]+'/'
server = sys.argv[5]
jobN = int(sys.argv[6])
if not os.path.exists(outputdir): os.makedirs(outputdir)

samples = filter(lambda x: not '.txt' in x and not '.out' in x and not '.log' in x, os.listdir(inputdir))
if specie == 'human':
	chrs = cf.hm_chrs
else:
	chrs = cf.ms_chrs
for i in xrange(len(samples)):
	sample = samples[i]
	bamdir = inputdir + sample + '/2.mapping/star/'
#	bamdir = inputdir + sample + '/'
	print bamdir
	if len(filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(bamdir))) == 0: continue
	unstranded = bamdir + filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(bamdir))[0]
	fastqdir = inputdir + sample + '/1.preprocess/sickle/'
	soutputdir = outputdir + sample + '/1.preprocess/rpds_fastqs/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	for j in xrange(len(chrs)):
		chrom = chrs[j]
		run_rpds = '/share/apps/programs/python/2.7.14/bin/python ' + rpdsPath + ' ' + unstranded + ' '+ specie + ' ' + stranded + ' ' + chrom + ' ' + fastqdir + ' ' + soutputdir
		
		queue = cf.queue(server, j)
		cf.qsub_time('rpds', jobN, 300)
		cf.qsub_execute(sample[:2]+'.rpds.'+chrom, queue, run_rpds, soutputdir+chrom+'.')


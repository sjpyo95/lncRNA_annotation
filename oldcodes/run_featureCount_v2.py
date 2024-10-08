import sys, os
import config as cf
#fcountPath = '/home/sjpyo/sjpyo/packages/subread-2.0.0-Linux-x86_64/bin/featureCounts'
fcountPath = '/export/home/sjpyo/new/packages/subread-1.6.4-source/bin/featureCounts'
#annotationDir = '/export/home/sjpyo/lncRNA_anno_profiling_exercise/hepg2/data/genome/'

#job = sys.argv[1]
inputdir = sys.argv[1] + '/'
anno = sys.argv[2]
strand = sys.argv[3] #1 or 2; 1 = fr-secondstrand, 2 = fr-firststrand
outputdir = sys.argv[4] + '/'
if not os.path.exists(outputdir): os.makedirs(outputdir)
rpds = sys.argv[5]
server = sys.argv[6]
jobN = int(sys.argv[7])

#samples = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(inputdir))
samples = filter(lambda x: not '.bam' in x and not '.t' in x and not '.out' in x and not 'total' in x and not '.log' in x, os.listdir(inputdir))
#bamfiles = []
for i in xrange(len(samples)):
	sample = samples[i]
	name = sample.split('_Aligned')[0]
#	name = sample.split('.bam')[0]
#	bamfile = inputdir + sample
	if rpds == 'rpds':
		sinputdir = inputdir + sample + '/2.mapping/rpds_star/'
	else:
		sinputdir = inputdir + sample + '/2.mapping/star/'
	bamfile = sinputdir + filter(lambda x: '.bam' in x and not '.bai' in x and name in x, os.listdir(sinputdir))[0]
#	soutputdir = inputdir + sample + '/quantification/featureCounts/'
#	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	fc_comm = fcountPath + ' -T 8 -p -s ' + strand + ' -t exon -g gene_id -M -a ' + anno + ' -o ' + outputdir + name + '.count.txt ' + bamfile
#soutputdir = outputdir + '/quantification/featureCounts/'
#if not os.path.exists(soutputdir): os.makedirs(soutputdir)
#fc_comm = fcountPath + ' -T 4 -p -t exon -g gene_id -a ' + anno + ' -o ' + soutputdir + 'total_quantity.txt ' + ' '.join(bamfiles)
	queue = cf.queue(server, i)
	cf.qsub_time('fcount', jobN, 60)
	cf.qsub_execute('fcount', queue, fc_comm, outputdir+sample+'.')
#	exit()
#	print fc_comm + '\n'
#	commands.getoutput(fc_comm)
